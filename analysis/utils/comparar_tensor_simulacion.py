#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compara el tensor de estres por grano que calcula el simulador en cada
atasco (atasco_tensor_*.dat, con la normal real de Box2D) con el que
reconstruye presion_tensor_arch.py a partir de atasco_*.xy y
atasco_contacts_*.dat (normal inferida de la geometria).

Usa las MISMAS funciones de presion_tensor_arch.py (parse_xy, parse_paredes,
fuerzas_contacto, normales_contacto, detectar_arco_primario), asi que lo que
se compara es lo que efectivamente se usa en el analisis.

Convenciones:
  simulador: s_ij = (1/A) sum_c F_i l_j  (F sobre el grano, compresion < 0)
  script:    S00 = sum r_x F_x, S11 = sum r_y F_y, S01 = sum r_x F_y,
             S10 = sum r_y F_x  ->  A s_xx = S00, A s_xy = S10,
             A s_yx = S01, A s_yy = S11.
  P = -(s_xx + s_yy)/2 en los dos (compresion > 0).

Pruebas:
  1) Contacto por contacto: la normal que infiere el script frente a la de
     Box2D (atasco_contacts_full_*.dat), por clase de contacto.
  2) Grano por grano: P, s_xx, s_yy, s_xy simetrizado y parte antisimetrica.
  3) Controles internos del tensor del simulador: equilibrio F + m g = 0,
     torque nulo, A (s_yx - s_xy) = tau y virial global.
  4) (con --r_salida) <P> pesado por area en el arco y en la ROI, por tipo
     de grano, con los dos tensores.

Uso:
  python3 analysis/utils/comparar_tensor_simulacion.py \\
      --frames /ruta/frames_<dirID> [--jam 1 2 3] [--r_salida 1.1] \\
      [--g 9.8] [--csv salida.csv]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "3_presion_y_mapas"))
from presion_tensor_arch import (parse_xy, parse_paredes, fuerzas_contacto,  # noqa: E402
                                 normales_contacto, detectar_arco_primario,
                                 leer_normales_box2d)

COLS_CONTACTOS = ['id_a', 'id_b', 'px', 'py', 'f_norm', 'f_tan']
DTYPES_CONTACTOS = {'id_a': np.int32, 'id_b': np.int32, 'px': np.float32,
                    'py': np.float32, 'f_norm': np.float32, 'f_tan': np.float32}


def leer_tabla(path):
    """Lee un .dat con cabecera '# col1 col2 ...' (la ultima linea '#')."""
    cols, dt = None, None
    with open(path) as f:
        for linea in f:
            if not linea.startswith('#'):
                break
            m = re.search(r'dt:\s*([-+0-9.eE]+)', linea)
            if m:
                dt = float(m.group(1))
            cols = linea[1:].split()
    df = pd.read_csv(path, sep=r'\s+', comment='#', names=cols,
                     na_values=['nan'])
    return df, dt


def tensor_script(xy_file, contact_file, dt, normales='auto'):
    """Tensor por grano tal como lo arma calcular_tensor_y_presion_vectorizado
    (mismas funciones y mismos tipos de datos al leer). normales: 'auto'
    (como process_one_jam_task: las de Box2D si hay atasco_contacts_full) o
    'inferidas' (siempre de la geometria, como con los datos viejos)."""
    particles = parse_xy(xy_file)
    dfp = pd.DataFrame.from_dict(particles, orient='index')
    dfp.index.name = 'id'
    dfc = pd.read_csv(contact_file, sep=r'\s+', names=COLS_CONTACTOS,
                      comment='#', dtype=DTYPES_CONTACTOS, engine='c')
    paredes = parse_paredes(xy_file)
    nb, fuente = (None, 'inferidas') if normales == 'inferidas' else \
        leer_normales_box2d(contact_file, dfc)
    df = fuerzas_contacto(dfp, dfc, dt, paredes=paredes, normales=nb)
    if nb is None:
        nx, ny = normales_contacto(df, paredes)
    else:
        nx, ny = nb[:, 0], nb[:, 1]
    df['nx_inf'], df['ny_inf'] = nx, ny
    print(f"  normales del analisis: {fuente}")

    a_ok, b_ok = df['a_valida'].to_numpy(), df['b_valida'].to_numpy()
    Fx, Fy = df['Fx_A'].to_numpy(), df['Fy_A'].to_numpy()
    rax, ray = (df['px'] - df['pos_x_A']).to_numpy(), (df['py'] - df['pos_y_A']).to_numpy()
    rbx, rby = (df['px'] - df['pos_x_B']).to_numpy(), (df['py'] - df['pos_y_B']).to_numpy()
    comps = {'S00': (rax, rbx, Fx), 'S11': (ray, rby, Fy),
             'S01': (rax, rbx, Fy), 'S10': (ray, rby, Fx)}
    tA = pd.DataFrame({k: np.where(a_ok, ra * F, 0.0) for k, (ra, rb, F) in comps.items()})
    tB = pd.DataFrame({k: np.where(b_ok, -rb * F, 0.0) for k, (ra, rb, F) in comps.items()})
    tA['id'], tB['id'] = df['id_a'].to_numpy(), df['id_b'].to_numpy()
    S = pd.concat([tA, tB]).groupby('id')[list(comps)].sum()
    dfp = dfp.join(S).fillna({k: 0.0 for k in comps})
    A = dfp['area'].to_numpy(np.float64)
    out = pd.DataFrame(index=dfp.index)
    out['A_py'] = A
    out['sxx_py'] = dfp['S00'] / A
    out['sxy_py'] = dfp['S10'] / A
    out['syx_py'] = dfp['S01'] / A
    out['syy_py'] = dfp['S11'] / A
    out['shape'] = dfp['shape']
    return out, df, dfp, dfc


def clase_contacto(df):
    sa, sb = df['shape_A'], df['shape_B']
    def tipo(s):
        return np.where(s.isna(), 'pared', np.where(s == 'Disco', 'disco', 'poligono'))
    ta, tb = tipo(sa), tipo(sb)
    pares = np.array([' - '.join(sorted(p)) for p in zip(ta, tb)])
    return pares


def resumen_normales(df, full):
    """Angulo entre la normal inferida y la de Box2D, por clase de contacto."""
    if len(df) != len(full):
        print(f"  ! filas distintas: contactos {len(df)}, contacts_full {len(full)}")
        return
    ok_ids = ((df['id_a'].to_numpy() == full['ID_A'].to_numpy()) &
              (df['id_b'].to_numpy() == full['ID_B'].to_numpy()))
    if not ok_ids.all():
        print(f"  ! {np.sum(~ok_ids)} filas con IDs distintos; no se comparan normales")
        return
    cosang = np.clip(df['nx_inf'] * full['nx'] + df['ny_inf'] * full['ny'], -1, 1)
    ang = np.degrees(np.arccos(cosang))
    fn = full['Fn'].to_numpy()
    clase = clase_contacto(df)
    tabla = pd.DataFrame({'clase': clase, 'npc': full['npc'], 'ang': ang, 'Fn': fn})
    print("  Normal del analisis vs Box2D (angulo, grados; solo puntos con Fn > 0):")
    print(f"    {'clase':22s} {'npc':>3s} {'N':>6s} {'>1 deg':>8s} {'>10 deg':>8s}"
          f" {'max':>8s} {'%Fn mal':>8s}")
    act = tabla[tabla['Fn'] > 0]
    for (c, npc), g in act.groupby(['clase', 'npc']):
        mal = g['ang'] > 1.0
        frac_fn = 100 * g.loc[mal, 'Fn'].sum() / g['Fn'].sum() if g['Fn'].sum() > 0 else 0
        print(f"    {c:22s} {npc:3d} {len(g):6d} {mal.sum():8d} {(g['ang'] > 10).sum():8d}"
              f" {g['ang'].max():8.2f} {frac_fn:8.2f}")
    return tabla


def controles_internos(ts, full, g):
    """Equilibrio, torque e identidad del virial con el tensor del simulador."""
    m, A = ts['m'].to_numpy(), ts['A'].to_numpy()
    w = m * abs(g)
    res_F = np.hypot(ts['Fx'], ts['Fy'] - w) / w
    R = np.sqrt(A / np.pi)
    res_tau = np.abs(ts['tau']) / (w * R)
    ident = np.abs(A * (ts['syx'] - ts['sxy']) - ts['tau']) / (w * R)
    print("  Controles del tensor del simulador:")
    print(f"    |F + m g| / (m g):      mediana {np.median(res_F):.2e}  "
          f"p95 {np.percentile(res_F, 95):.2e}  max {res_F.max():.2e}")
    print(f"    |tau| / (m g R):        mediana {np.median(res_tau):.2e}  "
          f"p95 {np.percentile(res_tau, 95):.2e}  max {res_tau.max():.2e}")
    print(f"    A(syx-sxy) - tau:       max {ident.max():.2e} (debe ser ~1e-7)")
    if 'vx' in ts:
        v = np.hypot(ts['vx'], ts['vy']) / np.sqrt(abs(g) * 2 * R)
        print(f"    |v| / sqrt(g d):        mediana {np.median(v):.2e}  "
              f"p95 {np.percentile(v, 95):.2e}  max {v.max():.2e}")

    # Virial: sum_p A s_p + sum_p F_p (x) x_p = sum_ext f (x) x_c, con "ext"
    # los puntos de contacto entre un grano activo y algo que no lo es.
    act = set(ts['gID'])
    a_in = full['ID_A'].isin(act).to_numpy()
    b_in = full['ID_B'].isin(act).to_numpy()
    nx, ny = full['nx'].to_numpy(), full['ny'].to_numpy()
    FBx = full['Fn'] * nx + full['Ft'] * ny
    FBy = full['Fn'] * ny - full['Ft'] * nx
    px, py = full['PosX'].to_numpy(), full['PosY'].to_numpy()
    ext = np.zeros((2, 2))
    for sel, sgn in ((a_in & ~b_in, -1.0), (b_in & ~a_in, 1.0)):
        fx, fy = sgn * FBx[sel], sgn * FBy[sel]
        ext += np.array([[np.sum(fx * px[sel]), np.sum(fx * py[sel])],
                         [np.sum(fy * px[sel]), np.sum(fy * py[sel])]])
    x, y = ts['x'].to_numpy(), ts['y'].to_numpy()
    lhs = np.array([[np.sum(A * ts['sxx']) + np.sum(ts['Fx'] * x),
                     np.sum(A * ts['sxy']) + np.sum(ts['Fx'] * y)],
                    [np.sum(A * ts['syx']) + np.sum(ts['Fy'] * x),
                     np.sum(A * ts['syy']) + np.sum(ts['Fy'] * y)]])
    escala = np.abs(np.sum(A * (ts['sxx'] + ts['syy'])))
    print(f"    virial global: |izq - der| / |sum A tr s| = "
          f"{np.abs(lhs - ext).max() / escala:.2e}")


def comparar_granos(ts, tp):
    """Diferencias grano por grano. La escala de cada grano es |P_sim| +
    |tau_dev_sim| (para no dividir por ~0 en granos casi sin carga)."""
    d = ts.set_index('gID').join(tp, how='inner')
    P_s = -0.5 * (d['sxx'] + d['syy'])
    P_p = -0.5 * (d['sxx_py'] + d['syy_py'])
    q_s = np.hypot(0.5 * (d['sxx'] - d['syy']), 0.5 * (d['sxy'] + d['syx']))
    esc = np.abs(P_s) + q_s
    cargados = esc > 1e-3 * np.median(esc[esc > 0])
    errs = {
        'P': (P_p - P_s),
        'sxx': d['sxx_py'] - d['sxx'],
        'syy': d['syy_py'] - d['syy'],
        'sxy_sim': 0.5 * (d['sxy_py'] + d['syx_py']) - 0.5 * (d['sxy'] + d['syx']),
        'antisim': 0.5 * (d['syx_py'] - d['sxy_py']) - 0.5 * (d['syx'] - d['sxy']),
    }
    print(f"  Grano por grano ({cargados.sum()} granos con carga de {len(d)}):"
          " error relativo a |P| + q del simulador")
    print(f"    {'magnitud':10s} {'mediana':>10s} {'p95':>10s} {'max':>10s} {'>1%':>6s} {'>10%':>6s}")
    for k, e in errs.items():
        r = (np.abs(e) / esc)[cargados]
        print(f"    {k:10s} {np.median(r):10.2e} {np.percentile(r, 95):10.2e} {r.max():10.2e}"
              f" {(r > 0.01).sum():6d} {(r > 0.1).sum():6d}")
    dA = np.abs(d['A_py'] - d['A']) / d['A']
    print(f"    area (xy vs Box2D): max error relativo {dA.max():.2e}")
    d['P_sim'], d['P_py'] = P_s, P_p
    return d


def presion_por_region(d, dfp, dfc, r_salida):
    """<P> pesado por area en arco y ROI, por tipo, con los dos tensores
    (la misma definicion de ROI que presion_tensor_arch.py)."""
    arco = detectar_arco_primario(dfp, dfc, r_salida)
    if not arco:
        print("  (no se detecto arco)")
        return
    y_arco = d.loc[d.index.isin(arco), 'y']
    ymin, ymax = y_arco.min(), y_arco.max()
    region = np.where(d.index.isin(arco), 'arco',
                      np.where((d['y'] >= ymin) & (d['y'] <= ymax), 'roi', 'bulk'))
    d = d.assign(region=region)
    print("  <P> pesado por area (simulador / script):")
    for (reg, tipo), g in d[d['region'] != 'bulk'].groupby(['region', 'tipo']):
        ps = np.sum(g['P_sim'] * g['A']) / g['A'].sum()
        pp = np.sum(g['P_py'] * g['A_py']) / g['A_py'].sum()
        print(f"    {reg:5s} tipo {tipo}: N={len(g):4d}  {ps:12.5g} / {pp:12.5g}"
              f"  ({100 * (pp - ps) / ps:+.3f} %)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--frames', required=True, help='carpeta frames_<dirID>')
    ap.add_argument('--jam', nargs='*', type=int, help='atascos (default: todos)')
    ap.add_argument('--g', type=float, default=9.8, help='|g| (default 9.8)')
    ap.add_argument('--r_salida', type=float, default=None,
                    help='radio del orificio, para detectar el arco')
    ap.add_argument('--csv', default=None, help='guarda la comparacion por grano')
    ap.add_argument('--normales', choices=['auto', 'inferidas'], default='auto',
                    help="auto: como presion_tensor_arch.py (Box2D si hay "
                         "atasco_contacts_full); inferidas: siempre de la geometria")
    args = ap.parse_args()

    if args.jam:
        jams = args.jam
    else:
        jams = sorted(int(re.search(r'(\d+)\.dat$', f).group(1)) for f in
                      glob.glob(os.path.join(args.frames, 'atasco_tensor_*.dat')))
    filas = []
    for j in jams:
        tag = f"{j:05d}"
        xy = os.path.join(args.frames, f"atasco_{tag}.xy")
        cf = os.path.join(args.frames, f"atasco_contacts_{tag}.dat")
        ff = os.path.join(args.frames, f"atasco_contacts_full_{tag}.dat")
        tf = os.path.join(args.frames, f"atasco_tensor_{tag}.dat")
        print(f"=== atasco {tag}")
        ts, dt = leer_tabla(tf)
        full, _ = leer_tabla(ff)
        tp, df, dfp, dfc = tensor_script(xy, cf, dt, args.normales)
        resumen_normales(df, full)
        controles_internos(ts, full, args.g)
        d = comparar_granos(ts, tp)
        if args.r_salida is not None:
            presion_por_region(d, dfp, dfc, args.r_salida)
        d['atasco'] = j
        filas.append(d)
    if args.csv and filas:
        pd.concat(filas).to_csv(args.csv)
        print(f"Comparacion por grano guardada en {args.csv}")


if __name__ == '__main__':
    main()
