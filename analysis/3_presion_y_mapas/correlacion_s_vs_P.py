#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correlacion entre el tamano medio de avalancha <s> y la presion media <P>
de los atascos, un punto por simulacion (sistema, chi, razon de areas,
orificio).

Figura: un solo panel, color y marcador por sistema (estilo_figuras.py),
puntos sin lineas. Se usa un solo orificio W (el _W_ de la carpeta, --W,
default 2.0) para que <s> sea comparable. --logy / --logx / --loglog: eje
y, eje x o los dos en escala logaritmica (el pdf lleva el sufijo -logy,
-logx o -loglog, asi no pisa la version lineal).
--paneles: un panel por sistema con el color segun chi.
--s_en_x: <s> en el eje x (variable independiente) y la presion en el eje y.
Los --logx / --logy siguen refiriendose a la presion y a <s> respectivamente
(--logy = <s> logaritmico, ahora en el eje x); el pdf se llama
16-P-vs-s-... con sufijo -logs, -logP o -loglog.

Entradas:
  resumen_presion_por_atasco.csv     (presion_tensor_arch.py)
  resultados_completos_esqueleto.csv (processor.py)
  master_sim_index.csv               (indexer.py: sim_hash <-> carpeta)
La simulacion se empareja por las dos ultimas carpetas de la ruta
(BARRIDO_.../Chi_..._W_...), asi no importa si el indice y la presion se
calcularon con el disco montado en rutas distintas.

Salida (en processed_data/5_graficos_finales/ por defecto):
  16-s-vs-P-<region>-ratio-<r>-W-<W>.pdf  y  ..._datos.csv con los valores graficados.
En la leyenda va la razon de areas R_a; el coeficiente de Pearson r de cada
sistema solo se imprime por pantalla (con --logy, se calcula con log<s>).

Uso:
    python3 analysis/3_presion_y_mapas/correlacion_s_vs_P.py \\
        --presion processed_data/presion/tensor_nuevo/resumen_presion_por_atasco.csv
    python3 analysis/3_presion_y_mapas/correlacion_s_vs_P.py --region roi --W 2.0 --radio_salida 1.0
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # rutas.py y estilo_figuras.py
import rutas  # noqa: E402
from estilo_figuras import (aplicar_estilo, estilo_sistema, guardar_figura,  # noqa: E402
                            ERRORBAR_KW, SYSTEM_NAME_MAP)
from comparar_presion_por_chi import filtrar_lotes_confirmados, _estadistica_pesada  # noqa: E402

COLUMNA_PRESION = {'P': 'SumP_reducida', 'tau': 'SumTau_reducida'}
# Presion reducida con la escala mg/d (fuerza por unidad de longitud en 2D),
# m = masa de un grano grande, d = diametro de un disco grande; en las
# simulaciones m = d = 1, asi que el valor numerico es P/g.
# Sufijo para no pisar las figuras con la etiqueta vieja (Mg).
SUFIJO_NUEVO = '_nuevo'
ETIQUETA_X = {'P': r'Reduced Mean Pressure $\langle P\rangle d / mg$',
              'tau': r'Reduced Mean Shear Stress $\langle \tau\rangle d / mg$'}
ETIQUETA_Y = {'big': r'Mean Avalanche Size $\langle s \rangle$',
              'total': r'Mean Avalanche Size $\langle s \rangle$ (all grains)'}
COLUMNA_S = {'big': 'tamano_big', 'total': 'total_mix'}
NOMBRE_REGION = {'arco': 'arch', 'roi': 'ROI', 'bulk': 'bulk'}


W_REGEX = re.compile(r'_W_([0-9]+(?:p[0-9]+)?)')
MARCADORES = ['o', 's', '^', 'D', 'v', 'p', 'h', '<', '>']


def ancho_W(clave):
    m = W_REGEX.search(str(clave))
    return float(m.group(1).replace('p', '.')) if m else np.nan


def clave_carpeta(ruta):
    """Las dos ultimas carpetas de una ruta: 'BARRIDO_X/Chi_0p2_W_2p0'."""
    partes = os.path.normpath(str(ruta)).replace('\\', '/').rstrip('/').split('/')
    return '/'.join(partes[-2:])


def presion_por_simulacion(path_csv, magnitud, region):
    df = pd.read_csv(path_csv, dtype={'lote': str, 'jam_id': str})
    df = filtrar_lotes_confirmados(df, path_csv)
    col_sum, col_area = f"{COLUMNA_PRESION[magnitud]}_{region}", f"Area_{region}"
    faltan = [c for c in (col_sum, col_area) if c not in df.columns]
    if faltan:
        sys.exit(f"[!] Al csv de presion le faltan las columnas {faltan}")
    filas = []
    for sim_path, g in df.groupby('sim_path'):
        w = g[col_area].to_numpy(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            x = g[col_sum].to_numpy(float) / w
        media, std, sem, n, _ = _estadistica_pesada(x, w)
        filas.append(dict(clave=clave_carpeta(sim_path), sim_path=sim_path, P=media, P_sem=sem,
                          n_atascos=n))
    return pd.DataFrame(filas)


def avalanchas_por_simulacion(path_av, path_idx, col_s, solo_con_arco):
    av = pd.read_csv(path_av)
    if solo_con_arco:
        if 'atasco_id' not in av.columns:
            sys.exit("[!] El esqueleto no tiene atasco_id: volve a correr processor.py (version nueva).")
        av = av[av['atasco_id'].notna()]
    idx = pd.read_csv(path_idx)
    if 'area_ratio_nominal' not in idx.columns and 'area_ratio' in idx.columns:
        idx['area_ratio_nominal'] = idx['area_ratio'].round(3)
    idx = idx.drop_duplicates(subset=['sim_hash'])
    idx['clave'] = idx['path'].map(clave_carpeta)
    g = av.groupby('sim_hash')[col_s].agg(['mean', 'sem', 'count']).reset_index()
    g.columns = ['sim_hash', 's', 's_sem', 'n_avalanchas']
    cols = [c for c in ('sim_hash', 'clave', 'system_type', 'chi', 'area_ratio_nominal', 'Radio_salida')
            if c in idx.columns]
    return g.merge(idx[cols], on='sim_hash', how='inner')


def agregar_puros(df):
    """Las simulaciones puras (chi = 0, system_type sin '-', ej. 'Cuad') se
    suman como punto chi = 0 de cada mezcla con el mismo grano grande
    ('Cuad-Disco'), igual que la referencia pura de graficos_finales.py."""
    st = df['system_type'].astype(str)
    puros, mezclas = df[~st.str.contains('-')], df[st.str.contains('-')]
    extra = []
    for sis in mezclas['system_type'].unique():
        base = str(sis).split('-')[0]
        p = puros[puros['system_type'] == base]
        if len(p) and not ((mezclas['system_type'] == sis) & (mezclas['chi'] == 0)).any():
            extra.append(p.assign(system_type=sis))
    return pd.concat([mezclas] + extra, ignore_index=True) if len(mezclas) else df


def ejes_log(ax, logx, logy):
    """Escala log con numeros comunes en las marcas (4, 10, 100 en vez de
    4x10^0). Las marcas se eligen segun el rango: si cubre menos de ~1.5
    decadas, 1..9 por decada (ej. 50, 60, 70); si no, 1, 2, 5."""
    from matplotlib.ticker import ScalarFormatter, NullFormatter, LogLocator
    for usar, eje, set_scale, get_lim in ((logx, ax.xaxis, ax.set_xscale, ax.get_xlim),
                                          (logy, ax.yaxis, ax.set_yscale, ax.get_ylim)):
        if not usar:
            continue
        set_scale('log')
        lo, hi = get_lim()
        decadas = np.log10(hi / lo) if lo > 0 and hi > lo else 2.0
        subs = np.arange(1, 10) if decadas < 1.5 else (1.0, 2.0, 5.0)
        eje.set_major_locator(LogLocator(base=10, subs=subs))
        fmt = ScalarFormatter()
        fmt.set_scientific(False)
        eje.set_major_formatter(fmt)
        eje.set_minor_formatter(NullFormatter())


def pearson(x, y, logy):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if logy:
        ok = y > 0
        x, y = x[ok], np.log(y[ok])
    return np.corrcoef(x, y)[0, 1] if len(x) > 2 else np.nan


def ejes_xy(magnitud, region, tipo_s, s_en_x):
    """(columna, error, etiqueta) del eje x y del eje y."""
    pres = ('P', 'P_sem', f"{ETIQUETA_X[magnitud]} ({NOMBRE_REGION[region]})")
    aval = ('s', 's_sem', ETIQUETA_Y[tipo_s])
    return (aval, pres) if s_en_x else (pres, aval)


def graficar_paneles(df, out, magnitud, region, tipo_s, logx, logy, error, chi_max, area_ratio,
                     s_en_x=False):
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    aplicar_estilo()
    (xc, xe, xl), (yc, ye, yl) = ejes_xy(magnitud, region, tipo_s, s_en_x)
    if s_en_x:
        logx, logy = logy, logx
    sistemas = sorted(df['sistema'].unique(), key=lambda s: (s != 'Disco(G)-Disco(C)', s))
    ncol = 2 if len(sistemas) > 1 else 1
    nfil = int(np.ceil(len(sistemas) / ncol))
    fig, axes = plt.subplots(nfil, ncol, figsize=(6.2 * ncol, 4.6 * nfil), sharex=True, sharey=True,
                             squeeze=False)
    norm, cmap = Normalize(0, chi_max), plt.cm.viridis
    familias = sorted(df[['Radio_salida', 'W']].drop_duplicates().itertuples(index=False, name=None),
                      key=lambda t: (np.nan_to_num(t[0]), np.nan_to_num(t[1])))
    marcador = {f: MARCADORES[i % len(MARCADORES)] for i, f in enumerate(familias)}
    vistos = set()
    for ax, sis in zip(axes.flat, sistemas):
        sub = df[df['sistema'] == sis]
        for fam, g in sub.groupby(['Radio_salida', 'W'], dropna=False):
            fam = tuple(fam)
            g = g.sort_values('chi')
            m = marcador.get(fam, 'o')
            if error:
                ax.errorbar(g[xc], g[yc], xerr=g[xe], yerr=g[ye], fmt='none', ecolor='0.55',
                            elinewidth=0.8, capsize=2, zorder=2)
            lbl = None
            if len(familias) > 1 and fam not in vistos:
                lbl = f"$R_{{out}}$ = {fam[0]:g}, W = {fam[1]:g}"
                vistos.add(fam)
            ax.scatter(g[xc], g[yc], c=g['chi'], cmap=cmap, norm=norm, marker=m, s=55,
                       edgecolors='black', linewidths=0.6, zorder=3, label=lbl)
        ax.set_title(rf"{estilo_sistema(sis)[2]}  ($R_{{\rm a}} = {area_ratio:g}$)", fontsize='medium')
        ejes_log(ax, logx, logy)
    for ax in axes.flat[len(sistemas):]:
        ax.set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel(xl)
    for ax in axes[:, 0]:
        ax.set_ylabel(yl)
    if len(familias) > 1:
        fig.legend(loc='lower center', ncol=min(len(familias), 4), frameon=False, fontsize='small',
                   bbox_to_anchor=(0.45, -0.02))
    fig.tight_layout(rect=(0, 0.05 if len(familias) > 1 else 0, 0.9, 1))
    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, label=r'$\chi$')
    return guardar_figura(out, fig)


def graficar(df, out, magnitud, region, tipo_s, logx, logy, error, area_ratio, s_en_x=False):
    aplicar_estilo()
    (xc, xe, xl), (yc, ye, yl) = ejes_xy(magnitud, region, tipo_s, s_en_x)
    if s_en_x:
        logx, logy = logy, logx
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    sistemas = sorted(df['sistema'].unique(), key=lambda s: (s != 'Disco(G)-Disco(C)', s))
    for sis in sistemas:
        sub = df[df['sistema'] == sis].sort_values('chi')
        color, marker, etiqueta = estilo_sistema(sis)
        lbl = etiqueta
        xerr = sub[xe] if error else None
        yerr = sub[ye] if error else None
        ax.errorbar(sub[xc], sub[yc], xerr=xerr, yerr=yerr, color=color, marker=marker, ls='none',
                    ms=7, label=lbl, **ERRORBAR_KW)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ejes_log(ax, logx, logy)
    ax.legend(title=rf"$R_{{\rm a}} = {area_ratio:g}$", frameon=False,
              loc='best', fontsize='small')
    return guardar_figura(out, fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--presion', default=os.path.join(rutas.PRESION_DIR, 'tensor_nuevo',
                                                      'resumen_presion_por_atasco.csv'),
                    help='resumen_presion_por_atasco.csv de presion_tensor_arch.py')
    ap.add_argument('--av', default=os.path.join(rutas.AVALANCHE_DIR, 'resultados_completos_esqueleto.csv'))
    ap.add_argument('--idx', default=rutas.MASTER_INDEX)
    ap.add_argument('--region', choices=['arco', 'roi', 'bulk'], default='arco')
    ap.add_argument('--magnitud', choices=['P', 'tau'], default='P', help='P (presion) o tau (deviador)')
    ap.add_argument('--s', choices=['big', 'total'], default='big',
                    help="big: granos grandes por avalancha (tamano_big, default); total: total_mix")
    ap.add_argument('--area_ratio', type=float, default=0.1, help='razon de areas nominal (default 0.1)')
    ap.add_argument('--radio_salida', type=float, default=None, help='filtra Radio_salida; por defecto todos')
    ap.add_argument('--chi_max', type=float, default=0.4)
    ap.add_argument('--sistemas', nargs='+', default=None,
                    help="system_type del indice (ej. Disco-Disco Cuad-Disco); por defecto todos")
    ap.add_argument('--solo_con_arco', action='store_true',
                    help='solo avalanchas que terminaron en un atasco guardado (con arco)')
    ap.add_argument('--W', type=float, default=2.0,
                    help='orificio W de la carpeta (_W_2p0 -> 2.0, default); un solo orificio por figura')
    ap.add_argument('--logx', action='store_true', help='eje x (presion) logaritmico')
    ap.add_argument('--logy', action='store_true', help='eje y (<s>) logaritmico (r se calcula con log<s>)')
    ap.add_argument('--loglog', action='store_true', help='los dos ejes logaritmicos (= --logx --logy)')
    ap.add_argument('--paneles', action='store_true',
                    help='un panel por sistema con el color segun chi (por defecto: un solo panel, '
                         'color y marcador por sistema)')
    ap.add_argument('--s_en_x', action='store_true',
                    help='<s> en el eje x (variable independiente) y la presion en el eje y; '
                         '--logy sigue siendo <s> logaritmico')
    ap.add_argument('--sin_error', action='store_true', help='sin barras de error')
    ap.add_argument('--out', default=None,
                    help='pdf de salida (default: 5_graficos_finales/16-s-vs-<P|tau>-<region>-ratio-<r>-W-<W>.pdf)')
    args = ap.parse_args()
    if args.loglog:
        args.logx = args.logy = True

    pres = presion_por_simulacion(args.presion, args.magnitud, args.region)
    aval = avalanchas_por_simulacion(args.av, args.idx, COLUMNA_S[args.s], args.solo_con_arco)
    df = aval.merge(pres, on='clave', how='inner')
    print(f"Simulaciones: {len(pres)} con presion, {len(aval)} con avalanchas, {len(df)} en ambas.")
    if df.empty:
        print("  ejemplos de claves (presion):  ", pres['clave'].head(3).tolist())
        print("  ejemplos de claves (avalanchas):", aval['clave'].head(3).tolist())
        sys.exit("[!] Ninguna simulacion coincide entre la presion y las avalanchas.")
    if 'area_ratio_nominal' in df.columns and args.area_ratio is not None:
        df = df[np.isclose(df['area_ratio_nominal'].fillna(1.0), args.area_ratio, atol=1e-3)
                | (df['chi'] == 0)]
    df['W'] = df['clave'].map(ancho_W)
    if args.radio_salida is not None:
        df = df[np.isclose(df['Radio_salida'], args.radio_salida)]
    df = df[np.isclose(df['W'], args.W)]
    df = df[df['chi'] <= args.chi_max]
    if args.sistemas:
        df = df[df['system_type'].isin(args.sistemas)]
    df = df.dropna(subset=['P', 's'])
    df = agregar_puros(df)
    if df.empty:
        sys.exit("[!] No quedo ninguna simulacion despues de filtrar (--W, --area_ratio, --radio_salida, --chi_max).")
    df['sistema'] = df['system_type'].map(lambda s: SYSTEM_NAME_MAP.get(str(s), str(s)))

    if args.s_en_x:
        nombre = f"16-{args.magnitud}-vs-s"
        sufijo = ('-loglog' if args.logx and args.logy else f'-log{args.magnitud}' if args.logx
                  else '-logs' if args.logy else '')
    else:
        nombre = f"16-s-vs-{args.magnitud}"
        sufijo = '-loglog' if args.logx and args.logy else '-logx' if args.logx else '-logy' if args.logy else ''
    out = args.out or os.path.join(
        rutas.GRAFICOS_FINALES_DIR,
        f"{nombre}-{args.region}-ratio-{args.area_ratio:g}-W-{args.W:g}{sufijo}{SUFIJO_NUEVO}.pdf")
    datos = os.path.splitext(out)[0] + '_datos.csv'
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    cols = ['sistema', 'chi', 'area_ratio_nominal', 'Radio_salida', 'W', 'clave', 'P', 'P_sem', 'n_atascos',
            's', 's_sem', 'n_avalanchas']
    df.sort_values(['sistema', 'chi'])[[c for c in cols if c in df.columns]].to_csv(datos, index=False)
    print("Puntos por sistema y familia (Radio_salida, W): valores de chi")
    for (sis, rs, w), g in df.groupby(['sistema', 'Radio_salida', 'W'], dropna=False):
        print(f"  {sis:20s} R_out={rs:<5g} W={w:<5g} chi={sorted(g['chi'].round(3).tolist())}")
    print("Pearson r (solo informativo, no va en la figura):")
    for sis, g in df.groupby('sistema'):
        print(f"  {sis:20s} r = {pearson(g['P'], g['s'], args.logy):+.2f}")
    if args.paneles:
        path = graficar_paneles(df, out, args.magnitud, args.region, args.s, args.logx, args.logy,
                                not args.sin_error, args.chi_max, args.area_ratio, args.s_en_x)
    else:
        path = graficar(df, out, args.magnitud, args.region, args.s, args.logx, args.logy,
                        not args.sin_error, args.area_ratio, args.s_en_x)
    print(f"Figura: {path}\nDatos:  {datos}")


if __name__ == '__main__':
    main()
