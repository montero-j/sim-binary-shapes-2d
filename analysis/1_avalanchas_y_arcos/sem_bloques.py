#!/usr/bin/env python3
"""
Error estandar de <s> por bloques (Flyvbjerg-Petersen) y autocorrelacion de
avalanchas sucesivas, para ver si el SEM ingenuo (std/sqrt(n), el que usa
flow_analysis.py) subestima el error por correlaciones dentro de una corrida.

Para cada simulacion se ordenan sus avalanchas por t_inicio y se agrupan de a
b = 1, 2, 4, 8, ... consecutivas. SEM(b) = std de las medias de bloque /
sqrt(n_bloques). Si las avalanchas son independientes, SEM(b) = SEM(1) para
todo b (dentro del ruido); si estan correlacionadas, SEM(b) crece hasta un
plateau, que es el error correcto. Se usan solo los b con al menos
--min_bloques bloques.

Ademas: autocorrelacion rho(k) de s a distancia k avalanchas y tiempo de
autocorrelacion integrado tau_int = 1/2 + sum_k rho(k) (hasta el primer
rho <= 0); con avalanchas independientes tau_int = 1/2 y el factor
sqrt(2 tau_int) que corrige el SEM es 1.

Salidas (en 5_graficos_finales/ por defecto):
  sem_bloques_<etiqueta>.csv   una fila por simulacion: n, SEM ingenuo,
                               SEM de bloques (plateau), cociente, rho(1), tau_int
  sem_bloques_<etiqueta>.pdf   SEM(b)/SEM(1) en funcion de b, una curva por chi

Uso (figura del paper, disco-disco R_a = 0.1):
  python3 sem_bloques.py --sistema Disco-Disco --area_ratio 0.1
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import rutas  # noqa: E402
from estilo_figuras import aplicar_estilo, guardar_figura, ERRORBAR_KW, REF_LINE_KW, SEQ_STYLE  # noqa: E402


def sem_por_bloques(x, min_bloques=16):
    """Lista de (b, n_bloques, SEM(b), error de SEM(b)) para b = 1, 2, 4, ...
    El error de SEM(b) es SEM(b) / sqrt(2 (n_bloques - 1))."""
    x = np.asarray(x, dtype=float)
    out, b = [], 1
    while len(x) // b >= min_bloques:
        nb = len(x) // b
        medias = x[:nb * b].reshape(nb, b).mean(axis=1)
        sem = medias.std(ddof=1) / np.sqrt(nb)
        out.append((b, nb, sem, sem / np.sqrt(2 * (nb - 1))))
        b *= 2
    return out


def autocorrelacion(x, k_max=50):
    """rho(k) para k = 1..k_max (estimador normalizado por la varianza)."""
    x = np.asarray(x, dtype=float) - np.mean(x)
    var = np.mean(x * x)
    if var == 0:
        return np.zeros(0)
    k_max = min(k_max, len(x) // 4)
    return np.array([np.mean(x[:-k] * x[k:]) / var for k in range(1, k_max + 1)])


def tau_integrado(rho):
    """1/2 + suma de rho(k) hasta el primer rho <= 0 (corte de ventana simple)."""
    tau = 0.5
    for r in rho:
        if r <= 0:
            break
        tau += r
    return tau


def analizar(df_av, min_bloques):
    filas, curvas = [], {}
    for h, g in df_av.groupby('sim_hash'):
        x = g.sort_values('t_inicio', kind='stable')['tamano_big'].to_numpy(float)
        if len(x) < 2 * min_bloques:
            continue
        bl = sem_por_bloques(x, min_bloques)
        sem1 = bl[0][2]
        # plateau: SEM del b mas grande con suficientes bloques (cota conservadora)
        sem_pl, err_pl = max(bl, key=lambda t: t[0])[2:4]
        rho = autocorrelacion(x)
        tau = tau_integrado(rho)
        meta = g.iloc[0]
        filas.append(dict(sim_hash=h, system_type=meta.get('system_type'), chi=meta.get('chi'),
                          area_ratio_nominal=meta.get('area_ratio_nominal'), n_avalanchas=len(x),
                          s_medio=x.mean(), sem_ingenuo=sem1, sem_bloques=sem_pl,
                          err_sem_bloques=err_pl, cociente=sem_pl / sem1,
                          b_max=max(t[0] for t in bl),
                          rho1=rho[0] if len(rho) else np.nan, tau_int=tau,
                          sem_tau=sem1 * np.sqrt(2 * tau)))
        curvas[h] = bl
    return pd.DataFrame(filas), curvas


def graficar(res, curvas, out):
    aplicar_estilo()
    fig, ax = plt.subplots()
    res = res.sort_values('chi')
    colores, marcadores = SEQ_STYLE(len(res))
    for (_, f), c, m in zip(res.iterrows(), colores, marcadores):
        bl = curvas[f['sim_hash']]
        b = np.array([t[0] for t in bl])
        sem1 = bl[0][2]
        ax.errorbar(b, [t[2] / sem1 for t in bl], yerr=[t[3] / sem1 for t in bl],
                    fmt=f'{m}-', color=c, label=rf"$\chi = {f['chi']:g}$", **ERRORBAR_KW)
    ax.axhline(1.0, **REF_LINE_KW)
    ax.set_xscale('log', base=2)
    ax.set_xlabel(r"Block size $b$ (consecutive avalanches)")
    ax.set_ylabel(r"$\mathrm{SEM}(b)/\mathrm{SEM}(1)$")
    ax.legend(loc='best', fontsize='small', ncol=2)
    return guardar_figura(out, fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--av', default=os.path.join(rutas.AVALANCHE_DIR, 'resultados_completos_esqueleto.csv'),
                    help='avalanchas una por fila (processor.py): sim_hash, t_inicio, tamano_big')
    ap.add_argument('--index', default=rutas.MASTER_INDEX)
    ap.add_argument('--sistema', nargs='+', default=None, help='system_type (ej. Disco-Disco Cuad-Disco)')
    ap.add_argument('--area_ratio', type=float, nargs='+', default=None)
    ap.add_argument('--chi', type=float, nargs='+', default=None)
    ap.add_argument('--min_bloques', type=int, default=16,
                    help='minimo de bloques para usar un tamaño de bloque (default 16)')
    ap.add_argument('--etiqueta', default=None, help='sufijo de los archivos de salida')
    ap.add_argument('--out_dir', default=rutas.GRAFICOS_FINALES_DIR)
    args = ap.parse_args()

    av = pd.read_csv(args.av)
    faltan = [c for c in ('sim_hash', 'tamano_big') if c not in av.columns]
    if faltan:
        sys.exit(f"[!] A {args.av} le faltan las columnas {faltan}")
    if 't_inicio' not in av.columns:
        # csv de una version anterior de processor.py: las avalanchas de cada
        # simulacion se escriben en el orden en que ocurren (rachas de flx.dat),
        # asi que el orden de las filas dentro de cada sim_hash es el temporal.
        print("[aviso] el csv no tiene t_inicio: se usa el orden de las filas dentro de cada simulacion.")
        av['t_inicio'] = av.groupby('sim_hash').cumcount().astype(float)
    idx = pd.read_csv(args.index)
    meta = [c for c in ('sim_hash', 'system_type', 'chi', 'area_ratio_nominal') if c in idx.columns]
    av = av.drop(columns=[c for c in meta if c != 'sim_hash' and c in av.columns]).merge(
        idx[meta].drop_duplicates('sim_hash'), on='sim_hash', how='left')
    if args.sistema:
        av = av[av['system_type'].isin(args.sistema)]
    if args.area_ratio:
        av = av[np.isclose(av['area_ratio_nominal'].to_numpy(float)[:, None], args.area_ratio).any(axis=1)]
    if args.chi:
        av = av[np.isclose(av['chi'].to_numpy(float)[:, None], args.chi).any(axis=1)]
    if av.empty:
        sys.exit("[!] No quedaron avalanchas despues de filtrar (--sistema / --area_ratio / --chi).")

    res, curvas = analizar(av, args.min_bloques)
    if res.empty:
        sys.exit(f"[!] Ninguna simulacion tiene al menos {2 * args.min_bloques} avalanchas.")
    etiqueta = args.etiqueta or '_'.join(filter(None, [
        '-'.join(args.sistema) if args.sistema else 'todos',
        f"Ra{'-'.join(f'{r:g}' for r in args.area_ratio)}" if args.area_ratio else None]))
    os.makedirs(args.out_dir, exist_ok=True)
    csv = os.path.join(args.out_dir, f"sem_bloques_{etiqueta}.csv")
    res.sort_values(['system_type', 'chi']).to_csv(csv, index=False)

    cols = ['sim_hash', 'system_type', 'chi', 'n_avalanchas', 's_medio', 'sem_ingenuo', 'sem_bloques',
            'cociente', 'b_max', 'rho1', 'tau_int']
    with pd.option_context('display.width', 200):
        print(res.sort_values(['system_type', 'chi'])[cols].to_string(index=False, float_format='%.3g'))
    print(f"\nCociente SEM(bloques)/SEM(ingenuo): mediana {res['cociente'].median():.2f}, "
          f"rango {res['cociente'].min():.2f}-{res['cociente'].max():.2f}")
    print(f"rho(1): mediana {res['rho1'].median():.3f}, rango {res['rho1'].min():.3f} a {res['rho1'].max():.3f}")
    print(f"Tabla: {csv}")
    if len(res) <= 20:
        pdf = graficar(res, curvas, os.path.join(args.out_dir, f"sem_bloques_{etiqueta}.pdf"))
        print(f"Figura: {pdf}")
    else:
        print("(mas de 20 simulaciones: no se grafica; filtrar con --sistema / --area_ratio / --chi)")


if __name__ == '__main__':
    main()
