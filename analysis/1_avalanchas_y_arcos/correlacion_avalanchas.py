#!/usr/bin/env python3
"""
Correlacion entre avalanchas sucesivas de una misma corrida: mapa de primer
retorno (s_{i+1} vs s_i) y funcion de autocorrelacion rho(k).

Si las avalanchas son independientes, el mapa de primer retorno es una nube
sin estructura (la distribucion de s_{i+1} no depende de s_i) y rho(k) cae
dentro de la banda +- 2/sqrt(n) para todo k >= 1.

Las avalanchas de cada simulacion se ordenan por t_inicio (o por el orden de
las filas si el csv no lo tiene) y las simulaciones duplicadas en el indice se
cuentan una sola vez, igual que en sem_bloques.py.

Salidas (en 5_graficos_finales/ por defecto):
  primer_retorno_<etiqueta>.pdf   un panel por chi: s_{i+1} vs s_i (log-log)
  autocorrelacion_<etiqueta>.pdf  rho(k) vs k, una curva por chi, banda 95 %
Por pantalla, una tabla con rho(1) (Pearson) y el coeficiente de Spearman de
(s_i, s_{i+1}), que no depende de la escala de s.

Uso (disco-disco R_a = 0.1, tres concentraciones):
  python3 correlacion_avalanchas.py --sistema Disco-Disco --area_ratio 0.1 --chi 0.01 0.1 0.2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import rutas  # noqa: E402
from estilo_figuras import aplicar_estilo, guardar_figura, REF_LINE_KW, SEQ_STYLE  # noqa: E402
from sem_bloques import autocorrelacion  # noqa: E402


def cargar(av_path, index_path, sistema, area_ratio, chi):
    """Avalanchas filtradas, ordenadas en el tiempo, sin simulaciones duplicadas.
    Devuelve {sim_hash: (fila de metadatos, array de s en orden)}."""
    av = pd.read_csv(av_path)
    if 't_inicio' not in av.columns:
        av['t_inicio'] = av.groupby('sim_hash').cumcount().astype(float)
    idx = pd.read_csv(index_path)
    meta = [c for c in ('sim_hash', 'system_type', 'chi', 'area_ratio_nominal') if c in idx.columns]
    av = av.drop(columns=[c for c in meta if c != 'sim_hash' and c in av.columns]).merge(
        idx[meta].drop_duplicates('sim_hash'), on='sim_hash', how='left')
    if sistema:
        av = av[av['system_type'].isin(sistema)]
    if area_ratio:
        av = av[np.isclose(av['area_ratio_nominal'].to_numpy(float)[:, None], area_ratio).any(axis=1)]
    if chi:
        av = av[np.isclose(av['chi'].to_numpy(float)[:, None], chi).any(axis=1)]
    sims, vistas = {}, set()
    for h, g in av.groupby('sim_hash'):
        x = g.sort_values('t_inicio', kind='stable')['tamano_big'].to_numpy(float)
        firma = hash(tuple(x))
        if firma in vistas or len(x) < 10:
            continue
        vistas.add(firma)
        sims[h] = (g.iloc[0], x)
    return sims


def spearman(a, b):
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


def figura_primer_retorno(sims, out):
    aplicar_estilo()
    orden = sorted(sims, key=lambda h: (sims[h][0]['system_type'], sims[h][0]['chi']))
    n = len(orden)
    ncol = min(3, n)
    nfil = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nfil, ncol, figsize=(3.6 * ncol, 3.5 * nfil), squeeze=False)
    for ax, h in zip(axes.flat, orden):
        meta, x = sims[h]
        a, b = x[:-1], x[1:]
        ok = (a > 0) & (b > 0)
        ax.loglog(a[ok], b[ok], 'o', ms=3, alpha=0.5, color='black', mec='none')
        lo, hi = min(a[ok].min(), b[ok].min()), max(a[ok].max(), b[ok].max())
        ax.plot([lo, hi], [lo, hi], **REF_LINE_KW)
        ax.set_title(rf"$\chi = {meta['chi']:g}$,  $r_S = {spearman(a, b):+.2f}$", fontsize='medium')
        ax.set_xlabel(r"$s_i$")
        ax.set_ylabel(r"$s_{i+1}$")
        ax.set_aspect('equal', adjustable='datalim')
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    fig.tight_layout()
    return guardar_figura(out, fig)


def figura_autocorrelacion(sims, out, k_max):
    aplicar_estilo()
    orden = sorted(sims, key=lambda h: sims[h][0]['chi'])
    colores, marcadores = SEQ_STYLE(len(orden))
    fig, ax = plt.subplots()
    n_min = min(len(sims[h][1]) for h in orden)
    banda = 2.0 / np.sqrt(n_min)
    ax.axhspan(-banda, banda, color='0.85', lw=0)
    ax.axhline(0.0, **REF_LINE_KW)
    for h, c, m in zip(orden, colores, marcadores):
        meta, x = sims[h]
        # rho(0) = 1 por definicion (cada avalancha consigo misma)
        rho = np.concatenate([[1.0], autocorrelacion(x, k_max)])
        ax.plot(np.arange(len(rho)), rho, f'{m}-', color=c, ms=4, label=rf"$\chi = {meta['chi']:g}$")
    ax.set_xlim(-0.5, k_max + 0.5)
    ax.set_xlabel(r"Lag $k$ (avalanches)")
    ax.set_ylabel(r"$\rho(k)$")
    ax.legend(loc='best', fontsize='small', ncol=2)
    return guardar_figura(out, fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--av', default=os.path.join(rutas.AVALANCHE_DIR, 'resultados_completos_esqueleto.csv'))
    ap.add_argument('--index', default=rutas.MASTER_INDEX)
    ap.add_argument('--sistema', nargs='+', default=None)
    ap.add_argument('--area_ratio', type=float, nargs='+', default=None)
    ap.add_argument('--chi', type=float, nargs='+', default=None,
                    help='concentraciones a graficar (para el mapa de primer retorno conviene elegir pocas)')
    ap.add_argument('--k_max', type=int, default=20, help='lag maximo de rho(k) (default 20)')
    ap.add_argument('--etiqueta', default=None)
    ap.add_argument('--out_dir', default=rutas.GRAFICOS_FINALES_DIR)
    args = ap.parse_args()

    sims = cargar(args.av, args.index, args.sistema, args.area_ratio, args.chi)
    if not sims:
        sys.exit("[!] No quedaron simulaciones despues de filtrar (--sistema / --area_ratio / --chi).")

    filas = []
    for h, (meta, x) in sorted(sims.items(), key=lambda t: (t[1][0]['system_type'], t[1][0]['chi'])):
        rho = autocorrelacion(x, 1)
        filas.append({'sistema': meta['system_type'], 'chi': f"{meta['chi']:g}",
                      'rho_1': f"{rho[0]:+.3f}", 'spearman': f"{spearman(x[:-1], x[1:]):+.3f}",
                      'banda_95': f"{2 / np.sqrt(len(x)):.3f}"})
    print(pd.DataFrame(filas).to_string(index=False))

    etiqueta = args.etiqueta or '_'.join(filter(None, [
        '-'.join(args.sistema) if args.sistema else 'todos',
        f"Ra{'-'.join(f'{r:g}' for r in args.area_ratio)}" if args.area_ratio else None]))
    os.makedirs(args.out_dir, exist_ok=True)
    p1 = figura_primer_retorno(sims, os.path.join(args.out_dir, f"primer_retorno_{etiqueta}.pdf"))
    p2 = figura_autocorrelacion(sims, os.path.join(args.out_dir, f"autocorrelacion_{etiqueta}.pdf"), args.k_max)
    print(f"\n{p1}\n{p2}")


if __name__ == '__main__':
    main()
