#!/usr/bin/env python3
"""
Figuras del Supplemental Material, con los filtros exactos del paper
(R_a = 0.1; chi = 0.2 en las figuras a concentracion fija) y la notacion del
texto (L_arch, h_arch, ell_arch). No toca las figuras del paper principal.

Genera en 5_graficos_finales/:
  SM-matrix-linear.pdf         <s> (escala lineal) vs chi, cuatro formas
  SM-clogging-zoom.pdf         2x2: ln(<s>+1) vs chi en 0 <= chi <= 0.2, una
                               curva por R_a, con el valor del sistema puro
  SM-avalanche-boxplot.pdf     distribucion de s a chi = 0.2
  SM-tortuosity-vs-chi.pdf     tortuosidad media del arco vs chi (con chi = 0)
  SM-arch-geometry-chi0.2.pdf  2x2 a chi = 0.2: (a) tortuosidad, (b) ell_arch
                               vs L_arch, (c) L_arch/2h_arch, (d) h_arch

Las demas figuras del Supplemental salen de otros scripts (ver
supplementary.tex): correlacion_avalanchas.py, comparar_presion_por_chi.py,
graficos_finales.py --solo_composicion, temporal_chi.py.

Uso:
  python3 figuras_suplementario.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import graficos_finales as gf  # noqa: E402
from graficos_finales import CFG, AR_COL  # noqa: E402

SISTEMAS = ['Tri(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Disco(G)-Disco(C)']
R_A, R_A_TOL = 0.1, 0.02
CHI_FIJO, CHI_TOL = 0.2, 1e-3
LETRAS = ['(a)', '(b)', '(c)', '(d)']
# chi intermedios de otra tanda de simulaciones: la Fig. 7 del paper
# (graficos_finales.plot_sweeps_by_system sin zoom) los excluye, y aca tambien
CHIS_EXCLUIDOS = [0.07, 0.075, 0.085, 0.09, 0.095, 0.105, 0.11, 0.15]


def sin_chis_excluidos(df):
    chi = df['chi'].astype(float).to_numpy()
    return df[~np.isclose(chi[:, None], CHIS_EXCLUIDOS, atol=1e-6).any(axis=1)]


def guardar(fig, nombre):
    os.makedirs(CFG.out_dir, exist_ok=True)
    path = os.path.join(CFG.out_dir, nombre)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print(f"  -> {path}")


def letra(ax, txt):
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha='left', va='top', fontsize='large')


def etiqueta(sis):
    return CFG.english_labels.get(sis, sis)


def color(sis):
    return CFG.system_style.get(sis, ('black', 'o'))


def con_ra(df, r_a=R_A):
    return df[np.isclose(df[AR_COL].astype(float), r_a, atol=R_A_TOL)]


# ---------------------------------------------------------------------------
# Avalanchas (estadisticas_flujo_agregadas.csv)
# ---------------------------------------------------------------------------
def fig_matrix_linear(df):
    fig, ax = plt.subplots()
    for sis in SISTEMAS:
        d = con_ra(df[(df['sistema'] == sis) & (df['chi'] > 0)])
        if d.empty:
            continue
        dd = d['D_d_sq'].value_counts().idxmax()
        d = d[d['D_d_sq'] == dd].groupby('chi').agg(S=('S_mean', 'mean'), e=('S_sem', 'mean')).reset_index()
        s0, e0 = gf.get_pure_ref(df, sis, dd, 'S_mean')
        x = ([0.0] if s0 is not None else []) + d['chi'].tolist()
        y = ([s0] if s0 is not None else []) + d['S'].tolist()
        e = ([e0] if s0 is not None else []) + d['e'].tolist()
        c, m = color(sis)
        ax.errorbar(x, y, yerr=e, fmt=f'{m}-', color=c, label=etiqueta(sis), **gf.ERRORBAR_KW)
    ax.set_xlabel(r"Concentration $\chi$")
    ax.set_ylabel(r"Mean Avalanche Size $\langle s \rangle$")
    ax.legend(loc='upper left')
    guardar(fig, "SM-matrix-linear.pdf")


def fig_clogging_zoom(df, chi_max=0.2):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    colores, marcadores = gf.SEQ_STYLE(3)
    for ax, sis, let in zip(axes.flat, SISTEMAS, LETRAS):
        d = sin_chis_excluidos(df[(df['sistema'] == sis) & (df['chi'] > 0) & (df['chi'] <= chi_max + 1e-9)])
        if d.empty:
            ax.set_visible(False)
            continue
        dd = d['D_d_sq'].value_counts().idxmax()
        d = d[d['D_d_sq'] == dd]
        p0, e0 = gf.get_pure_ref(df, sis, dd, 'log_S')
        if p0 is not None:
            ax.axhline(p0, **gf.REF_LINE_KW, label=r'$\chi = 0$ (pure)')
        for r_a, c, m in zip((0.1, 0.2, 0.4), colores, marcadores):
            g = con_ra(d, r_a)
            if g.empty:
                continue
            g = g.groupby('chi').agg(y=('log_S', 'mean'), e=('log_S_err', 'mean')).reset_index()
            x = ([0.0] if p0 is not None else []) + g['chi'].tolist()
            y = ([p0] if p0 is not None else []) + g['y'].tolist()
            e = ([e0] if p0 is not None else []) + g['e'].tolist()
            ax.errorbar(x, y, yerr=e, fmt=f'{m}-', color=c, label=rf"$R_{{\rm a}} = {r_a:g}$", **gf.ERRORBAR_KW)
        ax.set_xlim(-0.005, chi_max + 0.005)
        ax.set_xlabel(r"Concentration $\chi$")
        ax.set_ylabel(r"$\ln(\langle s \rangle + 1)$")
        ax.legend(title=etiqueta(sis), loc='lower right', fontsize='small', title_fontsize='small')
        letra(ax, let)
    fig.tight_layout()
    guardar(fig, "SM-clogging-zoom.pdf")


def fig_boxplot(av_path, index_path):
    av = pd.read_csv(av_path)
    idx = pd.read_csv(index_path)
    meta = [c for c in ('sim_hash', 'system_type', 'chi', AR_COL) if c in idx.columns]
    av = av.drop(columns=[c for c in meta if c != 'sim_hash' and c in av.columns]).merge(
        idx[meta].drop_duplicates('sim_hash'), on='sim_hash', how='left')
    av['sistema'] = av['system_type'].map(gf.SYSTEM_NAME_MAP)
    av = con_ra(av[av['sistema'].isin(SISTEMAS) & np.isclose(av['chi'].astype(float), CHI_FIJO, atol=CHI_TOL)])
    # simulaciones duplicadas en el indice (misma secuencia): una sola vez
    firmas = av.groupby('sim_hash')['tamano_big'].apply(lambda x: hash(tuple(x)))
    av = av[~av['sim_hash'].isin(firmas[firmas.duplicated()].index)]
    presentes = [s for s in SISTEMAS if s in av['sistema'].unique()]
    if not presentes:
        print("  [!] sin avalanchas a chi = 0.2 y R_a = 0.1")
        return
    fig, ax = plt.subplots()
    sns.boxplot(data=av, x='sistema', y='tamano_big', hue='sistema', order=presentes,
                palette={s: color(s)[0] for s in presentes}, legend=False, ax=ax)
    ax.set_yscale('log')
    ax.set_xticks(range(len(presentes)))
    ax.set_xticklabels([etiqueta(s) for s in presentes], rotation=15)
    ax.set_xlabel(None)
    ax.set_ylabel(r"Avalanche Size $s$")
    guardar(fig, "SM-avalanche-boxplot.pdf")


# ---------------------------------------------------------------------------
# Arcos (resultados_completos_arcos.csv)
# ---------------------------------------------------------------------------
def fig_tortuosity_vs_chi(df_arcos):
    mezclas = con_ra(df_arcos[df_arcos['sistema'].isin(SISTEMAS) & (df_arcos['chi'] > 0)])
    d = pd.concat([mezclas, gf.puros_como_mezcla(df_arcos, mezclas, SISTEMAS)], ignore_index=True)
    d = d[d['chi'] <= 0.4 + 1e-9]
    fig, ax = plt.subplots()
    for sis in SISTEMAS:
        g = d[d['sistema'] == sis]
        if g.empty:
            continue
        s = g.groupby('chi')['Tortuosidad'].agg(['mean', 'std', 'count']).reset_index()
        s = s[s['count'] >= 20]
        c, m = color(sis)
        ax.errorbar(s['chi'], s['mean'], yerr=s['std'] / np.sqrt(s['count']), fmt=f'{m}-', color=c,
                    label=etiqueta(sis), **gf.ERRORBAR_KW)
    ax.set_xlim(-0.01, 0.42)
    ax.set_xlabel(r"Concentration $\chi$")
    ax.set_ylabel("Arch Tortuosity")
    ax.legend(loc='upper right')
    guardar(fig, "SM-tortuosity-vs-chi.pdf")


def fig_geometria_chi_fijo(df_arcos):
    d = con_ra(df_arcos[df_arcos['sistema'].isin(SISTEMAS)
                        & np.isclose(df_arcos['chi'].astype(float), CHI_FIJO, atol=CHI_TOL)])
    presentes = [s for s in SISTEMAS if s in d['sistema'].unique()]
    if not presentes:
        print("  [!] sin arcos a chi = 0.2 y R_a = 0.1")
        return
    paleta = {s: color(s)[0] for s in presentes}
    ticks = [etiqueta(s).replace('-', '-\n') for s in presentes]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    (a, b), (c_, e) = axes

    sns.violinplot(data=d, x='sistema', y='Tortuosidad', hue='sistema', order=presentes, palette=paleta,
                   inner='quartile', legend=False, ax=a, cut=0)
    a.set_ylabel("Arch Tortuosity")

    for s in presentes:
        g = d[d['sistema'] == s][['Span', 'Longitud']].dropna()
        b.plot(g['Span'], g['Longitud'], 'o', ms=1.5, alpha=0.15, color=paleta[s], mec='none')
        if len(g) > 2:
            p = np.polyfit(g['Span'], g['Longitud'], 1)
            xx = np.linspace(g['Span'].min(), g['Span'].max(), 50)
            b.plot(xx, np.polyval(p, xx), '-', lw=2.5, color=paleta[s], label=etiqueta(s))
    b.set_xlabel(r"Arch Span $L_{\rm arch}/d$")
    b.set_ylabel(r"Arch Length $\ell_{\rm arch}/d$")
    b.legend(loc='lower right', fontsize='small')

    sns.violinplot(data=d, x='sistema', y='Relacion_Aspecto', hue='sistema', order=presentes, palette=paleta,
                   inner='quartile', legend=False, ax=c_, cut=0)
    c_.set_ylabel(r"Arch Aspect Ratio $L_{\rm arch}/2h_{\rm arch}$")

    sns.boxplot(data=d, x='sistema', y='Altura', hue='sistema', order=presentes, palette=paleta,
                legend=False, ax=e)
    e.set_ylabel(r"Arch Height $h_{\rm arch}/d$")

    for ax in (a, c_, e):
        ax.set_xticks(range(len(presentes)))
        ax.set_xticklabels(ticks, fontsize='small')
        ax.set_xlabel(None)
    for ax, let in zip((a, b, c_, e), LETRAS):
        letra(ax, let)
    fig.tight_layout()
    guardar(fig, "SM-arch-geometry-chi0.2.pdf")


def main():
    gf.aplicar_estilo()
    print("[+] Figuras del Supplemental (R_a = 0.1)")
    if os.path.exists(CFG.data_path):
        df = gf.prepare_flow_stats(pd.read_csv(CFG.data_path))
        fig_matrix_linear(df)
        fig_clogging_zoom(df)
    else:
        print(f"  [!] no esta {CFG.data_path}")
    if os.path.exists(CFG.raw_flow_path) and os.path.exists(gf.rutas.MASTER_INDEX):
        fig_boxplot(CFG.raw_flow_path, gf.rutas.MASTER_INDEX)
    if os.path.exists(CFG.arcos_path):
        df_arcos = gf.prepare_arcos(pd.read_csv(CFG.arcos_path))
        fig_tortuosity_vs_chi(df_arcos)
        fig_geometria_chi_fijo(df_arcos)
    else:
        print(f"  [!] no esta {CFG.arcos_path}")


if __name__ == '__main__':
    main()
