#!/usr/bin/env python3
"""
Sec. SIII.C del Supplemental: ¿son exponenciales las distribuciones de tamaño
de avalancha? Complementa SM-avalanche-boxplot.pdf (figuras_suplementario.py),
del que reutiliza la lectura de datos (resultados_completos_esqueleto.csv +
master_sim_index.csv, R_a = 0.1, simulaciones duplicadas una sola vez), los
colores/marcadores por forma y el estilo.

s = tamano_big: granos GRANDES descargados entre dos atascos.

Para cada mezcla (Tri-, Cuad-, Pent- y Disco-Disco) y chi = 0 (sistema puro,
'Tri', 'Cuad', ... en el indice, con el mismo orificio que las mezclas) y
chi = 0.2:
  - CCDF empirica P(S >= s), sin binning: para cada valor s distinto,
    (numero de avalanchas con tamaño >= s) / N. Se juntan las avalanchas de
    todas las corridas de esa (forma, chi), como en el boxplot.
  - Ajuste de la cola s > <s>/2 (--frac_cola) a P(S >= s) = A exp(-s/s*):
    cuadrados minimos de ln P vs s sobre los valores distintos de s con
    P >= --p_min_n / N (la punta, con 1-4 avalanchas, es puro ruido).
    Error de s* por bootstrap (remuestreo de avalanchas). Como control, el
    estimador de maxima verosimilitud de una exponencial (geometrica) por
    encima del corte: s*_MV = 1 / ln(1 + 1/<s - s_c>), que para s* >> 1 es
    ~ <s - s_c> + 1/2 (propiedad sin memoria).
  - Desviacion para avalanchas chicas:
      exceso = 1 - A exp(-s_min/s*): fraccion de avalanchas que NO explica la
               exponencial de la cola extrapolada hasta el s mas chico (s_min = 1).
               Exponencial pura -> 0; > 0 = sobran avalanchas chicas.
      f(s < <s>/2) medida, contra 1 - exp(-1/2) = 0.39 de una exponencial
               pura con la misma media.
    Un exceso de chicas baja <s> sin cambiar la cola, asi que tambien se ve
    como s*/<s> > 1.

Salidas (en --out_dir, por defecto processed_data/5_graficos_finales/; no se
sobrescribe nada que exista salvo con --sobrescribir):
  SM-avalanche-ccdf.pdf          (a) P(S>=s) vs s, chi = 0.2;  (b) vs s/<s>, chi = 0 y 0.2, con exp(-x)
  SM-avalanche-ccdf_tabla.csv/.md forma, chi, N, <s>, s*, s*/<s>, ...
  SM-avalanche-ccdf_parrafo.md   borrador del parrafo en ingles (revisar)

Uso:
  python3 ccdf_avalanchas.py
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
from figuras_suplementario import gf, CFG, AR_COL, SISTEMAS, con_ra, color, etiqueta  # noqa: E402

CHIS = (0.0, 0.2)
CHI_TOL = 1e-3


# ---------------------------------------------------------------------------
# Datos (misma lectura que fs.fig_boxplot, mas el sistema puro)
# ---------------------------------------------------------------------------
def leer_avalanchas(av_path, index_path):
    av = pd.read_csv(av_path)
    idx = pd.read_csv(index_path)
    meta = [c for c in ('sim_hash', 'system_type', 'chi', AR_COL, 'w_salida') if c in idx.columns]
    av = av.drop(columns=[c for c in meta if c != 'sim_hash' and c in av.columns]).merge(
        idx[meta].drop_duplicates('sim_hash'), on='sim_hash', how='left')
    av['sistema'] = av['system_type'].map(gf.SYSTEM_NAME_MAP)
    chi = av['chi'].astype(float)
    mez = con_ra(av[av['sistema'].isin(SISTEMAS) & (chi > 0)])
    partes = [mez[np.isclose(mez['chi'].astype(float), c, atol=CHI_TOL)] for c in CHIS if c > 0]
    for sis in SISTEMAS:
        base = sis.split('-')[0]                      # 'Tri(G)'
        p = av[(av['sistema'] == base) & np.isclose(chi, 0.0)]
        if p.empty:
            print(f"  [!] sin sistema puro {base} (chi = 0)")
            continue
        w = mez.loc[mez['sistema'] == sis, 'w_salida'].dropna() if 'w_salida' in av.columns else []
        if len(w):
            pw = p[np.isclose(p['w_salida'].astype(float), float(w.mode().iloc[0]))]
            p = pw if not pw.empty else p
        partes.append(p.assign(sistema=sis))
    av = pd.concat(partes, ignore_index=True)
    # simulaciones duplicadas en el indice (misma secuencia): una sola vez
    firmas = av.groupby('sim_hash')['tamano_big'].apply(lambda x: hash(tuple(x)))
    av = av[~av['sim_hash'].isin(firmas[firmas.duplicated()].index)]
    av['chi_nom'] = np.where(np.isclose(av['chi'].astype(float), 0.0), 0.0, av['chi'].astype(float).round(3))
    return av[av['tamano_big'] > 0]


# ---------------------------------------------------------------------------
# CCDF y ajustes
# ---------------------------------------------------------------------------
def ccdf(s):
    """Valores distintos de s y P(S >= s) empirica (sin binning)."""
    s = np.sort(np.asarray(s, float))
    u, primero = np.unique(s, return_index=True)
    return u, (len(s) - primero) / len(s)


def ajuste_cola(s, frac_cola=0.5, p_min_n=5):
    """(A, s*) de ln P(S>=s) = ln A - s/s* con s > frac_cola <s> y P >= p_min_n/N."""
    u, P = ccdf(s)
    sc = frac_cola * np.mean(s)
    ok = (u > sc) & (P >= p_min_n / len(s))
    if ok.sum() < 3:
        return np.nan, np.nan
    pend, ord0 = np.polyfit(u[ok], np.log(P[ok]), 1)
    return (np.exp(ord0), -1.0 / pend) if pend < 0 else (np.nan, np.nan)


def s_star_mv(s, frac_cola=0.5):
    s = np.asarray(s, float)
    sc = frac_cola * s.mean()
    exc = s[s > sc] - s[s > sc].min()
    m = exc.mean() if len(exc) else np.nan
    return 1.0 / np.log1p(1.0 / m) if m > 0 else np.nan


def analizar(s, frac_cola, p_min_n, n_boot, rng):
    s = np.asarray(s, float)
    media = s.mean()
    A, ss = ajuste_cola(s, frac_cola, p_min_n)
    boot = []
    for _ in range(n_boot):
        _, b = ajuste_cola(rng.choice(s, len(s)), frac_cola, p_min_n)
        boot.append(b)
    boot = np.array(boot)
    err = np.nanstd(boot, ddof=1) if np.isfinite(boot).sum() > 2 else np.nan
    sc = frac_cola * media
    f_emp = np.mean(s < sc)
    exceso = 1.0 - A * np.exp(-s.min() / ss) if np.isfinite(ss) else np.nan
    return dict(N=len(s), s_medio=media, sem=s.std(ddof=1) / np.sqrt(len(s)),
                s_star=ss, err_s_star=err, s_star_sobre_s=ss / media, err_cociente=err / media,
                s_star_MV=s_star_mv(s, frac_cola), A=A, s_corte=sc,
                f_chicas_emp=f_emp, f_chicas_exp=1 - np.exp(-frac_cola),
                exceso_chicas=exceso, cv=s.std(ddof=1) / media)


# ---------------------------------------------------------------------------
# Figura
# ---------------------------------------------------------------------------
def figura(av, out, presentes):
    gf.aplicar_estilo()
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12, 4.8))
    for sis in presentes:
        c, m = color(sis)
        for chi, ls, lleno in ((0.2, '-', True), (0.0, '--', False)):
            s = av.loc[(av['sistema'] == sis) & np.isclose(av['chi_nom'], chi), 'tamano_big'].to_numpy(float)
            if len(s) == 0:
                continue
            u, P = ccdf(s)
            kw = dict(color=c, marker=m, ms=5, markevery=max(1, len(u) // 12),
                      mfc=c if lleno else 'none', mec=c, lw=1.5, ls=ls)
            if chi == 0.2:
                axa.plot(u, P, label=etiqueta(sis), **kw)
            axb.plot(u / s.mean(), P, label=f"{etiqueta(sis)}, $\\chi = {chi:g}$", **kw)
    x_max = min(axb.get_xlim()[1], 8)
    xx = np.linspace(0, x_max, 200)
    axb.plot(xx, np.exp(-xx), color='black', ls=':', lw=2, label=r'$e^{-x}$', zorder=5)
    n_max = av.groupby(['sistema', 'chi_nom']).size().max()
    for ax in (axa, axb):
        ax.set_yscale('log')
        ax.set_ylabel(r"$P(S \geq s)$")
        ax.set_ylim(0.5 / n_max, 1.3)
        ax.set_xlim(left=0)
    axb.set_xlim(right=x_max)
    axa.set_xlabel(r"Avalanche Size $s$")
    axb.set_xlabel(r"$s/\langle s \rangle$")
    # las curvas bajan de arriba a la izquierda a abajo a la derecha: la
    # esquina de abajo a la izquierda y la de arriba a la derecha quedan libres
    axa.legend(title=r'$\chi = 0.2$', loc='lower left', fontsize='small', title_fontsize='small')
    axb.legend(loc='lower left', fontsize='x-small', ncol=2)
    for ax, let in ((axa, '(a)'), (axb, '(b)')):
        ax.text(0.97, 0.97, let, transform=ax.transAxes, ha='right', va='top', fontsize='large')
    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Tabla y parrafo
# ---------------------------------------------------------------------------
COLS_MD = [('forma', 'Shape', '{}'), ('chi', 'chi', '{:g}'), ('N', 'N', '{:d}'),
           ('s_medio', '<s>', '{:.1f}'), ('s_star', 's*', '{:.1f}'), ('err_s_star', 'err s*', '{:.1f}'),
           ('s_star_sobre_s', 's*/<s>', '{:.2f}'), ('s_star_MV', 's* (MV)', '{:.1f}'), ('A', 'A', '{:.2f}'),
           ('exceso_chicas', 'excess small', '{:+.2f}'), ('f_chicas_emp', 'f(s<<s>/2)', '{:.2f}'),
           ('cv', 'std/<s>', '{:.2f}')]


def tabla_md(df):
    def fmt(v, f):
        if pd.isna(v):
            return '-'
        return f.format(int(v)) if f == '{:d}' else f.format(v)
    L = ['| ' + ' | '.join(t for _, t, _ in COLS_MD) + ' |', '|' + '|'.join('---' for _ in COLS_MD) + '|']
    L += ['| ' + ' | '.join(fmt(f[c], fm) for c, _, fm in COLS_MD) + ' |' for _, f in df.iterrows()]
    return '\n'.join(L) + '\n'


def parrafo(df, frac_cola, tol=0.15, umbral_exceso=0.05):
    """Borrador en ingles armado con los numeros; hay que leerlo y ajustarlo.
    s*/<s> dentro de 1 +- tol se toma como exponencial; exceso de chicas
    > umbral_exceso como desviacion a señalar."""
    r = df['s_star_sobre_s']
    ok = df[(r - 1).abs() <= tol]
    malos = df[(r - 1).abs() > tol]
    txt = (f"Figure SM-avalanche-ccdf shows the empirical complementary cumulative distribution "
           f"$P(S \\geq s)$ of avalanche sizes for the four mixtures at $R_a = 0.1$, for the pure systems "
           f"($\\chi = 0$) and for $\\chi = 0.2$. Fitting the tail ($s > {frac_cola:g}\\langle s \\rangle$) "
           f"to $A\\,e^{{-s/s^*}}$ gives $s^*/\\langle s \\rangle$ between {r.min():.2f} and {r.max():.2f} "
           f"(Table SX).")
    if malos.empty:
        txt += (" The tails are straight lines in semilog scale and, once $s$ is rescaled by "
                "$\\langle s \\rangle$, all curves collapse onto $e^{-x}$: the distributions are exponential "
                "and fully characterized by $\\langle s \\rangle$.")
    else:
        txt += (f" The tails are exponential in all cases, and {len(ok)} of the {len(df)} distributions collapse "
                f"onto $e^{{-x}}$ when $s$ is rescaled by $\\langle s \\rangle$ "
                f"($|s^*/\\langle s \\rangle - 1| \\leq {tol:g}$).")
    exc = df[df['exceso_chicas'] > umbral_exceso].sort_values('exceso_chicas', ascending=False)
    if len(exc):
        quien = ', '.join(f"{f['forma']} at $\\chi = {f['chi']:g}$ ({100 * f['exceso_chicas']:.0f}\\%)"
                          for _, f in exc.iterrows())
        txt += (f" The deviations come from an excess of very small avalanches with respect to the "
                f"exponential tail extrapolated down to $s = 1$: {quien}. These extra small avalanches lower "
                f"$\\langle s \\rangle$ and raise $s^*/\\langle s \\rangle$ above one, without changing "
                f"the exponential form of the tail.")
    else:
        txt += " No significant excess of small avalanches is observed."
    return txt + "\n"


def no_pisar(path, sobrescribir):
    if os.path.exists(path) and not sobrescribir:
        sys.exit(f"[!] {path} ya existe: no se sobrescribe (usar otro --out_dir o --sobrescribir).")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--av', default=CFG.raw_flow_path)
    ap.add_argument('--index', default=gf.rutas.MASTER_INDEX)
    ap.add_argument('--frac_cola', type=float, default=0.5, help='la cola es s > frac_cola * <s> (default 0.5)')
    ap.add_argument('--p_min_n', type=float, default=5,
                    help='en el ajuste solo puntos con P >= p_min_n / N (default 5)')
    ap.add_argument('--n_boot', type=int, default=500)
    ap.add_argument('--out_dir', default=CFG.out_dir)
    ap.add_argument('--sobrescribir', action='store_true')
    args = ap.parse_args()

    for p in (args.av, args.index):
        if not os.path.exists(p):
            sys.exit(f"[!] No esta {p} (correr main_pipeline.py o pasar --av / --index).")
    pdf = no_pisar(os.path.join(args.out_dir, "SM-avalanche-ccdf.pdf"), args.sobrescribir)
    av = leer_avalanchas(args.av, args.index)
    presentes = [s for s in SISTEMAS if s in set(av['sistema'])]
    if not presentes:
        sys.exit("[!] sin avalanchas con R_a = 0.1 a chi = 0 / 0.2")

    rng = np.random.default_rng(12345)
    filas = []
    for sis in presentes:
        for chi in CHIS:
            g = av[(av['sistema'] == sis) & np.isclose(av['chi_nom'], chi)]
            if g.empty:
                print(f"  [!] {etiqueta(sis)}: sin avalanchas a chi = {chi:g}")
                continue
            f = analizar(g['tamano_big'].to_numpy(float), args.frac_cola, args.p_min_n, args.n_boot, rng)
            filas.append(dict(forma=etiqueta(sis), chi=chi, n_corridas=g['sim_hash'].nunique(), **f))
    df = pd.DataFrame(filas)

    print(f"[+] {figura(av, pdf, presentes)}")
    base = os.path.join(args.out_dir, "SM-avalanche-ccdf")
    csv = no_pisar(base + "_tabla.csv", args.sobrescribir)
    df.to_csv(csv, index=False)
    md = no_pisar(base + "_tabla.md", args.sobrescribir)
    with open(md, 'w') as fh:
        fh.write(f"# Avalanche-size distributions, R_a = 0.1\n\nTail fit P(S>=s) = A exp(-s/s*) for "
                 f"s > {args.frac_cola:g}<s> and P >= {args.p_min_n:g}/N; err s* from {args.n_boot} bootstrap "
                 f"resamples. s* (MV): maximum-likelihood exponential above the same cut. "
                 f"excess small = 1 - A exp(-s_min/s*): fraction of avalanches not accounted for by the "
                 f"tail exponential extrapolated down to s_min (0 for a pure exponential). "
                 f"f(s<<s>/2): measured fraction below the cut (pure exponential: "
                 f"{1 - np.exp(-args.frac_cola):.2f}). std/<s> = 1 for an exponential.\n\n")
        fh.write(tabla_md(df))
    par = no_pisar(base + "_parrafo.md", args.sobrescribir)
    with open(par, 'w') as fh:
        fh.write(parrafo(df, args.frac_cola))
    print()
    print(tabla_md(df))
    print(parrafo(df, args.frac_cola))
    print(f"Tablas: {csv}\n        {md}\nParrafo: {par}")


if __name__ == '__main__':
    main()
