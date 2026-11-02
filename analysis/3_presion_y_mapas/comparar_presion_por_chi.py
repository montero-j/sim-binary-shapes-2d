#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grafica en funcion de chi, un panel por cada geometria de grano grande,
comparando la region arco contra la roi (grandes en linea solida, chicas
en punteada), a partir de resumen_presion_por_atasco.csv (el crudo por
atasco, no el ya agregado por configuracion). Dos magnitudes, en figuras
separadas:

  - Presion media reducida  <P>d/mg    -> <out>_nuevo_<Forma>.pdf
  - Esfuerzo de corte medio <tau>d/mg  -> <out>_nuevo_tau_<Forma>.pdf
    (tau = radio del circulo de Mohr del tensor de cada particula, la
    parte deviatorica; promedio pesado por area, igual que P)

Por que desde el crudo y no desde resumen_presion_por_configuracion.csv:
si hay mas de una configuracion compartiendo el mismo chi (por ejemplo,
distintos anchos de orificio W), agregar solo por forma+chi aca evita que
esos puntos separados queden conectados en zigzag en vez de promediados
en un solo punto por chi.

Ademas, filtra por razon de areas (chica/grande) ANTES de agregar: las
carpetas "Ratio_X_..." codifican esa razon en el nombre, y las carpetas
"Chi_X_W_Y" (sin prefijo "Ratio_") son, confirmado con
check_ratio_tamano.py, la misma razon fisica que "Ratio_0p1" (0.158/0.5
al cuadrado ~ 0.0999) -- asi que se funden automaticamente dentro de ese
grupo.

Barras de error: dispersion ENTRE ATASCOS (cada atasco es una muestra;
las particulas de un mismo atasco no son independientes). Para cada
atasco se calcula su media en la region (Sum/Area de ese atasco) y, para
cada forma+chi, la dispersion de esas medias pesada por el area del
atasco en la region -- el mismo peso que usa el punto central, que es
Sum total / Area total:
  --error sem (default): error estandar = std / sqrt(n_eff),
                          n_eff = (sum w)^2 / sum w^2 (atascos efectivos)
  --error std           : desviacion estandar pesada
  --error ninguno       : sin barras
Ademas se guarda <out>_datos.csv con media, std, sem, n y n_eff de cada
punto.

Eje x (--xlim): mismo rango de chi en todas las figuras, default 0 a 0.4;
si una forma no tiene datos hasta el final, ese tramo queda vacio.

Leyenda (--leyenda): afuera (default, a la derecha del grafico, no tapa
datos), arriba (encima del grafico, 2 columnas) o adentro (loc='best').

Si el csv viene de presion_tensor_arch.py con lotes reanudables (columna
'lote'), solo se usan las filas de lotes confirmados en lotes_ok.txt
(misma carpeta que el csv): asi un lote a medio escribir de una corrida
cortada no aparece duplicado.

Uso:
    python3 comparar_presion_por_chi.py --csv resumen_presion_por_atasco.csv \\
        --out comparacion_chi.pdf
    # -> comparacion_chi_<Forma>.pdf (presion) y comparacion_chi_tau_<Forma>.pdf (corte)

    # Solo una magnitud, otra razon de areas, otras formas u otro rango de chi:
    python3 comparar_presion_por_chi.py --csv resumen_presion_por_atasco.csv \\
        --magnitud tau --ratio 0.3 --formas Disco Cuadrado --chi_max 0.4 --out comparacion_chi.pdf
"""
import argparse
import os
import re
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

# Estilo comun (estilo_figuras.py, base: graficos_finales.py / 5b-matrix-normalized)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # estilo_figuras.py vive en analysis/
from estilo_figuras import (aplicar_estilo, estilo_forma, guardar_figura,
                            REGION_STYLE, ERRORBAR_KW)

REGIONES = ('arco', 'roi')

# Sufijo de las figuras con la presion reducida <P>d/mg (m = masa de un grano
# grande, d = diametro de un disco grande; en las simulaciones m = d = 1, asi
# que el valor numerico es P/g)
SUFIJO_NUEVO = '_nuevo'

# magnitud -> (prefijo de las columnas del csv, etiqueta del eje y, sufijo del archivo)
MAGNITUDES = {
    'P':   ('SumP_reducida',   r'Reduced Mean Pressure $\langle P\rangle d / mg$', ''),
    'tau': ('SumTau_reducida', r'Reduced Mean Shear Stress $\langle \tau\rangle d / mg$', '_tau'),
}

RATIO_REGEX = re.compile(r'Ratio_([0-9]+p[0-9]+)')
# Carpetas "Chi_X_W_Y" sin prefijo "Ratio_" -> misma razon fisica que
# Ratio_0p1 (confirmado a mano: 0.158/0.5 al cuadrado = 0.0999). Si en el
# futuro aparece otra familia sin prefijo con otra razon, este valor por
# defecto ya no aplicaria y habria que revisarlo.
RATIO_SIN_PREFIJO = 0.1


def extraer_ratio_area(sim_path):
    m = RATIO_REGEX.search(str(sim_path))
    if m:
        return float(m.group(1).replace('p', '.'))
    return RATIO_SIN_PREFIJO


def filtrar_lotes_confirmados(df, csv_path):
    """Si el csv tiene la columna 'lote', se queda solo con los lotes de
    lotes_ok.txt (los que terminaron de guardarse)."""
    if 'lote' not in df.columns:
        return df
    f_ok = os.path.join(os.path.dirname(os.path.abspath(csv_path)), 'lotes_ok.txt')
    if not os.path.exists(f_ok):
        print(f"AVISO: el csv tiene columna 'lote' pero no encontre {f_ok}; se usan todas las filas.")
        return df
    with open(f_ok) as f:
        ok = {l.strip() for l in f if l.strip()}
    antes = len(df)
    df = df[df['lote'].astype(str).isin(ok)]
    if len(df) < antes:
        print(f"Se descartan {antes - len(df)} filas de lotes sin confirmar (corrida cortada).")
    return df


def _estadistica_pesada(x, w):
    """Media, desviacion estandar, error estandar y n de las medias por
    atasco x, pesadas por w (area del atasco en la region). Los atascos
    sin area en la region (w = 0, ej. sin chicas en el arco) no cuentan."""
    ok = np.isfinite(x) & (w > 0)
    x, w = x[ok], w[ok]
    n = len(x)
    if n == 0:
        return np.nan, np.nan, np.nan, 0, 0.0
    sw = w.sum()
    media = (w * x).sum() / sw
    n_eff = sw ** 2 / (w ** 2).sum()
    if n < 2 or n_eff <= 1:
        return media, np.nan, np.nan, n, n_eff
    # varianza pesada insesgada (pesos de "confiabilidad")
    var = (w * (x - media) ** 2).sum() / (sw - (w ** 2).sum() / sw)
    std = np.sqrt(var)
    return media, std, std / np.sqrt(n_eff), n, n_eff


def agregar_por_forma_y_chi(df, magnitudes, regiones=REGIONES):
    """Un punto por (forma_grande, chi), pooleando todos los atascos de
    cualquier configuracion/ancho que comparta ese chi. El valor central es
    Sum total / Area total (media pesada por area); la dispersion es la de
    las medias por atasco, con el mismo peso."""
    cols_necesarias = []
    for region in regiones:
        for t in ('B', 'S'):
            cols_necesarias.append(f'Area_{region}_{t}')
            for mag in magnitudes:
                cols_necesarias.append(f'{MAGNITUDES[mag][0]}_{region}_{t}')
    faltantes = [c for c in cols_necesarias if c not in df.columns]
    if faltantes:
        pista = (" Las columnas SumTau_* las escribe la version de presion_tensor_arch.py "
                 "con esfuerzo de corte: hay que volver a correrlo."
                 if any(c.startswith('SumTau') for c in faltantes) else "")
        raise SystemExit(f"Al csv le faltan columnas {faltantes}.{pista}")

    filas = []
    for (forma, chi), g in df.groupby(['forma_grande', 'chi'], dropna=False):
        fila = {'forma_grande': forma, 'chi': chi, 'n_atascos': len(g)}
        for mag in magnitudes:
            pref = MAGNITUDES[mag][0]
            for region in regiones:
                for t in ('B', 'S'):
                    w = g[f'Area_{region}_{t}'].to_numpy(dtype=float)
                    suma = g[f'{pref}_{region}_{t}'].to_numpy(dtype=float)
                    with np.errstate(invalid='ignore', divide='ignore'):
                        x = np.where(w > 0, suma / w, np.nan)
                    media, std, sem, n, n_eff = _estadistica_pesada(x, w)
                    clave = f'{mag}_{region}_{t}'
                    fila[f'mean_{clave}'] = media
                    fila[f'std_{clave}'] = std
                    fila[f'sem_{clave}'] = sem
                    fila[f'n_{clave}'] = n
                    fila[f'neff_{clave}'] = n_eff
        filas.append(fila)
    return pd.DataFrame(filas)


def limites_y(df_agg, formas, mags, error, regiones=REGIONES, margen=0.05):
    """(ymin, ymax) que abarcan todos los puntos y sus barras de error de
    las formas y magnitudes dadas, con un margen. Nunca por debajo de 0
    (presion y corte medios son positivos; una barra de std puede
    extenderse por debajo de 0 pero el eje no necesita mostrarlo)."""
    lo, hi = [], []
    sub = df_agg[df_agg['forma_grande'].isin(formas)]
    for mag in mags:
        for region in regiones:
            for t in ('B', 'S'):
                y = sub[f'mean_{mag}_{region}_{t}'].to_numpy(dtype=float)
                e = (np.zeros_like(y) if error == 'ninguno'
                     else np.nan_to_num(sub[f'{error}_{mag}_{region}_{t}'].to_numpy(dtype=float)))
                ok = np.isfinite(y)
                if ok.any():
                    lo.append(np.min(y[ok] - e[ok]))
                    hi.append(np.max(y[ok] + e[ok]))
    if not lo:
        return None
    ymin, ymax = min(lo), max(hi)
    pad = margen * (ymax - ymin if ymax > ymin else abs(ymax) or 1.0)
    return max(0.0, ymin - pad), ymax + pad


def graficar_una_forma(df_agg, forma, mag, etiqueta_chi, error='sem', desplazar=0.0,
                       ylim=None, regiones=REGIONES, leyenda='afuera', xlim=(0.0, 0.4)):
    """Una sola figura, para una sola forma de grano grande y una magnitud:
    arco vs roi, grandes vs chicas, en funcion de chi, con barras de error.

    Estilo de 5b-matrix-normalized: tamaño de figura por defecto, marcador +
    linea, capsize=4, leyenda adentro, sin titulo (la forma va como titulo de
    la leyenda, igual que en temporal_chi.py).
    Color = region (arco vermillion, roi azul). Grandes: marcador de la forma
    del grano, linea solida. Chicas (siempre discos): circulo vacio, linea
    discontinua."""
    sub = df_agg[df_agg['forma_grande'] == forma].sort_values('chi')
    chi = sub['chi'].to_numpy(dtype=float)
    paso = desplazar * (np.nanmax(chi) - np.nanmin(chi) if len(chi) > 1 else 1.0)
    label_forma, marker_forma = estilo_forma(forma)

    fig, ax = plt.subplots()
    series = [(region, t) for region in regiones for t in ('B', 'S')]
    for k, (region, t) in enumerate(series):
        color, label_region = REGION_STYLE[region]
        clave = f'{mag}_{region}_{t}'
        y = sub[f'mean_{clave}'].to_numpy(dtype=float)
        yerr = None if error == 'ninguno' else sub[f'{error}_{clave}'].to_numpy(dtype=float)
        # desplazamiento horizontal opcional, para que las barras no se tapen
        dx = (k - (len(series) - 1) / 2) * paso
        if t == 'B':
            estilo = dict(fmt=f'{marker_forma}-', color=color)
            label = f'{label_region}, large'
        else:
            estilo = dict(fmt='o--', color=color, markerfacecolor='white')
            label = f'{label_region}, small'
        # Las chicas no existen en chi=0 (Area=0 -> NaN); matplotlib corta
        # la linea ahi solo, no hace falta filtrar a mano.
        ax.errorbar(chi + dx, y, yerr=yerr, label=label, **estilo, **ERRORBAR_KW)

    if ylim is not None:
        ax.set_ylim(*ylim)
    # Mismo rango de chi en todas las figuras (aunque una forma no tenga
    # datos hasta el final), con un margen chico para no cortar los marcadores.
    if xlim is not None:
        margen = 0.025 * (xlim[1] - xlim[0])
        ax.set_xlim(xlim[0] - margen, xlim[1] + margen)
    ax.set_xlabel(etiqueta_chi)
    ax.set_ylabel(MAGNITUDES[mag][1])
    # Leyenda: 'afuera' (a la derecha del eje, nunca tapa datos), 'arriba'
    # (encima del eje, en 2 columnas) o 'adentro' (loc='best', como antes).
    if leyenda == 'afuera':
        ax.legend(title=label_forma, loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    elif leyenda == 'arriba':
        ax.legend(title=label_forma, loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=2,
                  borderaxespad=0.0)
    else:
        ax.legend(title=label_forma, loc='best')
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', default=None, help='resumen_presion_por_atasco.csv')
    parser.add_argument('--datos', default=None,
                        help='<out>_datos.csv de una corrida anterior (valores ya agregados por '
                             'forma y chi): grafica desde ahi sin necesitar --csv. --ratio no se '
                             'aplica (el archivo ya tiene una sola razon de areas); --formas, '
                             '--chi_min y --chi_max si.')
    parser.add_argument('--out', default='comparacion_chi.pdf',
                        help='nombre base; la extension define el formato (default pdf, '
                             'como graficos_finales.py)')
    parser.add_argument('--magnitud', choices=['P', 'tau', 'ambas'], default='ambas',
                        help="P (presion), tau (esfuerzo de corte) o ambas (default), "
                             "cada una en su propia figura")
    parser.add_argument('--ratio', type=float, default=0.1,
                        help='Razon de areas (chica/grande) a incluir (default: 0.1)')
    parser.add_argument('--formas', nargs='+', default=None,
                        help='Lista de forma_grande a incluir (default: todas las presentes en el csv)')
    parser.add_argument('--chi_min', type=float, default=None)
    parser.add_argument('--chi_max', type=float, default=None)
    parser.add_argument('--error', choices=['std', 'sem', 'ninguno'], default='sem',
                        help="Barras de error entre atascos: sem = error estandar de la media (default), "
                             "std = desviacion estandar, ninguno = sin barras")
    parser.add_argument('--xlim', type=float, nargs=2, default=[0.0, 0.4], metavar=('CHIMIN', 'CHIMAX'),
                        help='rango de chi del eje x, igual en todas las figuras (default 0 0.4)')
    parser.add_argument('--leyenda', choices=['afuera', 'arriba', 'adentro'], default='afuera',
                        help="Posicion de la leyenda: afuera (default, a la derecha del grafico), "
                             "arriba (encima, en 2 columnas) o adentro (donde matplotlib encuentre lugar)")
    parser.add_argument('--desplazar', type=float, default=0.0,
                        help="Desplazamiento horizontal entre series, como fraccion del rango de "
                             "chi (ej. 0.004), para que las barras no se superpongan. Default 0.")
    parser.add_argument('--escala_y', choices=['por_magnitud', 'comun', 'libre'],
                        default='por_magnitud',
                        help="por_magnitud (default): todas las figuras de P con el mismo eje y, y "
                             "todas las de tau con otro; comun: P y tau con el mismo eje y; "
                             "libre: cada figura con su propio eje y (como antes)")
    parser.add_argument('--ylim_P', type=float, nargs=2, default=None, metavar=('YMIN', 'YMAX'),
                        help='fuerza el eje y de las figuras de presion')
    parser.add_argument('--ylim_tau', type=float, nargs=2, default=None, metavar=('YMIN', 'YMAX'),
                        help='fuerza el eje y de las figuras de corte')
    parser.add_argument('--etiqueta_chi', default=r'Concentration $\chi$',
                        help='Texto del eje x')
    args = parser.parse_args()
    aplicar_estilo(verbose=True)

    magnitudes = ['P', 'tau'] if args.magnitud == 'ambas' else [args.magnitud]

    if (args.csv is None) == (args.datos is None):
        parser.error("hay que pasar --csv (resumen por atasco) o --datos (valores ya agregados), uno solo.")

    if args.datos is not None:
        # Valores ya agregados (media, std, sem, n por forma, chi, region y especie)
        df = pd.read_csv(args.datos)
        faltan = [f'mean_{m}_{r}_{t}' for m in magnitudes for r in REGIONES for t in ('B', 'S')
                  if f'mean_{m}_{r}_{t}' not in df.columns]
        if faltan:
            raise SystemExit(f"A {args.datos} le faltan las columnas {faltan} — revisa --magnitud.")
        print(f"Graficando desde valores ya agregados: {args.datos} (--ratio no se aplica).")
    else:
        df = pd.read_csv(args.csv, dtype={'lote': str, 'jam_id': str})
        df = filtrar_lotes_confirmados(df, args.csv)
        df['ratio_area'] = df['sim_path'].apply(extraer_ratio_area)
        df = df[np.isclose(df['ratio_area'], args.ratio, atol=1e-3)]
        if df.empty:
            raise SystemExit(f"No quedo ninguna fila con ratio_area ~= {args.ratio} — revisa --ratio.")
        print(f"Razon de areas = {args.ratio}: {df['sim_path'].nunique()} carpetas distintas incluidas "
              f"({df['config_id'].nunique()} config_id distintos).")

    if args.formas is not None:
        df = df[df['forma_grande'].isin(args.formas)]
    if args.chi_min is not None:
        df = df[df['chi'] >= args.chi_min]
    if args.chi_max is not None:
        df = df[df['chi'] <= args.chi_max + 1e-9]
    if df.empty:
        raise SystemExit("No quedaron filas despues de filtrar — revisa --formas/--chi_min/--chi_max.")

    df_agg = df.reset_index(drop=True) if args.datos is not None else agregar_por_forma_y_chi(df, magnitudes)
    formas = args.formas if args.formas is not None else sorted(df_agg['forma_grande'].dropna().unique())
    print("Atascos por forma y chi:")
    print(df_agg[['forma_grande', 'chi', 'n_atascos']].to_string(index=False))

    base, ext = os.path.splitext(args.out)
    # Sufijo para no pisar las figuras con la etiqueta vieja (<P>/g)
    base = base + SUFIJO_NUEVO
    if os.path.dirname(args.out):
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
    datos_path = f"{base}_datos.csv"
    df_agg.to_csv(datos_path, index=False)
    print(f"Valores graficados (media, std, sem, n, n_eff): {datos_path}")
    print(f"Barras de error: {args.error}")

    # Eje y compartido (incluye las barras de error de todas las formas)
    formas_ok = [f for f in formas if f in df_agg['forma_grande'].values]
    ylims = {}
    if args.escala_y == 'comun':
        lim = limites_y(df_agg, formas_ok, magnitudes, args.error)
        ylims = {m: lim for m in magnitudes}
    elif args.escala_y == 'por_magnitud':
        ylims = {m: limites_y(df_agg, formas_ok, [m], args.error) for m in magnitudes}
    if args.ylim_P is not None:
        ylims['P'] = tuple(args.ylim_P)
    if args.ylim_tau is not None:
        ylims['tau'] = tuple(args.ylim_tau)
    for m, lim in ylims.items():
        if lim is not None and m in magnitudes:
            print(f"Eje y {m}: [{lim[0]:.4g}, {lim[1]:.4g}]  "
                  f"(para fijarlo: --ylim_{m} {lim[0]:.4g} {lim[1]:.4g})")
    for mag in magnitudes:
        sufijo = MAGNITUDES[mag][2]
        for forma in formas:
            if forma not in df_agg['forma_grande'].values:
                print(f"(sin datos para {forma!r}, se salta)")
                continue
            fig = graficar_una_forma(df_agg, forma, mag, args.etiqueta_chi,
                                     error=args.error, desplazar=args.desplazar,
                                     ylim=ylims.get(mag), leyenda=args.leyenda,
                                     xlim=tuple(args.xlim))
            out_path = f"{base}{sufijo}_{forma}{ext}"
            guardar_figura(out_path, fig)
            print(f"Guardado: {out_path}")


if __name__ == '__main__':
    main()