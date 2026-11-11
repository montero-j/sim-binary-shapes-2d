"""chi real de cada simulacion, calculado desde params.in.

El generador de simulaciones (run_scripts/run_local.sh, sweep_mixure.sh)
calcula el area de un poligono como n r^2 tan(pi/n), es decir, toma r como
APOTEMA. Box2D (simulation/src/globalSetup.cpp) pone los vertices a distancia
r del centro, es decir, r es el radio CIRCUNSCRIPTO, y el area real es
(n/2) r^2 sin(2 pi/n) = cos^2(pi/n) * (area prevista).
  cuadrado 0.50, pentagono 0.65, triangulo 0.25, disco 1 (sin cambio).

Este script NO modifica nada: lee params.in, calcula las areas reales y
compara con el chi y el R_a de la carpeta.

  chi_real = m_chico / (m_chico + m_grande),  m = N * area_real * densidad
  Ra_real  = area_real_chico / area_real_grande

Uso (desde analysis/):
  python3 1_avalanchas_y_arcos/chi_real.py --dir /mnt/c/Users/Juli/Downloads/simulation_data
  python3 1_avalanchas_y_arcos/chi_real.py --dir RUTA1 RUTA2 --out_csv ../processed_data/chi_real.csv
"""
import argparse
import math
import os
import re

import pandas as pd

NOMBRE_FORMA = {1: 'Disk', 3: 'Triangle', 4: 'Square', 5: 'Pentagon', 6: 'Hexagon'}


def area_real(r, lados):
    """Area tal como la construye Box2D (vertices sobre el circulo de radio r)."""
    if lados <= 1:
        return math.pi * r ** 2
    return 0.5 * lados * r ** 2 * math.sin(2 * math.pi / lados)


def area_prevista(r, lados):
    """Area que supuso el generador (r como apotema)."""
    if lados <= 1:
        return math.pi * r ** 2
    return lados * r ** 2 * math.tan(math.pi / lados)


def leer_granos(params_path):
    """Lista de dicts (N, r, lados, dens) leida de params.in, en el mismo
    orden que globalSetup.cpp: tipo 0, tipo 1, ..."""
    with open(params_path) as f:
        tokens = f.read().split()
    try:
        i = tokens.index('noTipoGranos:')
    except ValueError:
        return None
    n_tipos = int(tokens[i + 1])
    granos, j = [], i + 2
    for _ in range(n_tipos):
        N, r, lados, dens = tokens[j:j + 4]
        granos.append(dict(N=int(N), r=float(r), lados=int(lados), dens=float(dens)))
        j += 6  # N r nLados dens fric rest
    return granos


def numero_de_ruta(path, claves):
    """Busca en la ruta un numero tipo Chi_0p1 / _X0p1 / _R0p1 / Ratio_0.1."""
    for clave in claves:
        m = re.search(clave + r'([0-9]+(?:[p.][0-9]+)?)', path)
        if m:
            return float(m.group(1).replace('p', '.'))
    return float('nan')


def analizar(sim_dir):
    granos = leer_granos(os.path.join(sim_dir, 'params.in'))
    if not granos:
        return None
    fila = dict(path=sim_dir,
                chi_carpeta=numero_de_ruta(sim_dir, [r'Chi_', r'_X', r'chi_']),
                Ra_carpeta=numero_de_ruta(sim_dir, [r'_R', r'Ratio_', r'ratio_']))
    if len(granos) == 1:
        g = granos[0]
        fila.update(sistema=NOMBRE_FORMA.get(g['lados'], f"{g['lados']}-gon"),
                    N_grande=g['N'], N_chico=0, chi_real=0.0, Ra_real=float('nan'),
                    Ra_prevista=float('nan'))
        return fila
    for g in granos:
        g['A'] = area_real(g['r'], g['lados'])
        g['A_prev'] = area_prevista(g['r'], g['lados'])
    grande, chico = sorted(granos[:2], key=lambda g: g['A'], reverse=True)
    m_g = grande['N'] * grande['A'] * grande['dens']
    m_c = chico['N'] * chico['A'] * chico['dens']
    m_g_prev = grande['N'] * grande['A_prev'] * grande['dens']
    m_c_prev = chico['N'] * chico['A_prev'] * chico['dens']
    fila.update(
        sistema=f"{NOMBRE_FORMA.get(grande['lados'], grande['lados'])}-"
                f"{NOMBRE_FORMA.get(chico['lados'], chico['lados'])}",
        N_grande=grande['N'], N_chico=chico['N'],
        chi_prevista=m_c_prev / (m_c_prev + m_g_prev),   # lo que quiso el generador
        chi_real=m_c / (m_c + m_g),
        Ra_prevista=chico['A_prev'] / grande['A_prev'],
        Ra_real=chico['A'] / grande['A'],
    )
    return fila


def buscar_simulaciones(raiz):
    for dirpath, _, files in os.walk(raiz):
        if 'params.in' in files:
            yield dirpath


def main():
    ap = argparse.ArgumentParser(description="chi y R_a reales (areas como en Box2D) de cada simulacion.")
    ap.add_argument('--dir', nargs='+', required=True, help="carpeta(s) a recorrer (busca params.in)")
    ap.add_argument('--out_csv', default=None, help="guardar la tabla (default: chi_real.csv en la primera --dir)")
    args = ap.parse_args()

    filas = [f for raiz in args.dir for d in buscar_simulaciones(raiz) if (f := analizar(d))]
    if not filas:
        raise SystemExit(f"No encontre params.in en {args.dir}")
    df = pd.DataFrame(filas)
    df['ruta_rel'] = [os.path.relpath(p, args.dir[0]) for p in df['path']]
    cols = [c for c in ('ruta_rel', 'sistema', 'N_grande', 'N_chico', 'chi_carpeta', 'chi_prevista',
                        'chi_real', 'Ra_carpeta', 'Ra_prevista', 'Ra_real') if c in df.columns]
    df = df.sort_values(['sistema', 'chi_carpeta', 'ruta_rel'])
    with pd.option_context('display.width', 250, 'display.max_rows', None,
                           'display.max_colwidth', 60):
        print(df[cols].to_string(index=False, float_format='%.4f'))

    out = args.out_csv or os.path.join(args.dir[0], 'chi_real.csv')
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    df.drop(columns='ruta_rel').to_csv(out, index=False)
    print(f"\n{len(df)} simulaciones. Tabla: {out}")


if __name__ == '__main__':
    main()
