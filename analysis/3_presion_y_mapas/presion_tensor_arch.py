#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import argparse
import json
import ast
import numpy as np
import pandas as pd
import networkx as nx
import gc
import time
import collections
import multiprocessing as mp
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # especies.py vive en analysis/
from especies import especie_por_tipo

# Altura extra de la ROI por encima del grano mas alto del arco (unidades de
# simulacion; con el grano grande de radio 0.5, en diametros de grano grande).
# Antes era 10.0; ahora 0 por defecto (--roi_extra), igual que force_arch_v4.py.
ROI_EXTRA = 0.0

CHI_REGEX = re.compile(r'(?:Chi_|_X)([0-9]+(?:p[0-9]+)?)')

def poly_area(coords):
    x = coords[:, 0]
    y = coords[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def nombre_forma(n_sides):
    if pd.isna(n_sides): return "Ninguna"
    n = int(n_sides)
    if n == 1: return "Disco"
    if n == 3: return "Triangulo"
    if n == 4: return "Cuadrado"
    if n == 5: return "Pentagono"
    if n == 6: return "Hexagono"
    if n == 8: return "Octogono"
    return f"Poli{n}"

def parse_xy(filename):
    particles = {}
    with open(filename, "r") as f:
        for line in f:
            p = line.split()
            if not p: continue
            pid = int(p[0])
            if pid < 0: continue
            n_sides = int(p[1])
            shape_name = nombre_forma(n_sides)
            if n_sides == 1:
                x, y, r, ptype = np.float32(p[2]), np.float32(p[3]), np.float32(p[4]), int(p[5])
                area = np.pi * (r**2)
                particles[pid] = {"pos_x": x, "pos_y": y, "r_eff": r, "type": ptype, "shape": shape_name, "area": area, "coords": None}
            else:
                coords = np.array(p[2:-1], dtype=np.float32).reshape(-1, 2)
                center = coords.mean(axis=0)
                ptype = int(p[-1])
                r_eff = np.float32(np.max(np.linalg.norm(coords - center, axis=1)))
                area = poly_area(coords)
                particles[pid] = {"pos_x": center[0], "pos_y": center[1], "r_eff": r_eff, "type": ptype, "shape": shape_name, "area": area, "coords": coords}
    return particles

def parse_paredes(filename):
    """Segmentos de las paredes/piso del silo (lineas con id < 0 del .xy):
    'id n x1 y1 ... xn yn tipo' es una cadena abierta de n vertices (o un
    borde si n = 2). Devuelve un array (k, 2, 2) con los k segmentos
    consecutivos (no se cierra la cadena)."""
    segs = []
    with open(filename, "r") as f:
        for line in f:
            p = line.split()
            if len(p) < 2 or int(p[0]) >= 0:
                continue
            n = int(p[1])
            if n < 2 or len(p) < 2 + 2 * n:
                continue
            v = np.array(p[2:2 + 2 * n], dtype=np.float64).reshape(-1, 2)
            segs.extend(np.stack([v[:-1], v[1:]], axis=1))
    return np.array(segs, dtype=np.float64).reshape(-1, 2, 2)


def parse_radio_salida(sim_path, default=2.5):
    """Lee Radio_salida desde params.in en la carpeta de la simulacion
    (mismo archivo/campo que usa graficar_arco_y_fuerzas.py). Si no se
    encuentra el archivo o el campo, devuelve el valor por defecto."""
    params_file = os.path.join(sim_path, "params.in")
    r_salida = default
    try:
        with open(params_file, 'r') as f:
            for line in f:
                if 'Radio_salida:' in line:
                    r_salida = float(line.split(':')[1].split('#')[0].strip())
                    break
    except Exception:
        pass
    return r_salida


def parse_g(sim_path, default=9.8):
    """Lee 'g' desde params.in en la carpeta de la simulacion, para poder
    expresar la presion en unidades reducidas (P/g). params.in trae g con
    signo negativo (ej. 'g: -9.8', convencion de gravedad hacia abajo) ->
    devolvemos el valor absoluto."""
    params_file = os.path.join(sim_path, "params.in")
    g = default
    try:
        with open(params_file, 'r') as f:
            for line in f:
                if line.strip().startswith('g:'):
                    g = float(line.split(':')[1].split('#')[0].strip())
                    break
    except Exception:
        pass
    return abs(g) if g != 0 else default


def detectar_arco_primario(df_particles, df_contacts, r_salida):
    """Detecta el arco de atasco como el camino de costo minimo entre un
    grano apoyado en el borde izquierdo del orificio y uno apoyado en el
    borde derecho, con el costo de cada arista = 10**y_max (favorece
    fuertisimo los caminos bajos), y despues recorta de los dos extremos
    los tramos que son solo cimientos estables apoyados en el piso (no
    parte real del puente que sostiene la carga).

    """
    if df_particles.empty or df_contacts.empty:
        return []

    # Agrupamos por par de particulas,
    # sumando f_norm de los (hasta 2) puntos de contacto de un mismo par,
    # para tener UNA arista por par tocandose (no una por punto de contacto).
    valid_c = df_contacts[(df_contacts['id_a'] >= 0) & (df_contacts['id_b'] >= 0)]
    if valid_c.empty:
        return []
    # Un par por particulas en contacto (no uno por punto de contacto),
    # ordenado por (pair_a, pair_b) como el groupby original, asi el orden
    # en que se agregan las aristas (y el desempate de Dijkstra) no cambia.
    pares = pd.DataFrame({
        'pair_a': np.minimum(valid_c['id_a'].to_numpy(), valid_c['id_b'].to_numpy()),
        'pair_b': np.maximum(valid_c['id_a'].to_numpy(), valid_c['id_b'].to_numpy()),
    }).drop_duplicates().sort_values(['pair_a', 'pair_b'])

    idx = df_particles.index
    pos_y = df_particles['pos_y']
    ya = pos_y.reindex(pares['pair_a'].to_numpy()).to_numpy(dtype=np.float64)
    yb = pos_y.reindex(pares['pair_b'].to_numpy()).to_numpy(dtype=np.float64)
    ok = ~np.isnan(ya) & ~np.isnan(yb)
    # 10**y en float de python (float64): con y ~ 100 en float32 desborda.
    costos = [10.0 ** float(y) for y in np.maximum(ya[ok], yb[ok])]

    G = nx.Graph()
    G.add_nodes_from(idx)
    G.add_weighted_edges_from(zip(pares['pair_a'].to_numpy()[ok].tolist(),
                                  pares['pair_b'].to_numpy()[ok].tolist(), costos),
                              weight='costo_arco')

    r_col = 'r_eff'  # radio geometrico real, no el nominal ajustado por forma

    def es_cimiento_estable(pid, path, i):
        x = df_particles.at[pid, 'pos_x']
        y = df_particles.at[pid, 'pos_y']
        r = df_particles.at[pid, r_col]
        if abs(x) < (r_salida - 0.2):
            return False
        if (y - r) < 0.2:
            return True
        vecinos_arco = []
        if i > 0: vecinos_arco.append(path[i - 1])
        if i < len(path) - 1: vecinos_arco.append(path[i + 1])
        apoyos_inferiores = 0
        for vecino in G.neighbors(pid):
            if vecino not in vecinos_arco and df_particles.at[vecino, 'pos_y'] < (y - 0.05):
                apoyos_inferiores += 1
        return apoyos_inferiores >= 2

    def purgar_cimientos(path):
        if not path:
            return []
        iL = 0
        while iL < len(path) and es_cimiento_estable(path[iL], path, iL):
            iL += 1
        iR = len(path) - 1
        while iR >= iL and es_cimiento_estable(path[iR], path, iR):
            iR -= 1
        if iL > iR:
            return []
        return path[iL:iR + 1]

    apoyado = ((df_particles['pos_y'] - df_particles[r_col]) < 0.2).to_numpy()
    xs = df_particles['pos_x'].to_numpy()
    nleft = idx[apoyado & (xs < 0.0)].tolist()
    nright = idx[apoyado & (xs > 0.0)].tolist()


    mejor_path = None
    min_costo_total = float('inf')
    for l in nleft:
        dist, paths = nx.single_source_dijkstra(G, l, weight='costo_arco')
        for rr in nright:
            costo_path = dist.get(rr)
            if costo_path is not None and costo_path < min_costo_total:
                min_costo_total = costo_path
                mejor_path = paths[rr]

    return purgar_cimientos(mejor_path) if mejor_path is not None else []

def obtener_vector_normal_poly(row):
    px, py = row['px'], row['py']
    cx, cy = row['pos_x_A'], row['pos_y_A']
    coords = row['coords']
    n_vertices = len(coords)
    min_dist = float('inf')
    mejor_normal = (0.0, 0.0)
    for i in range(n_vertices):
        x1, y1 = coords[i]
        x2, y2 = coords[(i+1) % n_vertices]
        dx, dy = x2 - x1, y2 - y1
        l2 = dx*dx + dy*dy
        if l2 == 0: continue
        t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / l2))
        proy_x, proy_y = x1 + t * dx, y1 + t * dy
        dist2 = (px - proy_x)**2 + (py - proy_y)**2
        if dist2 < min_dist:
            min_dist = dist2
            nx_v, ny_v = -dy, dx
            norma = np.hypot(nx_v, ny_v)
            if norma > 0:
                nx_v, ny_v = nx_v/norma, ny_v/norma
            if (nx_v * (proy_x - cx) + ny_v * (proy_y - cy)) < 0:
                nx_v, ny_v = -nx_v, -ny_v
            mejor_normal = (nx_v, ny_v)
    return mejor_normal

def obtener_vector_normal_poly_lado(row, lado):
    """Igual que obtener_vector_normal_poly, pero puede leer el lado A o B
    de una fila ya mergeada con sufijos (_A/_B), para poder usar el lado
    que corresponda cuando el otro lado es la pared (sin datos de geometria)."""
    px, py = row['px'], row['py']
    cx, cy = row[f'pos_x_{lado}'], row[f'pos_y_{lado}']
    coords = row[f'coords_{lado}']
    n_vertices = len(coords)
    min_dist = float('inf')
    mejor_normal = (0.0, 0.0)
    for i in range(n_vertices):
        x1, y1 = coords[i]
        x2, y2 = coords[(i+1) % n_vertices]
        dx, dy = x2 - x1, y2 - y1
        l2 = dx*dx + dy*dy
        if l2 == 0: continue
        t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / l2))
        proy_x, proy_y = x1 + t * dx, y1 + t * dy
        dist2 = (px - proy_x)**2 + (py - proy_y)**2
        if dist2 < min_dist:
            min_dist = dist2
            nx_v, ny_v = -dy, dx
            norma = np.hypot(nx_v, ny_v)
            if norma > 0:
                nx_v, ny_v = nx_v/norma, ny_v/norma
            if (nx_v * (proy_x - cx) + ny_v * (proy_y - cy)) < 0:
                nx_v, ny_v = -nx_v, -ny_v
            mejor_normal = (nx_v, ny_v)
    return mejor_normal


def normales_poligonos(px, py, cx, cy, coords, con_interior=False):
    """Version vectorizada de obtener_vector_normal_poly_lado para muchos
    contactos a la vez: para cada contacto busca el lado del poligono mas
    cercano al punto de contacto y devuelve su normal unitaria hacia afuera.
    Misma logica que la version fila por fila (proyeccion acotada al
    segmento, primer minimo, orientacion segun el centro), pero con numpy
    sobre todos los contactos de poligonos con igual cantidad de vertices."""
    m = len(px)
    out = np.zeros((m, 2))
    interior = np.zeros(m)
    if m == 0:
        return (out, interior) if con_interior else out
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    coords = list(coords)
    nverts = np.array([len(c) for c in coords])
    for n in np.unique(nverts):
        sel = np.flatnonzero(nverts == n)
        V = np.stack([np.asarray(coords[i], dtype=np.float64) for i in sel])  # (k, n, 2)
        d = np.roll(V, -1, axis=1) - V                                        # lados
        l2 = (d ** 2).sum(-1)
        P = np.stack([px[sel], py[sel]], -1)[:, None, :]
        with np.errstate(invalid='ignore', divide='ignore'):
            t = np.clip(((P - V) * d).sum(-1) / l2, 0.0, 1.0)
        proy = V + t[..., None] * d
        dist2 = ((P - proy) ** 2).sum(-1)
        dist2 = np.where(l2 == 0, np.inf, dist2)
        i = np.argmin(dist2, axis=1)          # primer minimo, igual que el '<' original
        filas = np.arange(len(sel))
        ds = d[filas, i]
        ps = proy[filas, i]
        nv = np.stack([-ds[:, 1], ds[:, 0]], -1)
        norma = np.hypot(nv[:, 0], nv[:, 1])
        nv = np.where(norma[:, None] > 0, nv / np.where(norma > 0, norma, 1.0)[:, None], nv)
        C = np.stack([cx[sel], cy[sel]], -1)
        flip = ((nv * (ps - C)).sum(-1)) < 0
        nv[flip] *= -1
        nv[np.isinf(dist2).all(axis=1)] = 0.0   # poligono degenerado -> (0, 0) como antes
        out[sel] = nv
        # distancia de la proyeccion al vertice mas cercano del lado elegido:
        # ~0 si el contacto esta sobre un vertice (normal ambigua), grande si
        # esta en el medio de una cara (cara de referencia de Box2D)
        ts = t[filas, i]
        interior[sel] = np.minimum(ts, 1.0 - ts) * np.sqrt(l2[filas, i])
    return (out, interior) if con_interior else out


def normales_paredes(px, py, cx, cy, segs):
    """Para cada contacto con una pared: segmento de pared mas cercano al
    punto de contacto, su normal unitaria orientada hacia el cuerpo de centro
    (cx, cy), y la 'interioridad' del punto en ese segmento (igual que en
    normales_poligonos)."""
    m = len(px)
    out = np.zeros((m, 2))
    interior = np.zeros(m)
    if m == 0 or segs is None or len(segs) == 0:
        return out, np.full(m, -np.inf)
    P = np.stack([np.asarray(px, float), np.asarray(py, float)], -1)[:, None, :]   # (m, 1, 2)
    V = segs[None, :, 0, :]                                                         # (1, k, 2)
    d = (segs[:, 1, :] - segs[:, 0, :])[None]                                       # (1, k, 2)
    l2 = (d ** 2).sum(-1)
    with np.errstate(invalid='ignore', divide='ignore'):
        t = np.clip(((P - V) * d).sum(-1) / l2, 0.0, 1.0)
    proy = V + t[..., None] * d
    dist2 = np.where(l2 == 0, np.inf, ((P - proy) ** 2).sum(-1))
    i = np.argmin(dist2, axis=1)
    filas = np.arange(m)
    ds = d[0, i]
    nv = np.stack([-ds[:, 1], ds[:, 0]], -1)
    nv /= np.where(np.hypot(nv[:, 0], nv[:, 1]) > 0, np.hypot(nv[:, 0], nv[:, 1]), 1.0)[:, None]
    C = np.stack([np.asarray(cx, float), np.asarray(cy, float)], -1)
    flip = ((nv * (C - proy[filas, i])).sum(-1)) < 0
    nv[flip] *= -1
    out[:] = nv
    ts = t[filas, i]
    interior[:] = np.minimum(ts, 1.0 - ts) * np.sqrt(l2[0, i])
    return out, interior


def normales_contacto(df, paredes=None, metodo='nuevo'):
    """Normal unitaria de cada contacto, orientada de A hacia B (convencion
    de Box2D). df es el merge de contactos con las particulas (columnas
    _A/_B); una pared (id < 0) tiene shape NaN.

    metodo='nuevo': si alguno de los dos es disco, la direccion radial del
      disco (exacta). Si no (poligono-poligono o poligono-pared), se calcula
      la normal de la superficie de cada lado y se usa la del lado donde el
      punto de contacto cae en el interior de una cara (cara de referencia),
      no sobre un vertice.
    metodo='viejo': la de la version anterior (lado A si A es grano, si no
      lado B; en poligonos, siempre el lado mas cercano de ese cuerpo)."""
    m = len(df)
    px = df['px'].to_numpy(np.float64)
    py = df['py'].to_numpy(np.float64)
    ax_, ay_ = df['pos_x_A'].to_numpy(np.float64), df['pos_y_A'].to_numpy(np.float64)
    bx_, by_ = df['pos_x_B'].to_numpy(np.float64), df['pos_y_B'].to_numpy(np.float64)
    shA, shB = df['shape_A'], df['shape_B']
    a_ok, b_ok = shA.notna().to_numpy(), shB.notna().to_numpy()
    discoA = (shA == 'Disco').to_numpy()
    discoB = (shB == 'Disco').to_numpy()
    polyA = a_ok & ~discoA
    polyB = b_ok & ~discoB

    def radial(cx, cy, hacia_afuera):
        # hacia_afuera=True: del centro al contacto; False: del contacto al centro
        vx, vy = (px - cx, py - cy) if hacia_afuera else (cx - px, cy - py)
        nr = np.hypot(vx, vy)
        nr = np.where(nr > 0, nr, 1.0)
        return np.stack([vx / nr, vy / nr], -1)

    # Candidatos por lado, todos ya orientados de A hacia B, con su interioridad
    nA = np.zeros((m, 2)); iA = np.full(m, -np.inf)
    nB = np.zeros((m, 2)); iB = np.full(m, -np.inf)
    if discoA.any():
        nA[discoA] = radial(ax_, ay_, True)[discoA]; iA[discoA] = np.inf
    if discoB.any():
        nB[discoB] = radial(bx_, by_, False)[discoB]; iB[discoB] = np.inf
    for lado, mask, cx, cy, signo, nn, ii in (('A', polyA, ax_, ay_, +1.0, nA, iA),
                                              ('B', polyB, bx_, by_, -1.0, nB, iB)):
        if mask.any():
            nv, it = normales_poligonos(px[mask], py[mask], cx[mask], cy[mask],
                                        df[f'coords_{lado}'].to_numpy()[mask], con_interior=True)
            nn[mask] = signo * nv          # normal hacia afuera de B -> se invierte
            ii[mask] = it
    # Paredes: normal del segmento orientada hacia el otro cuerpo (= de A hacia B si la pared es A)
    pared_A = ~a_ok & b_ok
    pared_B = a_ok & ~b_ok
    if paredes is not None and len(paredes) > 0:
        if pared_A.any():
            nv, it = normales_paredes(px[pared_A], py[pared_A], bx_[pared_A], by_[pared_A], paredes)
            nA[pared_A] = nv; iA[pared_A] = it
        if pared_B.any():
            nv, it = normales_paredes(px[pared_B], py[pared_B], ax_[pared_B], ay_[pared_B], paredes)
            nB[pared_B] = -nv; iB[pared_B] = it

    if metodo == 'viejo':
        usar_A = a_ok
        n = np.where(usar_A[:, None], nA, nB)
        # en la version vieja, si A era pared se usaba solo el lado B
        return n[:, 0], n[:, 1]
    usar_A = iA >= iB
    n = np.where(usar_A[:, None], nA, nB)
    return n[:, 0], n[:, 1]


def leer_normales_box2d(contact_file, df_contacts):
    """Normales de Box2D (de A hacia B) de atasco_contacts_full_<jam>.dat, el
    archivo que escribe el simulador desde la rama tensor-en-simulacion
    (simulation/README_tensor_atasco.md). Tiene las mismas filas, en el mismo
    orden, que atasco_contacts_<jam>.dat.

    Devuelve (normales, motivo): normales es un array (m, 2), o None si el
    archivo no existe o no coincide fila por fila con df_contacts (en ese
    caso se usa la normal inferida de la geometria, como con los datos
    viejos); motivo explica cual de los dos casos se dio."""
    d, nombre = os.path.split(contact_file)
    full = os.path.join(d, nombre.replace('atasco_contacts_', 'atasco_contacts_full_', 1))
    if not os.path.exists(full):
        return None, 'inferidas (sin atasco_contacts_full)'
    cols = None
    with open(full) as f:
        for linea in f:
            if not linea.startswith('#'):
                break
            cols = linea[1:].split()
    df_full = pd.read_csv(full, sep=r'\s+', comment='#', names=cols,
                          usecols=['ID_A', 'ID_B', 'nx', 'ny'], engine='c')
    if (len(df_full) != len(df_contacts)
            or not np.array_equal(df_full['ID_A'].to_numpy(), df_contacts['id_a'].to_numpy())
            or not np.array_equal(df_full['ID_B'].to_numpy(), df_contacts['id_b'].to_numpy())):
        return None, 'inferidas (atasco_contacts_full no coincide fila por fila)'
    return df_full[['nx', 'ny']].to_numpy(np.float64), 'box2d'


def fuerzas_contacto(df_particles, df_contacts, dt, paredes=None, metodo='nuevo',
                     signo_tan=1.0, dividir_dt=True, normales=None):
    """Fuerza de cada contacto sobre el cuerpo A (Fx_A, Fy_A), con la normal
    de A hacia B. La fuerza sobre B es -F_A (3ra ley). Box2D aplica sobre B
    el impulso fn*n + ft*t_b2d, t_b2d = (n_y, -n_x) = -t con t = (-n_y, n_x);
    sobre A: -fn*n - ft*t_b2d = -fn*n + ft*t. Verificado con el equilibrio
    estatico de cada grano (analysis/utils/verificar_contactos.py).
    normales: array (m, 2) con la normal de Box2D de cada fila de
    df_contacts (ver leer_normales_box2d); si es None, se infiere de la
    geometria con normales_contacto. La inferencia falla en los contactos
    esquina contra esquina entre poligonos (simulation/README_tensor_atasco.md).
    signo_tan y dividir_dt solo existen para el diagnostico."""
    part_cols = df_particles[['pos_x', 'pos_y', 'shape', 'coords']]
    df = df_contacts.merge(part_cols, left_on='id_a', right_index=True, how='left')
    df = df.merge(part_cols, left_on='id_b', right_index=True, how='left', suffixes=('_A', '_B'))
    if normales is not None:
        nx, ny = normales[:, 0], normales[:, 1]
    else:
        nx, ny = normales_contacto(df, paredes, metodo)
    tx, ty = -ny, nx
    fn = df['f_norm'].to_numpy(np.float64)
    ft = df['f_tan'].to_numpy(np.float64)
    esc = 1.0 / dt if dividir_dt else 1.0
    df['Fx_A'] = (-fn * nx + signo_tan * ft * tx) * esc
    df['Fy_A'] = (-fn * ny + signo_tan * ft * ty) * esc
    df['a_valida'] = df['shape_A'].notna()
    df['b_valida'] = df['shape_B'].notna()
    return df


def calcular_tensor_y_presion_vectorizado(df_particles, df_contacts, arco_path, dt, umbral_r, g=9.8,
                                          paredes=None, normales=None):
    # Fuerza de cada contacto sobre A, con la normal de A hacia B (ver
    # normales_contacto y fuerzas_contacto); la fuerza sobre B es -F_A.
    # normales: las de Box2D si estan disponibles (leer_normales_box2d).
    # Los contactos contra pared (id < 0, shape NaN) se mantienen: aportan
    # solo al tensor del grano, no se le "inventa" area a la pared.
    df = fuerzas_contacto(df_particles, df_contacts, dt, paredes=paredes, normales=normales)
    a_valida = df['a_valida'].to_numpy()
    b_valida = df['b_valida'].to_numpy()
    Fx, Fy = df['Fx_A'].to_numpy(), df['Fy_A'].to_numpy()
    rax, ray = (df['px'] - df['pos_x_A']).to_numpy(), (df['py'] - df['pos_y_A']).to_numpy()
    rbx, rby = (df['px'] - df['pos_x_B']).to_numpy(), (df['py'] - df['pos_y_B']).to_numpy()
    # S_ij = sum_c r_i F_j  (r = punto de contacto - centro de la particula):
    #   S00 = r_x F_x, S11 = r_y F_y, S01 = r_x F_y, S10 = r_y F_x
    # Sobre A: F_A ; sobre B: -F_A.
    for nombre, (ra, rb, F) in {'S00': (rax, rbx, Fx), 'S11': (ray, rby, Fy),
                                'S01': (rax, rbx, Fy), 'S10': (ray, rby, Fx)}.items():
        df[f'{nombre}_A'] = np.where(a_valida, ra * F, 0.0).astype(np.float32)
        df[f'{nombre}_B'] = np.where(b_valida, -rb * F, 0.0).astype(np.float32)

    comps_A = {'S00_A': 'S00', 'S11_A': 'S11', 'S01_A': 'S01', 'S10_A': 'S10'}
    comps_B = {'S00_B': 'S00', 'S11_B': 'S11', 'S01_B': 'S01', 'S10_B': 'S10'}
    tensor_A = df.groupby('id_a')[list(comps_A)].sum().rename(columns=comps_A)
    tensor_B = df.groupby('id_b')[list(comps_B)].sum().rename(columns=comps_B)
    tensor_total = tensor_A.add(tensor_B, fill_value=0)
    
    del df, tensor_A, tensor_B
    
    p_df = df_particles.join(tensor_total).fillna({'S00': 0, 'S11': 0, 'S01': 0, 'S10': 0})
    p_df['P'] = np.where(p_df['area'] > 0, -0.5 * (p_df['S00'] + p_df['S11']) / p_df['area'], 0.0)
    # Presion en unidades reducidas (P/g), para comparar directamente entre
    # configuraciones sin que la escala dependa de la aceleracion de la
    # gravedad usada en cada corrida.
    p_df['P_reducida'] = p_df['P'] / g

    # Tensor de esfuerzos por particula, con la MISMA convencion de signo que
    # P (compresion positiva): sigma_ij = -S_ij / area, asi P = (sxx+syy)/2.
    # La parte cruzada se simetriza: el tensor de una particula armado con
    # r x F no es exactamente simetrico si los torques no estan balanceados.
    area_ok = p_df['area'] > 0
    area_seg = p_df['area'].where(area_ok, 1.0)
    sxx = np.where(area_ok, -p_df['S00'] / area_seg, 0.0)
    syy = np.where(area_ok, -p_df['S11'] / area_seg, 0.0)
    sxy = np.where(area_ok, -0.5 * (p_df['S01'] + p_df['S10']) / area_seg, 0.0)
    # Esfuerzo de corte maximo (radio del circulo de Mohr, parte
    # deviatorica en 2D). No depende de la convencion de signo.
    p_df['tau'] = np.sqrt((0.5 * (sxx - syy)) ** 2 + sxy ** 2)
    p_df['tau_reducida'] = p_df['tau'] / g
    p_df['sigma_xx_reducida'] = sxx / g
    p_df['sigma_yy_reducida'] = syy / g
    p_df['sigma_xy_reducida'] = sxy / g
    
    # ROI: granos que no son del arco con el centro en la franja horizontal
    # [y_min, y_max + ROI_EXTRA] del arco (y_min, y_max = centros mas bajo y
    # mas alto del arco), de pared a pared. Con ROI_EXTRA = 0 (default) es la
    # misma franja que usa force_arch_v4.py para los contactos.
    y_arco = p_df.loc[p_df.index.isin(arco_path), 'pos_y'] if arco_path else []
    roi_ymin = y_arco.min() if len(y_arco) > 0 else 0
    roi_ymax = y_arco.max() + ROI_EXTRA if len(y_arco) > 0 else 20.0
    
    p_df['region'] = np.where(p_df.index.isin(arco_path), 'arco',
                     np.where((p_df['pos_y'] >= roi_ymin) & (p_df['pos_y'] <= roi_ymax), 'roi', 'bulk'))
    
    if 'tipo' not in p_df.columns:
        # Sin especie precalculada (process_one_jam_task la trae de la columna
        # 'type' del .xy, ver analysis/especies.py): umbral de radio.
        r_col = 'r_nominal' if 'r_nominal' in p_df.columns else 'r_eff'
        p_df['tipo'] = np.where(p_df[r_col] >= umbral_r, 'B', 'S')
    p_df['P_A'] = p_df['P'] * p_df['area']
    p_df['P_A_reducida'] = p_df['P_reducida'] * p_df['area']
    # Mismo esquema para el esfuerzo de corte: suma de tau*area, para poder
    # promediar pesando por area igual que P (media de |tau| por particula).
    p_df['T_A'] = p_df['tau'] * p_df['area']
    p_df['T_A_reducida'] = p_df['tau_reducida'] * p_df['area']

    res = {}
    for region in ['arco', 'roi', 'bulk']:
        mask_reg = p_df['region'] == region
        sub_df = p_df[mask_reg]
        
        # Guardamos numeradores y denominadores totales de la región
        res[f"SumP_{region}"] = sub_df['P_A'].sum()
        res[f"SumP_reducida_{region}"] = sub_df['P_A_reducida'].sum()
        res[f"SumTau_{region}"] = sub_df['T_A'].sum()
        res[f"SumTau_reducida_{region}"] = sub_df['T_A_reducida'].sum()
        res[f"Area_{region}"] = sub_df['area'].sum()
        
        for tipo in ['B', 'S']:
            mask_tipo = sub_df['tipo'] == tipo
            # Guardamos numeradores y denominadores por componente
            res[f"SumP_{region}_{tipo}"] = sub_df.loc[mask_tipo, 'P_A'].sum()
            res[f"SumP_reducida_{region}_{tipo}"] = sub_df.loc[mask_tipo, 'P_A_reducida'].sum()
            res[f"SumTau_{region}_{tipo}"] = sub_df.loc[mask_tipo, 'T_A'].sum()
            res[f"SumTau_reducida_{region}_{tipo}"] = sub_df.loc[mask_tipo, 'T_A_reducida'].sum()
            res[f"Area_{region}_{tipo}"] = sub_df.loc[mask_tipo, 'area'].sum()

    # Detalle por particula (para poder guardar la presion individual de
    # cada grano, ademas del resumen agregado por region de arriba).
    particulas = p_df[['pos_x', 'pos_y', 'shape', 'tipo', 'region', 'area', 'P', 'P_reducida',
                       'tau', 'tau_reducida', 'sigma_xx_reducida', 'sigma_yy_reducida',
                       'sigma_xy_reducida', 'coords']].reset_index()
    # 'coords' es None para discos, o un array (N,2) de vertices absolutos
    # para poligonos -> lo serializamos a JSON para que sobreviva el CSV.
    # Sin esto, graficar_presion_atasco.py no tiene forma de distinguir un
    # cuadrado/pentagono/triangulo de un circulo equivalente al dibujar.
    particulas['coords'] = particulas['coords'].apply(
        lambda c: json.dumps(np.asarray(c, dtype=float).round(5).tolist()) if c is not None else ''
    )

    del p_df
    return res, particulas

def _fijar_roi_extra(valor):
    # los procesos del Pool no siempre heredan lo que main() cambia: se pasa explicito
    global ROI_EXTRA
    ROI_EXTRA = float(valor)


def worker_wrapper(args):
    return process_one_jam_task(*args)

def process_one_jam_task(xy_file, contact_file, jam_id, config_id, dt, r_large_m, r_small_m, sides_large_m, sides_small_m, master_chi, sim_path):
    try:
        particles = parse_xy(xy_file)
        df_particles = pd.DataFrame.from_dict(particles, orient='index')
        df_particles.index.name = 'id'
        
        # Grande/chico por especie (columna 'type' del .xy), sin umbral de
        # radio: el inradio/circunradio de un poligono no se compara bien con
        # el radio de un disco. Ver analysis/especies.py.
        especie = especie_por_tipo(df_particles['type'], df_particles['area'], master_chi)
        df_particles['tipo'] = df_particles['type'].map(especie)
        hay_chicas = (df_particles['tipo'] == 'S').any()
        if not hay_chicas:
            r_s_val = np.nan
        elif pd.notna(r_small_m):
            r_s_val = r_small_m
        else:
            # Radio de igual area (la misma definicion que los radios de params.in)
            r_s_val = float(np.sqrt(df_particles.loc[df_particles['tipo'] == 'S', 'area'].median() / np.pi))
        umbral_r = 0.0   # no se usa: df_particles ya trae 'tipo'

        cols = ['id_a', 'id_b', 'px', 'py', 'f_norm', 'f_tan']
        dtypes_dict = {'id_a': np.int32, 'id_b': np.int32, 'px': np.float32, 'py': np.float32, 'f_norm': np.float32, 'f_tan': np.float32}
        df_contacts = pd.read_csv(contact_file, sep=r'\s+', names=cols, comment='#', dtype=dtypes_dict, engine='c')
        
        r_salida = parse_radio_salida(sim_path)
        g_abs = parse_g(sim_path)
        arco_path = detectar_arco_primario(df_particles, df_contacts, r_salida)
        if not arco_path:
            return {"ok": False, "error": "No arco detectado"}
            
        # Normal real de Box2D si el simulador escribio atasco_contacts_full;
        # si no (datos viejos), se infiere de la geometria.
        normales, fuente_normales = leer_normales_box2d(contact_file, df_contacts)
        summary, particulas = calcular_tensor_y_presion_vectorizado(df_particles, df_contacts, arco_path, dt, umbral_r, g=g_abs,
                                                                    paredes=parse_paredes(xy_file),
                                                                    normales=normales)
        summary['normales'] = fuente_normales
        
        if pd.notna(sides_large_m):
            forma_grande = nombre_forma(sides_large_m)
        else:
            grandes_shapes = df_particles.loc[df_particles['tipo'] == 'B', 'shape']
            forma_grande = grandes_shapes.mode()[0] if not grandes_shapes.empty else "Ninguna"

        if not hay_chicas:
            forma_chica = "Ninguna"
        elif pd.notna(sides_small_m):
            forma_chica = nombre_forma(sides_small_m)
        else:
            chicas_shapes = df_particles.loc[df_particles['tipo'] == 'S', 'shape']
            forma_chica = chicas_shapes.mode()[0] if not chicas_shapes.empty else "Ninguna"
        
        if pd.notna(master_chi):
            chi_val = float(master_chi)
        else:
            chi_match = CHI_REGEX.search(config_id)
            chi_val = float(chi_match.group(1).replace('p', '.')) if chi_match else np.nan

        summary['jam_id'] = jam_id
        summary['config_id'] = config_id
        summary['sim_path'] = sim_path 
        summary['forma_grande'] = forma_grande
        summary['forma_chica'] = forma_chica
        summary['chi'] = chi_val
        summary['r_chica'] = r_s_val 

        particulas['jam_id'] = jam_id
        particulas['config_id'] = config_id
        particulas['sim_path'] = sim_path
        particulas['chi'] = chi_val
        # Forma de la MEZCLA (grano grande / grano chico), igual en todas las
        # particulas del atasco: dividir_particulas_por_forma.py separa por
        # 'forma_grande', asi los discos chicos de una mezcla con triangulos
        # quedan en el archivo de Triangulo y no en el de Disco.
        particulas['forma_grande'] = forma_grande
        particulas['forma_chica'] = forma_chica
        
        del df_particles, df_contacts
        return {"ok": True, "summary": summary, "particulas": particulas}
    except Exception as e:
        return {"ok": False, "error": str(e)}

def parse_grain_column(val):
    if pd.isna(val):
        return {}
    val_str = str(val).strip()
    if val_str.startswith('"') and val_str.endswith('"'):
        val_str = val_str[1:-1]
    val_str = val_str.replace('""', '"')
    try:
        return json.loads(val_str)
    except:
        try:
            return ast.literal_eval(val_str)
        except:
            return {}


def cargar_master(path_master):
    df_master = pd.read_csv(path_master)
    df_master = df_master.drop_duplicates(subset=['path'])
    master_dict = {}
    for _, row in df_master.iterrows():
        path_val = row.get('path')
        if pd.isna(path_val):
            continue
        grain_data = parse_grain_column(row.get('grain_json'))
        large_radius = large_sides = small_radius = small_sides = np.nan
        if isinstance(grain_data, list) and len(grain_data) > 0:
            large_radius = grain_data[0].get('radius', np.nan)
            large_sides = grain_data[0].get('sides', np.nan)
            if len(grain_data) > 1:
                small_radius = grain_data[1].get('radius', np.nan)
                small_sides = grain_data[1].get('sides', np.nan)
        if pd.isna(large_radius) and pd.notna(row.get('r_large')):
            large_radius = row.get('r_large')
        if pd.isna(small_radius) and pd.notna(row.get('r_small')):
            small_radius = row.get('r_small')
        info = {'r_large': large_radius, 'r_small': small_radius,
                'sides_large': large_sides, 'sides_small': small_sides, 'chi': row.get('chi')}
        path_abs = os.path.abspath(os.path.normpath(str(path_val)))
        master_dict[path_abs] = info
        parts = path_abs.split(os.sep)
        if len(parts) >= 2:
            master_dict[os.path.join(parts[-2], parts[-1])] = info
    return master_dict


# ---------------------------------------------------------------------------
# Listado de archivos
# ---------------------------------------------------------------------------

def listar_atascos(dirs, incluir=None):
    """Una sola recorrida por el arbol (os.walk): junta los atasco_*.xy y
    sus atasco_contacts_*.dat con el listado de cada carpeta, sin un
    os.path.exists() por archivo (cada consulta al disco USB desde WSL es
    lenta)."""
    pares = []
    for d in dirs:
        for root, _, files in os.walk(d):
            nombres = set(files)
            for f in files:
                if f.startswith('atasco_') and f.endswith('.xy') and not f.startswith('atasco_contacts_'):
                    c = f.replace('atasco_', 'atasco_contacts_', 1)[:-3] + '.dat'
                    if c in nombres:
                        pares.append((os.path.join(root, f), os.path.join(root, c)))
    return _filtrar_y_ordenar(pares, incluir)


def leer_lista(path_lista, incluir=None):
    """Lee una lista de rutas (una por linea) con los atasco_*.xy y
    atasco_contacts_*.dat, por ejemplo generada desde Windows con
    'dir /s /b' (mucho mas rapido que recorrer el disco desde WSL).
    Acepta rutas de Windows (E:\\...) y las pasa a /mnt/e/..."""
    xys, contactos = [], set()
    with open(path_lista, encoding='utf-8', errors='replace') as f:
        for linea in f:
            p = linea.strip().replace('\r', '')
            if not p:
                continue
            m = re.match(r'^([A-Za-z]):[\\/](.*)$', p)
            if m:
                p = f"/mnt/{m.group(1).lower()}/" + m.group(2).replace('\\', '/')
            nombre = os.path.basename(p)
            if nombre.startswith('atasco_contacts_') and nombre.endswith('.dat'):
                contactos.add(p)
            elif nombre.startswith('atasco_') and nombre.endswith('.xy'):
                xys.append(p)
    pares = []
    for xy in xys:
        d, n = os.path.split(xy)
        c = os.path.join(d, n.replace('atasco_', 'atasco_contacts_', 1)[:-3] + '.dat')
        if c in contactos:
            pares.append((xy, c))
    return _filtrar_y_ordenar(pares, incluir)


def _filtrar_y_ordenar(pares, incluir):
    if incluir:
        rx = re.compile(incluir)
        pares = [p for p in pares if rx.search(p[0])]
    # orden por ruta: los archivos de una misma simulacion se leen juntos
    # (en un disco rotativo eso reduce saltos del cabezal)
    return sorted(pares)


def generar_tareas(pares, master_dict, dt):
    for xy_file, contact_file in pares:
        jam_id = os.path.basename(xy_file).split('_')[-1].replace('.xy', '')
        config_id = os.path.basename(os.path.dirname(xy_file))
        abs_xy = os.path.abspath(os.path.normpath(xy_file))
        sim_path = os.path.dirname(os.path.dirname(abs_xy))  # raiz de la simulacion
        config_info = master_dict.get(sim_path)
        if not config_info:
            parts = sim_path.split(os.sep)
            if len(parts) >= 2:
                config_info = master_dict.get(os.path.join(parts[-2], parts[-1]), {})
            else:
                config_info = {}
        yield (xy_file, contact_file, jam_id, config_id, dt,
               config_info.get('r_large', np.nan), config_info.get('r_small', np.nan),
               config_info.get('sides_large', np.nan), config_info.get('sides_small', np.nan),
               config_info.get('chi', np.nan), sim_path)


# ---------------------------------------------------------------------------
# Escritura por lotes, reanudable
# ---------------------------------------------------------------------------
#
# Cada lote (unos cientos de miles de filas de particulas) se guarda asi:
#   1. out/partes/<forma>__<lote>.parquet   (un archivo por forma_grande)
#   2. filas del lote en resumen_presion_por_atasco.csv (columna 'lote')
#   3. los .xy del lote en procesados.txt ("<lote>\t<ruta>")
#   4. el id del lote en lotes_ok.txt   <- recien aca el lote cuenta
# Si la corrida se corta en cualquier punto, al reanudar se descartan los
# restos de lotes que no llegaron al paso 4 y esos atascos se reprocesan:
# no se pierde ni se duplica nada.

META_PART = ['forma_grande', 'forma_chica', 'config_id', 'sim_path', 'jam_id', 'chi', 'id']
META_ATASCO = ['forma_grande', 'forma_chica', 'chi', 'r_chica', 'config_id', 'sim_path', 'jam_id']


class Registro:
    def __init__(self, out_dir, reanudar):
        self.out = out_dir
        self.dir_partes = os.path.join(out_dir, 'partes')
        self.f_ok = os.path.join(out_dir, 'lotes_ok.txt')
        self.f_proc = os.path.join(out_dir, 'procesados.txt')
        self.f_resumen = os.path.join(out_dir, 'resumen_presion_por_atasco.csv')
        os.makedirs(self.dir_partes, exist_ok=True)
        if not reanudar:
            for f in (self.f_ok, self.f_proc, self.f_resumen):
                if os.path.exists(f):
                    os.remove(f)
            for f in os.listdir(self.dir_partes):
                os.remove(os.path.join(self.dir_partes, f))
        self.lotes_ok = self._leer_ok()
        # restos de lotes sin confirmar (corrida anterior cortada a la mitad)
        for f in os.listdir(self.dir_partes):
            lote = f.rsplit('__', 1)[-1].replace('.parquet', '')
            if lote not in self.lotes_ok:
                os.remove(os.path.join(self.dir_partes, f))
        self.run = 1 + max([int(l.split('_')[0]) for l in self.lotes_ok] or [0])
        self.n_lote = 0

    def _leer_ok(self):
        if not os.path.exists(self.f_ok):
            return set()
        with open(self.f_ok) as f:
            return {l.strip() for l in f if l.strip()}

    def ya_procesados(self):
        hechos = set()
        if os.path.exists(self.f_proc):
            with open(self.f_proc) as f:
                for l in f:
                    lote, _, ruta = l.rstrip('\n').partition('\t')
                    if lote in self.lotes_ok:
                        hechos.add(ruta)
        return hechos

    def guardar_lote(self, particulas, resumenes, xy_files):
        import pyarrow as pa
        import pyarrow.parquet as pq
        self.n_lote += 1
        lote = f"{self.run:03d}_{self.n_lote:05d}"
        if particulas:
            df = pd.concat(particulas, ignore_index=True)
            df = df[[c for c in META_PART if c in df.columns] +
                    [c for c in df.columns if c not in META_PART]]
            for forma, sub in df.groupby('forma_grande', sort=False):
                tmp = os.path.join(self.dir_partes, f".tmp_{forma}__{lote}.parquet")
                pq.write_table(pa.Table.from_pandas(sub, preserve_index=False), tmp, compression='zstd')
                os.replace(tmp, os.path.join(self.dir_partes, f"{forma}__{lote}.parquet"))
            del df
        if resumenes:
            dr = pd.DataFrame(resumenes)
            dr['lote'] = lote
            dr = dr[[c for c in META_ATASCO if c in dr.columns] +
                    [c for c in dr.columns if c not in META_ATASCO]]
            dr.to_csv(self.f_resumen, mode='a', index=False, header=not os.path.exists(self.f_resumen))
        with open(self.f_proc, 'a') as f:
            f.writelines(f"{lote}\t{x}\n" for x in xy_files)
        with open(self.f_ok, 'a') as f:
            f.write(lote + "\n")
            f.flush()
            os.fsync(f.fileno())
        self.lotes_ok.add(lote)


def combinar_partes(out_dir, borrar=False):
    """Junta out/partes/<forma>__*.parquet en un
    presion_por_particula_<forma>.parquet por forma (el formato que lee
    barrido_heatmap.py). Lee y escribe de a un archivo parcial por vez."""
    import pyarrow.parquet as pq
    dir_partes = os.path.join(out_dir, 'partes')
    por_forma = collections.defaultdict(list)
    for f in sorted(os.listdir(dir_partes)):
        if f.endswith('.parquet') and not f.startswith('.tmp_'):
            por_forma[f.rsplit('__', 1)[0]].append(os.path.join(dir_partes, f))
    for forma, archivos in por_forma.items():
        destino = os.path.join(out_dir, f"presion_por_particula_{forma}.parquet")
        tmp = destino + ".tmp"
        escritor, esquema, filas = None, None, 0
        for a in archivos:
            t = pq.read_table(a)
            if escritor is None:
                esquema = t.schema
                escritor = pq.ParquetWriter(tmp, esquema, compression='zstd')
            elif not t.schema.equals(esquema):
                t = t.cast(esquema)
            escritor.write_table(t)
            filas += t.num_rows
        escritor.close()
        os.replace(tmp, destino)
        print(f"  {os.path.basename(destino)}: {filas:,} filas ({len(archivos)} partes)", flush=True)
        if borrar:
            for a in archivos:
                os.remove(a)


def resumen_por_configuracion(out_dir, lotes_ok):
    atasco_csv_path = os.path.join(out_dir, "resumen_presion_por_atasco.csv")
    config_csv_path = os.path.join(out_dir, "resumen_presion_por_configuracion.csv")
    if not os.path.exists(atasco_csv_path):
        return
    df_all = pd.read_csv(atasco_csv_path, dtype={'lote': str, 'jam_id': str})
    if 'lote' in df_all.columns:
        df_all = df_all[df_all['lote'].isin(lotes_ok)]
    group_cols = ['forma_grande', 'forma_chica', 'chi', 'r_chica', 'config_id', 'sim_path']
    grouped = df_all.groupby(group_cols, dropna=False)
    df_config = grouped.size().reset_index(name='n_atascos')
    variables_sumar = [c for c in df_all.columns
                       if c.startswith('SumP_') or c.startswith('SumTau_') or c.startswith('Area_')]
    cols_to_agg = variables_sumar
    if not cols_to_agg:
        return
    df_sums = grouped[cols_to_agg].sum().reset_index(drop=True)
    df_config = pd.concat([df_config, df_sums], axis=1)
    # mean_<magnitud>_<region>[_B|_S] = Sum<magnitud>_<region>[_B|_S] / Area_<region>[_B|_S]
    for c in cols_to_agg:
        for pref, nombre in (('SumP_', 'mean_P_'), ('SumTau_', 'mean_tau_')):
            if c.startswith(pref):
                resto = c[len(pref):]                    # ej. 'reducida_arco_B' o 'arco_B'
                area = 'Area_' + resto.replace('reducida_', '', 1)
                if area in df_config.columns:
                    df_config[nombre + resto] = df_config[c] / df_config[area]
    df_config = df_config.drop(columns=cols_to_agg)
    df_config.to_csv(config_csv_path, index=False)
    print(f"Resumen por configuracion: {config_csv_path}", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Tensor de esfuerzos, presion y corte por particula en cada atasco.")
    parser.add_argument("--dir", nargs='+', default=None,
                        help="carpeta(s) raiz de las simulaciones (se recorren recursivamente)")
    parser.add_argument("--lista", default=None,
                        help="archivo con la lista de atasco_*.xy y atasco_contacts_*.dat (una ruta "
                             "por linea, sirve la salida de 'dir /s /b' de Windows); evita recorrer "
                             "el disco desde WSL")
    parser.add_argument("--incluir", default=None,
                        help="regex sobre la ruta: procesa solo los atascos que coinciden "
                             "(ej. 'Ratio_0p1_Chi_0p(0|2)_W')")
    parser.add_argument("--out", required=True)
    parser.add_argument("--master", required=True)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--roi_extra", type=float, default=ROI_EXTRA,
                        help="altura extra de la ROI por encima del grano mas alto del arco "
                             "(unidades de simulacion; default 0: solo la franja del arco). Antes era 10.")
    parser.add_argument("--filas_lote", type=int, default=500_000,
                        help="filas de particulas por lote guardado (default 500000)")
    parser.add_argument("--reanudar", action="store_true",
                        help="sigue una corrida cortada: saltea los atascos ya guardados en --out")
    parser.add_argument("--solo_combinar", action="store_true",
                        help="no procesa nada: solo junta out/partes/ en un parquet por forma")
    parser.add_argument("--borrar_partes", action="store_true",
                        help="borra out/partes/ despues de combinar")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.solo_combinar:
        print("Combinando partes...", flush=True)
        combinar_partes(args.out, borrar=args.borrar_partes)
        return

    if not args.dir and not args.lista:
        raise SystemExit("Hace falta --dir o --lista.")

    registro = Registro(args.out, args.reanudar)

    print("Cargando master index y parseando JSON de granos...", flush=True)
    master_dict = cargar_master(args.master)

    t_busq = time.time()
    if args.lista:
        print(f"Leyendo lista de archivos {args.lista}...", flush=True)
        pares = leer_lista(args.lista, args.incluir)
    else:
        print("Recorriendo carpetas en busca de atasco_*.xy...", flush=True)
        pares = listar_atascos(args.dir, args.incluir)
    print(f"Atascos encontrados: {len(pares)}  (busqueda: {time.time() - t_busq:.0f} s)", flush=True)

    if args.reanudar:
        hechos = registro.ya_procesados()
        antes = len(pares)
        pares = [p for p in pares if p[0] not in hechos]
        print(f"--reanudar: {antes - len(pares)} ya procesados, quedan {len(pares)}", flush=True)

    total = len(pares)
    part_buf, res_buf, xy_buf = [], [], []
    filas_buf = 0
    procesados = 0
    errores = collections.Counter()
    t0 = time.time()

    def guardar():
        nonlocal part_buf, res_buf, xy_buf, filas_buf
        registro.guardar_lote(part_buf, res_buf, xy_buf)
        part_buf, res_buf, xy_buf = [], [], []
        filas_buf = 0
        gc.collect()

    tareas = generar_tareas(pares, master_dict, args.dt)
    _fijar_roi_extra(args.roi_extra)
    print(f"ROI: franja del arco + {args.roi_extra:g} por encima del grano mas alto", flush=True)
    with mp.Pool(processes=args.workers, maxtasksperchild=50,
                 initializer=_fijar_roi_extra, initargs=(args.roi_extra,)) as pool:
        for xy, res in pool.imap_unordered(worker_wrapper_con_ruta, tareas, chunksize=4):
            procesados += 1
            xy_buf.append(xy)   # tambien los fallidos: no tiene sentido reintentarlos
            if res["ok"]:
                res_buf.append(res["summary"])
                part_buf.append(res["particulas"])
                filas_buf += len(res["particulas"])
            else:
                errores[res.get("error", "?")] += 1
            if filas_buf >= args.filas_lote:
                guardar()
            if procesados % 200 == 0:
                dtr = time.time() - t0
                ritmo = procesados / dtr
                resta = (total - procesados) / ritmo if ritmo > 0 else float('nan')
                print(f"Procesados {procesados} de {total}  ({ritmo:.1f} atascos/s, "
                      f"faltan ~{resta / 60:.0f} min, {sum(errores.values())} fallidos)", flush=True)
    if xy_buf:
        guardar()

    print(f"\n{procesados} atascos en {(time.time() - t0) / 60:.1f} min, "
          f"{sum(errores.values())} fallidos.", flush=True)
    for msg, n in errores.most_common(10):
        print(f"  {n:>6}  {msg}")

    print("\nCombinando partes en un parquet por forma_grande...", flush=True)
    combinar_partes(args.out, borrar=args.borrar_partes)
    resumen_por_configuracion(args.out, registro.lotes_ok)


def worker_wrapper_con_ruta(args):
    return args[0], process_one_jam_task(*args)


if __name__ == "__main__":
    mp.freeze_support()
    main()