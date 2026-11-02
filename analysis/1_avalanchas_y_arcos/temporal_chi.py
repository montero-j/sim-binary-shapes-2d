import argparse
import json
import os
from collections import Counter
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # backend sin display; evita fallos silenciosos de
                        # plt.figure()/savefig() bajo multiprocessing en
                        # entornos headless (server/cluster sin X11).
import matplotlib.pyplot as plt

# =============================================================================
# ESTILO -- comun a todos los scripts (estilo_figuras.py, base: graficos_finales.py)
# =============================================================================
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # estilo_figuras.py y rutas.py viven en analysis/
import rutas
from estilo_figuras import aplicar_estilo, estilo_sistema, guardar_figura, REF_LINE_KW

# Letras mas grandes y sin grilla que el estilo comun (estas figuras son
# anchas y van a ancho de columna)
FUENTES = {'axes.labelsize': 26, 'xtick.labelsize': 22, 'ytick.labelsize': 22,
           'legend.fontsize': 18, 'legend.title_fontsize': 18}


def _estilo():
    aplicar_estilo()
    plt.rcParams.update(FUENTES)
    plt.rcParams['axes.grid'] = False


# a nivel de modulo, asi los procesos del Pool tambien lo tienen
_estilo()

# Colores para la figura comparativa (paleta Wong, la misma de estilo_figuras)
COLORES_CASOS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#000000']
SUAVIZADO_BATCHES = 5   # media movil (en ventanas) para medir segregacion

# =============================================================================
# CONFIGURACIÓN (valores por defecto; se cambian por argumentos)
# =============================================================================
INDEX_FILE = rutas.MASTER_INDEX
WINDOW_SIZE = 10.0   # tamaño del batch en segundos
CPU_CORES = cpu_count()
FORMATO = 'pdf'      # igual que save_plot() de graficos_finales.py
# Reemplazo de prefijo en la columna 'path' del index (ej. el index tiene
# /home/juli/work/... pero en esta PC el disco esta en /mnt/h/work/...).
PREFIJO_VIEJO = None
PREFIJO_NUEVO = None
# Intervalo de tiempo a analizar (s); None = desde el principio / hasta el final
T_MIN = None
T_MAX = None


def area_grano(radio, n_lados):
    """Area de un grano de radio r: disco (n_lados = 1) o poligono regular
    de n_lados con los vertices sobre la circunferencia de radio r."""
    n = int(n_lados) if n_lados else 1
    if n <= 2:
        return np.pi * radio ** 2
    return 0.5 * n * radio ** 2 * np.sin(2 * np.pi / n)


def calcular_chi_temporal(row):
    """
    chi de salida (fraccion de area de chicos) por ventana de WINDOW_SIZE s.
    Devuelve (info, df) con info['status'] = 'ok' o el motivo exacto por el
    que no se pudo calcular; df tiene columnas tiempo, chi_eff (o None).
    """
    path = row['path']
    if PREFIJO_VIEJO and str(path).startswith(PREFIJO_VIEJO):
        path = PREFIJO_NUEVO + str(path)[len(PREFIJO_VIEJO):]
    flux_filename = row['fluxFile']
    sim_hash = row['sim_hash']
    tirada = row['tirada']
    base = {'sim_hash': sim_hash, 'tirada': tirada, 'path': path}

    if not os.path.isdir(path):
        return {**base, 'status': 'dir_no_existe'}, None

    try:
        grain_info = json.loads(row['grain_json'])
    except Exception as e:
        return {**base, 'status': 'json_invalido', 'detail': str(e)}, None

    # chi_nominal viene DIRECTO del index (columna 'chi'), no se recalcula
    # desde grain_json. grain_json solo trae radius/sides, nunca tuvo
    # counts de población -- por eso el KeyError('count') anterior.
    # 'chi' ya es la concentración objetivo del diseño experimental
    # (confirmado: a r_ratio_nominal fijo, 'chi' es la variable barrida).
    chi_nominal = row['chi']

    if len(grain_info) < 2:
        # Sistema monodisperso (sin aditivo): no hay dos poblaciones que
        # comparar en el tiempo, se omite.
        return {**base, 'status': 'monodisperso_omitido'}, None

    radios = [g['radius'] for g in grain_info]
    # Area real de cada especie: los poligonos tienen los n vertices sobre la
    # circunferencia de radio r (globalSetup.cpp), area (n/2) r^2 sin(2pi/n);
    # los discos pi r^2. Antes se usaba pi r^2 para todos, lo que sesgaba
    # chi_out en mezclas de formas distintas (ej. cuadrado-disco).
    areas_unitarias = {i: area_grano(g['radius'], g.get('sides', 1)) for i, g in enumerate(grain_info)}
    tipo_chica = int(np.argmin(radios))

    flux_path = os.path.join(path, flux_filename)
    if not os.path.exists(flux_path):
        return {**base, 'status': 'flux_no_existe', 'detail': flux_path}, None

    try:
        df_flux = pd.read_csv(flux_path, sep=r'\s+', comment='#', header=None, usecols=[1, 2])
    except Exception as e:
        return {**base, 'status': 'error_lectura_flux', 'detail': str(e)}, None

    if df_flux.empty:
        return {**base, 'status': 'flux_vacio'}, None
    df_flux.columns = ['tipo', 'tiempo']
    # Lineas incompletas o corruptas (ej. la ultima, cortada al interrumpir la
    # simulacion) se leen como NaN: se descartan y se informa cuantas eran.
    df_flux = df_flux.apply(pd.to_numeric, errors='coerce')
    n_total = len(df_flux)
    df_flux = df_flux.dropna()
    lineas_descartadas = n_total - len(df_flux)
    if df_flux.empty:
        return {**base, 'status': 'flux_vacio', 'detail': f'{lineas_descartadas} lineas invalidas'}, None
    df_flux['tipo'] = df_flux['tipo'].astype(int)
    if lineas_descartadas:
        base['lineas_descartadas'] = lineas_descartadas

    # Chequeo de consistencia: los 'tipo' de flx.dat deben ser subset de
    # los índices 0..N-1 derivados de grain_json.
    tipos_en_flux = set(df_flux['tipo'].unique())
    tipos_esperados = set(areas_unitarias.keys())
    if not tipos_en_flux.issubset(tipos_esperados):
        return {**base, 'status': 'tipo_mismatch',
                'detail': f'flux tiene {sorted(tipos_en_flux)}, se esperaba subset de {sorted(tipos_esperados)}'}, None

    # recorte temporal: las ventanas arrancan en T_MIN (si se dio)
    if T_MIN is not None:
        df_flux = df_flux[df_flux['tiempo'] >= T_MIN]
    if T_MAX is not None:
        df_flux = df_flux[df_flux['tiempo'] <= T_MAX]
    if df_flux.empty:
        return {**base, 'status': 'sin_datos_en_intervalo', 'detail': f'[{T_MIN}, {T_MAX}] s'}, None
    df_flux = df_flux.copy()

    t_start = T_MIN if T_MIN is not None else df_flux['tiempo'].min()
    t_end = df_flux['tiempo'].max()
    bins = np.arange(t_start, t_end + WINDOW_SIZE, WINDOW_SIZE)
    if len(bins) < 2:
        return {**base, 'status': 'pocos_bins', 'detail': f't_start={t_start}, t_end={t_end}'}, None

    df_flux['batch'] = pd.cut(df_flux['tiempo'], bins=bins, labels=bins[:-1])

    batch_results = []
    for label, group in df_flux.groupby('batch', observed=True):
        counts = group['tipo'].value_counts()
        a_chicas_out = counts.get(tipo_chica, 0) * areas_unitarias[tipo_chica]
        a_total_out = sum(counts.get(t, 0) * areas_unitarias[t] for t in areas_unitarias.keys())
        if a_total_out > 0:
            batch_results.append({'tiempo': float(label) + WINDOW_SIZE / 2,
                                  'chi_eff': a_chicas_out / a_total_out})

    df_plot = pd.DataFrame(batch_results)
    if df_plot.empty:
        return {**base, 'status': 'sin_batches_validos'}, None
    return {**base, 'status': 'ok', 'chi_nominal': chi_nominal}, df_plot


def medir_segregacion(df_plot, chi_nominal):
    """Cuanto se aparta chi_out del nominal, sobre la curva suavizada con una
    media movil de SUAVIZADO_BATCHES ventanas (asi no cuenta el ruido de
    cada ventana). Relativo a chi_nominal:
      desvio_medio_rel  media de |chi_suave - chi_nom| / chi_nom
      desvio_max_rel    maximo de lo mismo
      deriva_rel        (ultimo 20 % - primer 20 % del tiempo) / chi_nom
    """
    chi = df_plot['chi_eff'].to_numpy()
    suave = pd.Series(chi).rolling(SUAVIZADO_BATCHES, center=True, min_periods=1).mean().to_numpy()
    n = max(1, len(chi) // 5)
    dev = np.abs(suave - chi_nominal) / chi_nominal
    return {'desvio_medio_rel': float(dev.mean()), 'desvio_max_rel': float(dev.max()),
            'deriva_rel': float((chi[-n:].mean() - chi[:n].mean()) / chi_nominal),
            'chi_out_medio': float(chi.mean()), 'n_ventanas': len(chi)}


def analizar_flujo_temporal(row):
    """
    Grafico chi_out(t) de una simulacion (en su carpeta). Devuelve un dict
    con 'status' explicando el resultado (ok, o el motivo exacto de por
    que no se genero grafico para esta fila) y las medidas de segregacion.
    """
    info, df_plot = calcular_chi_temporal(row)
    if df_plot is None:
        return info
    base = {k: info[k] for k in ('sim_hash', 'tirada', 'path', 'lineas_descartadas') if k in info}
    path, sim_hash, tirada, chi_nominal = info['path'], info['sim_hash'], info['tirada'], info['chi_nominal']
    if chi_nominal > 0:
        base.update(medir_segregacion(df_plot, chi_nominal))

    try:
        color, marker, label_en = estilo_sistema(row.get('system_type', ''))

        plt.figure(figsize=(10, 5))
        # referencia: linea negra discontinua, como el chi = 0 de graficos_finales.py
        plt.axhline(y=chi_nominal, **REF_LINE_KW,
                    label=rf'$\chi = {chi_nominal:.3g}$')
        # solo linea (sin simbolos)
        plt.plot(df_plot['tiempo'], df_plot['chi_eff'], '-', color=color, lw=1.2,
                 label=rf'$\chi_{{out}}$ ({WINDOW_SIZE:g} s windows)', zorder=2)

        # diferencia relativa entre el ultimo batch y el valor nominal
        last_row = df_plot.iloc[-1]
        t_final, chi_final = last_row['tiempo'], last_row['chi_eff']
        plt.plot([t_final, t_final], [chi_final, chi_nominal], color='gray', linestyle=':',
                 linewidth=1.5, zorder=1)
        if chi_nominal > 0:
            diff_pct = (chi_nominal - chi_final) / chi_nominal * 100
            plt.annotate(f"{diff_pct:.1f}%", xy=(t_final, (chi_final + chi_nominal) / 2),
                         xytext=(6, 0), textcoords='offset points', va='center',
                         color='gray', fontsize=14)

        plt.grid(False)
        plt.xlabel(r"Simulation Time (s)")
        plt.ylabel(r"Outflow Concentration $\chi_{out}$")
        plt.ylim(0, min(1.0, max(chi_nominal * 2.5, df_plot['chi_eff'].max() * 1.3, 0.05)))
        _xlim_intervalo(plt.gca())
        plt.legend(title=f"{label_en}, run {tirada}", loc='upper right')

        graph_name = f"dinamica_flujo_final_{sim_hash}_T{tirada}.{FORMATO}"
        plt.savefig(os.path.join(path, graph_name), bbox_inches='tight')
        plt.close()
    except Exception as e:
        plt.close('all')
        return {**base, 'status': 'error_al_graficar', 'detail': str(e)}

    return {**base, 'status': 'ok'}


def _xlim_intervalo(ax):
    """Si se pidio un intervalo, el eje x va exactamente de T_MIN a T_MAX."""
    if T_MIN is not None or T_MAX is not None:
        lo, hi = ax.get_xlim()
        ax.set_xlim(T_MIN if T_MIN is not None else lo, T_MAX if T_MAX is not None else hi)


def figura_comparativa(df_idx, casos, out, normalizar=False):
    """chi_out(t) de varias simulaciones en un mismo grafico (solo lineas),
    cada una con su chi nominal como linea discontinua del mismo color.
    casos: lista de sim_hash (se aceptan prefijos)."""
    _estilo()
    fig, ax = plt.subplots(figsize=(10, 5))
    y_max = 0.0
    for k, h in enumerate(casos):
        fila = df_idx[df_idx['sim_hash'].astype(str).str.startswith(str(h))]
        if fila.empty:
            print(f"  [!] sim_hash {h} no esta en el index")
            continue
        row = fila.iloc[0]
        info, d = calcular_chi_temporal(row)
        if d is None:
            print(f"  [!] {h}: {info['status']} {info.get('detail', '')}")
            continue
        chi_nom = info['chi_nominal']
        y = d['chi_eff'] / chi_nom if normalizar else d['chi_eff']
        color = COLORES_CASOS[k % len(COLORES_CASOS)]
        _, _, nombre = estilo_sistema(row.get('system_type', ''))
        extra = ''
        if 'area_ratio_nominal' in row and pd.notnull(row['area_ratio_nominal']):
            extra = rf', $A_r = {row["area_ratio_nominal"]:g}$'
        ax.plot(d['tiempo'], y, '-', color=color, lw=1.2, zorder=2,
                label=rf'{nombre}, $\chi = {chi_nom:.3g}${extra}')
        ax.axhline(1.0 if normalizar else chi_nom, color=color, ls='--', lw=1.5, zorder=3)
        y_max = max(y_max, float(y.max()))
        if chi_nom > 0:
            seg = medir_segregacion(d, chi_nom)
            print(f"  {row['sim_hash']}: {nombre} chi={chi_nom:g}  desvio medio {100 * seg['desvio_medio_rel']:.1f} %"
                  f"  max {100 * seg['desvio_max_rel']:.1f} %  deriva {100 * seg['deriva_rel']:+.1f} %")
    if not ax.lines:
        plt.close(fig)
        print("[!] Ningun caso tiene datos: no se guarda la figura.")
        return
    ax.grid(False)
    ax.set_xlabel(r'Simulation Time (s)')
    ax.set_ylabel(r'$\chi_{out} / \chi$' if normalizar else r'Outflow Concentration $\chi_{out}$')
    # margen arriba para que la leyenda (adentro) no tape las curvas
    n_casos = len(ax.get_legend_handles_labels()[1])
    ax.set_ylim(0, y_max * (1.2 + 0.12 * n_casos) if y_max > 0 else 1)
    _xlim_intervalo(ax)
    ax.legend(loc='upper right')
    guardar_figura(out, fig)
    plt.close(fig)
    print(f"Figura: {out if os.path.splitext(out)[1] else out + '.' + FORMATO}")


def _inicializar_worker(viejo, nuevo, ventana, formato, t_min=None, t_max=None):
    # los workers no heredan lo que main() cambia despues de importar el
    # modulo (con 'spawn' seguro, con 'fork' depende): se pasa explicito
    global PREFIJO_VIEJO, PREFIJO_NUEVO, WINDOW_SIZE, FORMATO, T_MIN, T_MAX
    PREFIJO_VIEJO, PREFIJO_NUEVO, WINDOW_SIZE, FORMATO = viejo, nuevo, ventana, formato
    T_MIN, T_MAX = t_min, t_max


def _segregacion_fila(row):
    info, d = calcular_chi_temporal(row)
    if d is None or not info.get('chi_nominal'):
        return None
    meta = {k: row[k] for k in ('system_type', 'chi', 'area_ratio_nominal', 'w_salida') if k in row}
    return {**{k: info[k] for k in ('sim_hash', 'tirada', 'path', 'lineas_descartadas') if k in info}, **meta,
            **medir_segregacion(d, info['chi_nominal'])}


def main():
    parser = argparse.ArgumentParser(
        description="Grafica chi de salida en funcion del tiempo para cada simulacion del index.")
    parser.add_argument('--index', default=INDEX_FILE, help='master_sim_index.csv')
    parser.add_argument('--dir', nargs='+', default=None,
                        help="carpeta(s) de simulaciones a recorrer (ej. .../R_2.0/CUADRADOS .../TRIANGULOS); "
                             "arma el indice en el momento con indexer.py y no usa --index")
    parser.add_argument('--incluir', default=None,
                        help="regex sobre la ruta de la simulacion: solo las que coinciden")
    parser.add_argument('--out_csv', default=None,
                        help="csv del --ranking (default: 5_graficos_finales/chi_temporal_ranking_segregacion.csv)")
    parser.add_argument('--top', type=int, default=20, help="cuantas mostrar por pantalla en --ranking")
    parser.add_argument('--ventana', type=float, default=WINDOW_SIZE,
                        help='tamaño de la ventana temporal en segundos (default 10)')
    parser.add_argument('--formato', choices=['pdf', 'svg', 'png'], default=FORMATO,
                        help='formato de salida (default pdf, como graficos_finales.py)')
    parser.add_argument('--prefijo_viejo', default=None,
                        help="prefijo de las rutas del index a reemplazar (ej. /home/juli/work)")
    parser.add_argument('--prefijo_nuevo', default=None,
                        help="prefijo nuevo (ej. /mnt/h/work)")
    parser.add_argument('--workers', type=int, default=CPU_CORES)
    parser.add_argument('--casos', nargs='+', default=None,
                        help="sim_hash (o prefijos) para UNA figura comparativa con todas las curvas; "
                             "no genera las figuras por simulacion")
    parser.add_argument('--normalizar', action='store_true',
                        help="figura comparativa con chi_out / chi nominal (referencia en 1)")
    parser.add_argument('--t_min', type=float, default=None,
                        help="inicio del intervalo de tiempo (s); vale para todos los modos")
    parser.add_argument('--t_max', type=float, default=None,
                        help="fin del intervalo de tiempo (s); vale para todos los modos")
    parser.add_argument('--out', default=os.path.join(rutas.GRAFICOS_FINALES_DIR, 'chi_temporal_comparacion'),
                        help="figura comparativa (sin extension -> --formato)")
    parser.add_argument('--ranking', action='store_true',
                        help="no grafica: calcula la segregacion de cada simulacion y lista las que mas segregan")
    parser.add_argument('--sistema', nargs='+', default=None,
                        help="filtra system_type (ej. Disco-Disco Cuad-Disco)")
    parser.add_argument('--chi', type=float, nargs='+', default=None, help="filtra chi nominal")
    parser.add_argument('--area_ratio', type=float, nargs='+', default=None,
                        help="filtra area_ratio_nominal (relacion de areas chico/grande, ej. 0.1)")
    args = parser.parse_args()
    if (args.prefijo_viejo is None) != (args.prefijo_nuevo is None):
        raise SystemExit("--prefijo_viejo y --prefijo_nuevo van juntos.")

    if args.dir:
        from indexer import indexar
        partes = []
        for d in args.dir:
            print(f"Indexando {d} ...")
            parte = indexar(d)
            n = 0 if parte is None else len(parte)
            print(f"  {n} simulaciones encontradas.")
            if parte is not None:
                partes.append(parte)
        if not partes:
            raise SystemExit(f"No encontre simulaciones (carpetas con params.in) en {args.dir}")
        df_idx = pd.concat(partes, ignore_index=True).drop_duplicates(subset=['path'])
    else:
        if not os.path.exists(args.index):
            print(f"Error: No se encontró {args.index} para chi_temporal.")
            return
        df_idx = pd.read_csv(args.index)
    if args.incluir:
        df_idx = df_idx[df_idx['path'].astype(str).str.contains(args.incluir, regex=True)]
        print(f"--incluir {args.incluir!r}: quedan {len(df_idx)} simulaciones.")
        if df_idx.empty:
            raise SystemExit("Ninguna simulacion coincide con --incluir.")
    # este proceso tambien usa las opciones (figura comparativa)
    _inicializar_worker(args.prefijo_viejo, args.prefijo_nuevo, args.ventana, args.formato, args.t_min, args.t_max)
    if args.sistema:
        df_idx = df_idx[df_idx['system_type'].isin(args.sistema)]
    if args.chi:
        df_idx = df_idx[np.isclose(df_idx['chi'].to_numpy()[:, None], args.chi).any(axis=1)]
    if args.area_ratio:
        if 'area_ratio_nominal' not in df_idx.columns:
            raise SystemExit("El index no tiene area_ratio_nominal: volve a correr indexer.py.")
        df_idx = df_idx[np.isclose(df_idx['area_ratio_nominal'].to_numpy(float)[:, None],
                                   args.area_ratio).any(axis=1)]
    if args.sistema or args.chi or args.area_ratio:
        print(f"Filtros: quedan {len(df_idx)} simulaciones.")
        if df_idx.empty:
            raise SystemExit("Ninguna simulacion cumple los filtros (revisar --sistema / --chi / --area_ratio).")

    if args.casos:
        figura_comparativa(df_idx, args.casos, args.out, normalizar=args.normalizar)
        return

    if args.ranking:
        tasks = [row for _, row in df_idx.iterrows()]
        with Pool(args.workers, initializer=_inicializar_worker,
                  initargs=(args.prefijo_viejo, args.prefijo_nuevo, args.ventana, args.formato,
                            args.t_min, args.t_max)) as p:
            res = p.map(_segregacion_fila, tasks)
        r = pd.DataFrame([x for x in res if x is not None])
        if r.empty:
            print("Ninguna simulacion se pudo calcular (revisar rutas / --prefijo_viejo).")
            return
        r = r.sort_values('desvio_medio_rel', ascending=False)
        out_csv = args.out_csv or os.path.join(rutas.GRAFICOS_FINALES_DIR, 'chi_temporal_ranking_segregacion.csv')
        os.makedirs(os.path.dirname(out_csv), exist_ok=True)
        r.to_csv(out_csv, index=False)
        cols = [c for c in ('sim_hash', 'system_type', 'chi', 'area_ratio_nominal', 'tirada',
                            'desvio_medio_rel', 'desvio_max_rel', 'deriva_rel') if c in r.columns]
        print(f"Segregacion de chi_out respecto del nominal (media movil de {SUAVIZADO_BATCHES} ventanas), "
              f"las {args.top} que mas segregan de {len(r)}:")
        with pd.option_context('display.width', 200):
            print(r[cols].head(args.top).to_string(index=False, float_format='%.3f'))
        print(f"\nTabla completa: {out_csv}")
        print("Mediana de desvio_medio_rel: " f"{r['desvio_medio_rel'].median():.3f}")
        return

    print(f"Analizando flujos temporales (Chi) sobre {len(df_idx)} filas...")

    tasks = [row for _, row in df_idx.iterrows()]
    with Pool(args.workers, initializer=_inicializar_worker,
              initargs=(args.prefijo_viejo, args.prefijo_nuevo, args.ventana, args.formato,
                            args.t_min, args.t_max)) as p:
        resultados = p.map(analizar_flujo_temporal, tasks)

    conteo = Counter(r['status'] for r in resultados)
    print("\n--- Resumen por motivo ---")
    for status, n in conteo.most_common():
        print(f"  {status}: {n}")

    print("\n--- Ejemplos por motivo de fallo (hasta 3 c/u) ---")
    for status in conteo:
        if status == 'ok':
            continue
        ejemplos = [r for r in resultados if r['status'] == status][:3]
        print(f"\n[{status}]")
        for e in ejemplos:
            print(f"  sim_hash={e['sim_hash']} tirada={e['tirada']} detail={e.get('detail', '')}")

    con_descarte = [r for r in resultados if r.get('lineas_descartadas')]
    if con_descarte:
        print(f"\n--- flx.dat con lineas invalidas descartadas ({len(con_descarte)} simulaciones) ---")
        for r in con_descarte:
            print(f"  sim_hash={r['sim_hash']} tirada={r['tirada']}: {r['lineas_descartadas']} lineas  ({r['path']})")

    print(f"\n✅ Proceso terminado. {conteo.get('ok', 0)} gráficos generados de {len(df_idx)} filas.")


if __name__ == "__main__":
    main()