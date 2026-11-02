import os
import re
from dataclasses import dataclass, field
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # estilo_figuras.py y rutas.py viven en analysis/
import rutas
from estilo_figuras import (SYSTEM_STYLE, ENGLISH_LABELS, SYSTEM_NAME_MAP,
                            set_plot_fonts, aplicar_estilo, ERRORBAR_KW,
                            REF_LINE_KW, SEQ_STYLE, anotar_chi)

# =============================================================================
# CONFIGURACIÓN GLOBAL
# =============================================================================

@dataclass
class AnalysisConfig:
    """Rutas y estilos comunes para todas las figuras."""
    data_path: str = os.path.join(rutas.AVALANCHE_DIR, "estadisticas_flujo_agregadas.csv")
    arcos_path: str = os.path.join(rutas.AVALANCHE_DIR, "resultados_completos_arcos.csv")
    raw_flow_path: str = os.path.join(rutas.AVALANCHE_DIR, "resultados_completos_esqueleto.csv")
    out_dir: str = rutas.GRAFICOS_FINALES_DIR

    # Estilos: vienen de estilo_figuras.py (comun a todos los scripts)
    system_style: Dict[str, Tuple[str, str]] = field(default_factory=lambda: dict(SYSTEM_STYLE))
    english_labels: Dict[str, str] = field(default_factory=lambda: dict(ENGLISH_LABELS))

CFG = AnalysisConfig()

# Esquema NUEVO de los CSV corregidos.
AR_COL = "area_ratio_nominal"
RR_COL = "r_ratio_nominal"

# Sufijo de las figuras regeneradas con la notacion nueva de arcos
# (chi_arch, L_arch, h_arch) o con etiquetas corregidas, para no pisar
# las versiones anteriores.
SUFIJO_NUEVO = "_nuevo"

# SYSTEM_NAME_MAP viene de estilo_figuras.py

# =============================================================================
# HELPERS
# =============================================================================

def save_plot(filename: str):
    os.makedirs(CFG.out_dir, exist_ok=True)
    plt.savefig(f"{CFG.out_dir}/{filename}.pdf", bbox_inches='tight')
    plt.close()


def require_columns(df: pd.DataFrame, cols: List[str], context: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"Faltan columnas en {context}: {missing}\n"
            f"Columnas disponibles: {df.columns.tolist()}"
        )


def prepare_flow_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Valida el CSV estadisticas_flujo_agregadas.csv generado por flow_analysis.py.

    Esquema esperado:
        sistema, chi, r_ratio_nominal, area_ratio_nominal, D_d_sq, grano_base,
        S_mean, S_sem, log_S, log_S_err, ...
    """
    require_columns(
        df,
        ['sistema', 'chi', AR_COL, RR_COL, 'D_d_sq', 'grano_base', 'S_mean', 'S_sem', 'log_S', 'log_S_err'],
        'estadisticas_flujo_agregadas.csv',
    )
    return df.copy()


def prepare_arcos(df_arcos: pd.DataFrame) -> pd.DataFrame:
    """
    Valida y normaliza resultados_completos_arcos.csv generado por arch_consolidator.py.

    Ese archivo debe traer columnas reales y nominales desde master_sim_index.csv.
    Para agrupar/graficar se usan SOLO las nominales.
    """
    require_columns(
        df_arcos,
        ['sim_hash', 'system_type', 'chi', AR_COL, RR_COL,
         'Total_Particulas_Arco', 'Particulas_Tipo_0', 'Particulas_Tipo_1',
         'Tortuosidad', 'Relacion_Aspecto', 'Altura', 'Span', 'Longitud'],
        'resultados_completos_arcos.csv',
    )

    out = df_arcos.copy()
    # Relacion de aspecto = (L/2)/H: mitad del ancho (Span, distancia horizontal
    # entre los centros de los granos de los extremos) sobre la altura (Altura,
    # entre los centros mas alto y mas bajo). Se recalcula aca para no depender
    # de la definicion con que se genero el csv (antes era H/L). Arcos con
    # H = 0 (todos los centros a la misma altura) quedan sin valor.
    altura = out['Altura'].astype(float)
    out['Relacion_Aspecto'] = np.where(altura > 0, 0.5 * out['Span'].astype(float) / altura.where(altura > 0, 1.0),
                                       np.nan)
    out['sistema'] = out['system_type'].map(SYSTEM_NAME_MAP)
    out = out.dropna(subset=['sistema'])
    return out


def get_pure_ref(df: pd.DataFrame, system_name: str, orifice: float, col: str = 'log_S'):
    """Referencia del sistema puro de la matriz correspondiente."""
    matrix_type = system_name.split('(')[0]
    ref = df[(df['chi'] == 0.0) &
             (df['sistema'].str.startswith(matrix_type)) &
             (df['D_d_sq'] == orifice)]
    if ref.empty:
        return None, None

    if col == 'log_S':
        err_col = 'log_S_err'
    elif col == 'S_mean':
        err_col = 'S_sem'
    else:
        err_col = col + '_err'

    if err_col not in ref.columns:
        return ref[col].mean(), None

    return ref[col].mean(), ref[err_col].mean()


def filter_standard_size(df: pd.DataFrame, target_ar: float, atol: float, system_col: str = 'sistema') -> pd.DataFrame:
    """
    Filtro estándar para comparar sistemas a tamaño aproximadamente fijo.

    Usa area_ratio_nominal. Esto es importante: no usar area_ratio real para seleccionar
    curvas porque separa artificialmente 0.0998 y 0.0999.
    """
    require_columns(df, [system_col, AR_COL], 'filter_standard_size')

    mask_disco = (
        (df[system_col] == 'Disco(G)-Disco(C)') &
        np.isclose(df[AR_COL].astype(float), 0.1, atol=0.03)
    )
    mask_otros = (
        (df[system_col] != 'Disco(G)-Disco(C)') &
        np.isclose(df[AR_COL].astype(float), target_ar, atol=atol)
    )
    return df[mask_disco | mask_otros].copy()


def select_by_nominal_ar(df: pd.DataFrame, target: float, atol: float = 0.05) -> pd.DataFrame:
    require_columns(df, [AR_COL], 'select_by_nominal_ar')
    return df[np.isclose(df[AR_COL].astype(float), target, atol=atol)].copy()

# =============================================================================
# GRÁFICOS DE FLUJO
# =============================================================================

def plot_puros_vs_orificio(df: pd.DataFrame):
    print("  -> Plotting 1: Pure Systems vs Orifice...")
    df_puros = df[(df['chi'] == 0.0) & (df['sistema'].str.startswith('Disco'))].sort_values('D_d_sq')
    if df_puros.empty:
        return

    plt.figure()
    plt.errorbar(
        df_puros['D_d_sq'], df_puros['log_S'], yerr=df_puros['log_S_err'],
        fmt='o-', color='black', label='Pure Disks', **ERRORBAR_KW
    )
    plt.xlabel(r"Orifice Ratio $(D/d)^2$")
    plt.ylabel(r"$\ln(\langle s \rangle + 1)$")
    plt.legend(title="System")
    save_plot("1-pure-vs-orifice")


def plot_disco_vs_area_ratio(df: pd.DataFrame):
    print("  -> Plotting 2: Disco(G)-Disco(C) vs Area Ratio...")

    df_sys = df[df['sistema'] == 'Disco(G)-Disco(C)'].copy()
    if df_sys.empty:
        return

    main_dd = df_sys['D_d_sq'].value_counts().idxmax()
    df_sys = df_sys[df_sys['D_d_sq'] == main_dd]

    plt.figure()

    val_ref, _ = get_pure_ref(df, 'Disco(G)-Disco(C)', main_dd, 'log_S')
    if val_ref is not None:
        plt.axhline(y=val_ref, **REF_LINE_KW, label=r'$\chi = 0$ (Pure)')

    target_chis = [0.01, 0.1, 0.14, 0.2]
    colors, markers = SEQ_STYLE(len(target_chis))

    for chi, color, marker in zip(target_chis, colors, markers):
        grp = df_sys[np.isclose(df_sys['chi'], chi, atol=0.005)].sort_values(AR_COL)
        if grp.empty:
            continue

        plot_data = (
            grp.groupby(AR_COL, dropna=False)
               .agg(log_S=('log_S', 'mean'), log_S_err=('log_S_err', 'mean'))
               .reset_index()
               .sort_values(AR_COL)
        )

        plt.errorbar(
            plot_data[AR_COL], plot_data['log_S'], yerr=plot_data['log_S_err'],
            fmt=f'{marker}-', color=color, label=rf'$\chi = {chi}$',
            **ERRORBAR_KW,
        )

    plt.xlim(left=0.0)
    plt.xlabel(r"Area Ratio $R_{\rm a}$")
    plt.ylabel(r"$\ln(\langle s \rangle + 1)$")

    plt.legend(loc='best')

    save_plot("2-discog-discoc-vs-area-ratio" + SUFIJO_NUEVO)


def _puntos_de(curvas, n_seg=30):
    """Puntos de las curvas (R_a, x, y, err): los datos, los extremos de las
    barras de error y puntos intermedios sobre las lineas que los unen."""
    puntos = []
    for _, x, y, err in curvas:
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        e = np.nan_to_num(np.asarray(err, dtype=float), nan=0.0)
        for xi, yi, ei in zip(x, y, e):
            puntos += [(xi, yi), (xi, yi - ei), (xi, yi + ei)]
        t = np.linspace(0.0, 1.0, n_seg)
        for i in range(len(x) - 1):
            puntos += list(zip(x[i] + t * (x[i + 1] - x[i]), y[i] + t * (y[i + 1] - y[i])))
    return np.asarray(puntos, dtype=float)


def _puntos_tapados(ax, puntos, titulo, margen_px=10.0):
    """Cuantos puntos (en coordenadas de datos) quedan debajo de la leyenda
    arriba a la derecha, con un margen en pixeles por marcadores y barras."""
    leg = ax.legend(loc='upper right', title=titulo)
    ax.figure.canvas.draw()
    caja = leg.get_window_extent()
    disp = ax.transData.transform(puntos)
    return int(np.sum((disp[:, 0] >= caja.x0 - margen_px) & (disp[:, 0] <= caja.x1 + margen_px) &
                      (disp[:, 1] >= caja.y0 - margen_px) & (disp[:, 1] <= caja.y1 + margen_px)))


def margen_para_leyenda(ax, curvas, titulo, paso=0.05, max_extra=1.5):
    """Fraccion del rango de y que hay que agregar arriba para que la leyenda
    entre arriba a la derecha sin tapar curvas (0 si ya entra). No cambia los ejes."""
    puntos = _puntos_de(curvas)
    y0, y1 = ax.get_ylim()
    extra = 0.0
    while extra <= max_extra:
        ax.set_ylim(y0, y1 + extra * (y1 - y0))
        if _puntos_tapados(ax, puntos, titulo) == 0:
            break
        extra += paso
    ax.set_ylim(y0, y1)
    return min(extra, max_extra)


def plot_sweeps_by_system(
    df: pd.DataFrame,
    target_systems: List[Tuple[str, str]],
    zoom: bool = False,
    xlim: Tuple[float, float] = (0.05, 0.15),
):
    if zoom:
        print("  -> Plotting 4 zoom: Clogging Curves (chi 0.05 - 0.15)...")
    else:
        print("  -> Plotting 4: Clogging Curves (Sweeps)...")

    # Chis que NO querés mostrar en el gráfico general.
    # La versión zoom sí los conserva.
    excluded_chis_no_zoom = [0.07, 0.075, 0.085, 0.09, 0.095, 0.105, 0.11, 0.15]

    # =========================================================================
    # PRIMERA PASADA: arma las curvas de cada sistema (sin dibujar todavia) y
    # junta, de paso, el rango global de x e y (contemplando barra de error)
    # entre TODOS los sistemas de esta llamada. Con eso, las 4 figuras que
    # salen de una misma llamada (la tanda sin zoom, o la tanda con zoom)
    # quedan con el mismo rango de ejes y son directamente comparables entre
    # si. Las dos llamadas (zoom=False y zoom=True) siguen siendo rangos
    # independientes entre si, a proposito: forzar el mismo rango ancho de la
    # tanda sin zoom sobre la tanda con zoom anularia el propio sentido del
    # zoom.
    # =========================================================================
    datos_por_sistema = {}
    x_all, y_lo_all, y_hi_all = [], [], []

    for system_name, _ in target_systems:
        df_mix = df[(df['sistema'] == system_name) & (df['chi'] > 0.0)].copy()

        if df_mix.empty:
            continue

        # En el gráfico SIN zoom, quitamos los chi intermedios.
        # En el gráfico CON zoom, los dejamos para ver el detalle.
        if not zoom:
            mask_excluded = np.zeros(len(df_mix), dtype=bool)
            for chi_excl in excluded_chis_no_zoom:
                mask_excluded |= np.isclose(df_mix['chi'].astype(float), chi_excl, atol=1e-6)

            df_mix = df_mix[~mask_excluded].copy()

        if df_mix.empty:
            continue

        main_dd = df_mix['D_d_sq'].value_counts().idxmax()
        p0, e0 = get_pure_ref(df, system_name, main_dd, 'log_S')

        if p0 is None:
            continue

        # R_a = 0.5 no se incluye: solo existe para Disco-Disco y no se puede
        # comparar entre sistemas.
        targets = [0.1, 0.2, 0.4]
        curvas = []

        for target in targets:
            grp_raw = select_by_nominal_ar(df_mix, target, atol=0.06)

            if grp_raw.empty:
                continue

            available = np.sort(grp_raw[AR_COL].dropna().unique())
            chosen = available[np.argmin(np.abs(available - target))]
            grp_raw = grp_raw[np.isclose(grp_raw[AR_COL], chosen, atol=1e-12)]

            grp = (
                grp_raw.groupby('chi')
                       .agg(log_S=('log_S', 'mean'), log_S_err=('log_S_err', 'mean'))
                       .reset_index()
                       .sort_values('chi')
            )

            if grp.empty:
                continue

            x = [0.0] + grp['chi'].tolist()
            y = [p0] + grp['log_S'].tolist()
            err = [e0] + grp['log_S_err'].tolist()
            curvas.append((chosen, x, y, err))

            print(f"\nSistema: {system_name}\n  R_a {target} -> {AR_COL} real {chosen:.4f}")

            x_all.extend(x)
            y_arr = np.asarray(y, dtype=float)
            err_arr = np.nan_to_num(np.asarray(err, dtype=float), nan=0.0)
            y_lo_all.extend((y_arr - err_arr).tolist())
            y_hi_all.extend((y_arr + err_arr).tolist())

        if curvas:
            datos_por_sistema[system_name] = curvas

    if not datos_por_sistema:
        return

    # Rango de x compartido: en la version zoom se respeta el xlim explicito
    # que ya pasa el caller (es la razon de ser del zoom); en la version sin
    # zoom se calcula de los datos, con un margen del 5%.
    if zoom:
        x_range = xlim
    else:
        x_min, x_max = min(x_all), max(x_all)
        pad_x = (x_max - x_min) * 0.05 if x_max > x_min else 0.02
        x_range = (x_min - pad_x, x_max + pad_x)

    # Rango de y compartido, con margen del 8% incluyendo las barras de error.
    y_min, y_max = min(y_lo_all), max(y_hi_all)
    pad_y = (y_max - y_min) * 0.08 if y_max > y_min else 0.1
    y_range = (y_min - pad_y, y_max + pad_y)

    # =========================================================================
    # SEGUNDA PASADA: una figura por sistema, todas con el mismo x_range/y_range.
    # =========================================================================
    figuras = []
    for system_name, _ in target_systems:
        if system_name not in datos_por_sistema:
            continue
        system_id = system_name.lower().replace('(', '').replace(')', '').replace('_', '-').replace(' ', '-')

        fig = plt.figure()
        ax = plt.gca()
        curvas = datos_por_sistema[system_name]
        colors, markers = SEQ_STYLE(len(curvas))
        for (chosen, x, y, err), color, marker in zip(curvas, colors, markers):
            plt.errorbar(x, y, yerr=err, fmt=f'{marker}-', color=color,
                         label=rf"$R_{{\rm a}} = {chosen:.3g}$", **ERRORBAR_KW)

        plt.xlabel(r"Concentration $\chi$")
        plt.ylabel(r"$\ln(\langle s \rangle + 1)$")
        plt.xlim(*x_range)
        plt.ylim(*y_range)
        titulo = CFG.english_labels.get(system_name, system_name)
        nombre = (f"4-{system_id}-clogging-curves-zoom-chi-0p05-0p15" if zoom
                  else f"4-{system_id}-clogging-curves")
        figuras.append((fig, ax, curvas, titulo, nombre))

    # Cartel con el sistema (ej. "Triangle-Disk") como titulo de la leyenda,
    # adentro, arriba a la derecha. Si en alguna figura tapa curvas, se agrega
    # el mismo margen arriba en TODAS (siguen con el mismo eje y).
    extra = max(margen_para_leyenda(ax, curvas, titulo) for _, ax, curvas, titulo, _ in figuras)
    for fig, ax, curvas, titulo, nombre in figuras:
        ax.set_ylim(y_range[0], y_range[1] + extra * (y_range[1] - y_range[0]))
        ax.legend(loc='upper right', title=titulo)
        plt.figure(fig.number)
        save_plot(nombre)


def plot_geometric_effects(df: pd.DataFrame, target_ar: float = 0.15, total: bool = False, masa: bool = False):
    """5a/5b: tamano medio de avalancha normalizado por el del sistema puro.
    total=False: <s> = granos grandes por avalancha (S_mean).
    total=True:  <s_T> = todos los granos, grandes + chicos (ST_mean); solo 5b,
                 archivo 5b-matrix-normalized-total.
    masa=True:   lo mismo en masa (MB_mean / MT_mean de flow_analysis.py); solo 5b,
                 archivos 5b-matrix-normalized-mass y 5b-matrix-normalized-total-mass.
                 Para los grandes, M/M0 = s/s0 (la masa del grano grande se cancela)."""
    if masa:
        mean_col, sem_col = ('MT_mean', 'MT_sem') if total else ('MB_mean', 'MB_sem')
    else:
        mean_col, sem_col = ('ST_mean', 'ST_sem') if total else ('S_mean', 'S_sem')
    if mean_col not in df.columns:
        print(f"  [!] {mean_col} no esta en estadisticas_flujo_agregadas.csv: correr flow_analysis.py de nuevo.")
        return
    print(f"  -> Plotting 5{'b' if total or masa else 'a/b'}{' total' if total else ''}{' mass' if masa else ''}: "
          f"Normalized Effects (R_a ~ {target_ar})...")

    df_mix = df[df['chi'] > 0].copy()
    if df_mix.empty:
        return

    main_dd = df_mix['D_d_sq'].value_counts().idxmax()
    df_puros = (
        df[df['chi'] == 0][['grano_base', 'D_d_sq', mean_col, sem_col]]
        .rename(columns={mean_col: 'S0', sem_col: 'S0_err'})
    )

    df_norm = pd.merge(df_mix[df_mix['D_d_sq'] == main_dd], df_puros, on=['grano_base', 'D_d_sq'])
    if df_norm.empty:
        return

    df_norm['log_S_norm'] = np.log((df_norm[mean_col] / df_norm['S0']) + 1)
    df_norm['log_S_norm_err'] = df_norm[sem_col] / (df_norm[mean_col] + df_norm['S0'])
    v = 'M' if masa else 's'
    ylabel = (rf"$\ln(\langle {v}_{{\rm T}} \rangle / \langle {v}_{{\rm T,0}} \rangle + 1)$" if total
              else rf"$\ln(\langle {v} \rangle / \langle {v}_0 \rangle + 1)$")

    def draw(sistemas: List[str], filename: str):
        plt.figure()
        for sis in sistemas:
            grp = filter_standard_size(df_norm[df_norm['sistema'] == sis], target_ar, 0.07)
            if grp.empty:
                continue
            grp_sis = (
                grp.groupby('chi')
                   .agg(log_S_norm=('log_S_norm', 'mean'), log_S_norm_err=('log_S_norm_err', 'mean'))
                   .reset_index()
                   .sort_values('chi')
            )
            if grp_sis.empty:
                continue

            color, marker = CFG.system_style.get(sis, ('black', 'o'))
            label_en = CFG.english_labels.get(sis, sis)
            plt.errorbar(
                [0.0] + grp_sis['chi'].tolist(),
                [np.log(2.0)] + grp_sis['log_S_norm'].tolist(),
                yerr=[0.0] + grp_sis['log_S_norm_err'].tolist(),
                fmt=f'{marker}-', color=color, label=label_en, **ERRORBAR_KW,
            )

        plt.xlabel(r"Concentration $\chi$")
        plt.ylabel(ylabel)
        if not total and not masa:
            plt.ylim(0.4, 2.5)
        plt.legend(loc='upper left')
        save_plot(filename)

    matriz = ['Disco(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Tri(G)-Disco(C)', 'Oct(G)-Disco(C)']
    if total or masa:
        draw(matriz, "5b-matrix-normalized" + ("-total" if total else "") + ("-mass" if masa else ""))
        return
    draw(['Disco(G)-Disco(C)', 'Disco(G)-Cuad(C)', 'Disco(G)-Pent(C)', 'Disco(G)-Tri(C)', 'Disco(G)-Oct(C)'], "5a-additive-normalized")
    draw(matriz, "5b-matrix-normalized")


def plot_geometric_effects_raw(df: pd.DataFrame, target_ar: float = 0.15):
    print(f"  -> Plotting 5c: Raw Log-Size Matrix Effect (R_a ~ {target_ar})...")
    main_dd = df['D_d_sq'].value_counts().idxmax()

    plt.figure()
    for sis in ['Disco(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Tri(G)-Disco(C)', 'Oct(G)-Disco(C)']:
        grp = filter_standard_size(df[df['D_d_sq'] == main_dd], target_ar, 0.07)
        grp_sis = (
            grp[grp['sistema'] == sis]
            .groupby('chi')
            .agg(log_S=('log_S', 'mean'), log_S_err=('log_S_err', 'mean'))
            .reset_index()
            .sort_values('chi')
        )
        if grp_sis.empty:
            continue

        p0, e0 = get_pure_ref(df, sis, main_dd, 'log_S')
        if p0 is None:
            continue

        color, marker = CFG.system_style.get(sis, ('black', 'o'))
        label_en = CFG.english_labels.get(sis, sis)
        plt.errorbar(
            [0.0] + grp_sis['chi'].tolist(),
            [p0] + grp_sis['log_S'].tolist(),
            yerr=[e0] + grp_sis['log_S_err'].tolist(),
            fmt=f'{marker}-', color=color, label=label_en, **ERRORBAR_KW,
        )

    plt.xlabel(r"Concentration $\chi$")
    plt.ylabel(r"$\ln(\langle s \rangle + 1)$")
    plt.ylim(2.2, 6.5)
    plt.legend()
    save_plot("5c-matrix-raw-log")


def plot_geometric_effects_linear(df: pd.DataFrame, target_ar: float = 0.15):
    print(f"  -> Plotting 11: Linear Scale (R_a ~ {target_ar})...")
    main_dd = df['D_d_sq'].value_counts().idxmax()

    def draw_lin(sistemas: List[str], filename: str):
        plt.figure()
        for sis in sistemas:
            grp = filter_standard_size(df[df['D_d_sq'] == main_dd], target_ar, 0.07)
            grp_sis = (
                grp[grp['sistema'] == sis]
                .groupby('chi')
                .agg(S_mean=('S_mean', 'mean'), S_sem=('S_sem', 'mean'))
                .reset_index()
                .sort_values('chi')
            )
            if grp_sis.empty:
                continue

            p0, e0 = get_pure_ref(df, sis, main_dd, 'S_mean')
            if p0 is None:
                continue

            color, marker = CFG.system_style.get(sis, ('black', 'o'))
            label_en = CFG.english_labels.get(sis, sis)
            plt.errorbar(
                [0.0] + grp_sis['chi'].tolist(),
                [p0] + grp_sis['S_mean'].tolist(),
                yerr=[e0] + grp_sis['S_sem'].tolist(),
                fmt=f'{marker}-', color=color, label=label_en, **ERRORBAR_KW,
            )

        plt.xlabel(r"Concentration $\chi$")
        plt.ylabel(r"Mean Size $\langle s \rangle$")
        plt.legend(loc='upper left')
        save_plot(filename)

    draw_lin(['Disco(G)-Disco(C)', 'Disco(G)-Cuad(C)', 'Disco(G)-Pent(C)', 'Disco(G)-Tri(C)', 'Disco(G)-Oct(C)'], "11a-additive-linear")
    draw_lin(['Disco(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Tri(G)-Disco(C)', 'Oct(G)-Disco(C)'], "11b-matrix-linear")

# =============================================================================
# GRÁFICOS DE ARCOS
# =============================================================================

def plot_arch_analysis(df_stats: pd.DataFrame):
    """Plot 6, 9 y 10: composición, tortuosidad y tamaños de avalancha."""
    if not os.path.exists(CFG.arcos_path):
        return

    print("\n[+] Starting Arch Analysis (Target chi = 0.2)...")
    df_arcos = prepare_arcos(pd.read_csv(CFG.arcos_path))

    matrix_systems = ['Disco(G)-Disco(C)', 'Tri(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Oct(G)-Disco(C)']
    df_filtered = df_arcos[df_arcos['sistema'].isin(matrix_systems)].copy()

    # Filtro general: disco-disco usa R_a nominal ~0.1; el resto target_ar ~0.15.
    df_mix = filter_standard_size(df_filtered[df_filtered['chi'] > 0], 0.15, 0.1, system_col='sistema')
    if df_mix.empty:
        print("  [!] No hay datos de arcos para el filtro estándar de tamaño.")
        return

    # --- Plot 6: Composición ---
    # Filtrar el DataFrame para conservar solo los valores de chi entre 0 y 0.4
    df_filtered = df_mix[(df_mix['chi'] >= 0) & (df_mix['chi'] <= 0.4)].copy()
    
    df_filtered['frac_arco'] = df_filtered['Particulas_Tipo_1'] / df_filtered['Total_Particulas_Arco']
    plt.figure()
    for sys_name in matrix_systems:
        group = df_filtered[df_filtered['sistema'] == sys_name]
        if group.empty:
            continue
        grp = group.groupby('chi')['frac_arco'].agg(['mean', 'std', 'count']).reset_index()
        grp['sem'] = grp['std'] / np.sqrt(grp['count'])
        color, marker = CFG.system_style.get(sys_name, ('black', 'o'))
        plt.errorbar(
            grp['chi'], grp['mean'], yerr=grp['sem'],
            fmt=f'{marker}-', color=color,
            label=CFG.english_labels.get(sys_name, sys_name), **ERRORBAR_KW,
        )

    # Ajustar la línea de referencia homogénea para que vaya de 0 a 0.4
    plt.plot([0, 0.5], [0, 0.5], **REF_LINE_KW, label='Homogeneous')
    
    # Opcional: fijar explícitamente los límites del eje x
    plt.xlim(0, 0.45)
    
    plt.xlabel(r"Global Concentration $\chi$")
    plt.ylabel(r"Arch Additive Fraction $\chi_{\rm arch}$")
    plt.legend(loc='upper left')
    save_plot("6-arch-composition" + SUFIJO_NUEVO)

    # --- Filtro chi = 0.2 ---
    target_chi = 0.2
    df_target = df_mix[np.isclose(df_mix['chi'], target_chi, atol=0.05)].copy()
    if df_target.empty:
        print(f"  [!] Advertencia: No hay datos de arcos para chi ~ {target_chi}")
        return

    color_map = {sys: CFG.system_style[sys][0] for sys in matrix_systems}
    present = [s for s in matrix_systems if s in df_target['sistema'].unique()]
    new_ticks = [CFG.english_labels.get(s, s) for s in present]

    # --- Plot 9: Tortuosidad ---
    plt.figure()
    ax = sns.violinplot(
        data=df_target, x='sistema', y='Tortuosidad', hue='sistema',
        palette=color_map, inner='quartile', order=present, legend=False,
    )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(new_ticks, rotation=15)
    anotar_chi(ax, target_chi)
    plt.xlabel(None)
    plt.ylabel("Tortuosity")
    save_plot(f"9-arch-tortuosity-chi-{target_chi}")

    # --- Plot 10: Avalanche Boxplot ---
    if os.path.exists(CFG.raw_flow_path):
        df_raw = pd.read_csv(CFG.raw_flow_path)
        require_columns(df_raw, ['sim_hash', 'tamano_big'], 'resultados_completos_esqueleto.csv')
        df_plot = pd.merge(
            df_raw,
            df_target[['sim_hash', 'sistema', 'chi']].drop_duplicates(),
            on='sim_hash',
            how='inner',
        )
        if not df_plot.empty:
            plt.figure()
            ax_box = sns.boxplot(
                data=df_plot, x='sistema', y='tamano_big', hue='sistema',
                palette=color_map, order=present, legend=False,
            )
            ax_box.set_xticks(range(len(present)))
            ax_box.set_xticklabels(new_ticks, rotation=15)
            anotar_chi(ax_box, target_chi)
            plt.xlabel(None)
            plt.yscale('log')
            plt.ylabel(r"Avalanche Size $s$")
            save_plot(f"10-avalanche-boxplot-chi-{target_chi}")


def plot_arch_geometry_deep_dive(df_stats: pd.DataFrame):
    """Gráficos 12, 13 y 14 de geometría de arcos para chi = 0.2."""
    if not os.path.exists(CFG.arcos_path):
        return

    print("[+] Generando Geometría de Arcos (chi = 0.2)...")
    df_arcos = prepare_arcos(pd.read_csv(CFG.arcos_path))

    matrix_systems = ['Disco(G)-Disco(C)', 'Tri(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Oct(G)-Disco(C)']
    target_chi = 0.2
    df_plot = df_arcos[
        (df_arcos['sistema'].isin(matrix_systems)) &
        (np.isclose(df_arcos['chi'], target_chi, atol=0.05))
    ].copy()
    if df_plot.empty:
        return

    color_map = {sys: CFG.system_style[sys][0] for sys in matrix_systems}
    present = [s for s in matrix_systems if s in df_plot['sistema'].unique()]
    new_ticks = [CFG.english_labels.get(s, s) for s in present]

    # --- Plot 12: Longitud vs Span ---
    plt.figure()
    for sys_name in present:
        data = df_plot[df_plot['sistema'] == sys_name]
        sns.regplot(
            data=data, x='Span', y='Longitud', label=CFG.english_labels[sys_name],
            color=color_map[sys_name], scatter_kws={'alpha': 0.05, 's': 2},
            line_kws={'lw': 2.5},
        )
    plt.xlabel("Orifice Span ($L$)")
    plt.ylabel("Arch Length ($s$)")
    plt.legend(loc='upper left')
    save_plot(f"12-arch-length-vs-span-chi-{target_chi}")

    # --- Plot 13: Relación de aspecto ---
    plt.figure()
    ax = sns.violinplot(
        data=df_plot, x='sistema', y='Relacion_Aspecto', hue='sistema',
        palette=color_map, inner='quartile', order=present, legend=False,
    )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(new_ticks, rotation=15)
    plt.xlabel(None)
    plt.ylabel("Arch Aspect Ratio $L/2H$")
    anotar_chi(ax, target_chi)
    save_plot(f"13-arch-aspect-ratio-chi-{target_chi}")

    # --- Plot 14: Altura ---
    plt.figure()
    ax = sns.boxplot(
        data=df_plot, x='sistema', y='Altura', hue='sistema',
        palette=color_map, order=present, legend=False,
    )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(new_ticks, rotation=15)
    plt.xlabel(None)
    plt.ylabel("Arch Height ($H$)")
    anotar_chi(ax, target_chi)
    save_plot(f"14-arch-height-boxplot-chi-{target_chi}")


def aspecto_medio(g, n_boot=500, seed=0):
    """Relacion de aspecto del arco medio, <L>/(2<H>), y su error estandar por
    bootstrap sobre los arcos (remuestreo con reposicion). Usa Span (L) y
    Altura (H) de cada arco."""
    L = g['Span'].to_numpy(float)
    H = g['Altura'].to_numpy(float)
    ok = np.isfinite(L) & np.isfinite(H)
    L, H = L[ok], H[ok]
    n = len(L)
    if n == 0 or H.mean() <= 0:
        return {'mean': np.nan, 'sem': np.nan, 'count': n}
    valor = L.mean() / (2.0 * H.mean())
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = L[idx].mean(axis=1) / (2.0 * H[idx].mean(axis=1))
    return {'mean': valor, 'sem': float(np.std(boot, ddof=1)), 'count': n}


W_REGEX = re.compile(r'_W_([0-9]+(?:p[0-9]+)?)')


def _w_de_path(path) -> float:
    """Ancho/radio de salida del nombre de carpeta (..._W_2p0 -> 2.0); NaN si no esta."""
    m = W_REGEX.search(str(path))
    return float(m.group(1).replace('p', '.')) if m else np.nan


def puros_como_mezcla(df_arcos: pd.DataFrame, df_mezclas: pd.DataFrame, sistemas: List[str]) -> pd.DataFrame:
    """Arcos de los sistemas puros (chi = 0) renombrados como la mezcla con el
    mismo grano grande ('Tri(G)' -> 'Tri(G)-Disco(C)'), para graficarlos como
    el punto chi = 0 de esa curva. Si el csv trae 'path', solo se usan los
    puros con el mismo orificio (_W_ de la carpeta) que las mezclas de ese
    sistema, por si hay puros de barridos de orificio."""
    partes = []
    for sis in sistemas:
        base = sis.split('-')[0]
        p = df_arcos[(df_arcos['sistema'] == base) & np.isclose(df_arcos['chi'].astype(float), 0.0)]
        if p.empty:
            print(f"  [!] Sin arcos del sistema puro {base} para el punto chi = 0 de {sis}.")
            continue
        if 'path' in df_arcos.columns:
            w_mezcla = df_mezclas.loc[df_mezclas['sistema'] == sis, 'path'].map(_w_de_path).dropna().unique()
            if len(w_mezcla):
                p = p[p['path'].map(_w_de_path).isin(w_mezcla)]
        print(f"  -> {sis}: {len(p)} arcos del sistema puro {base} como punto chi = 0")
        if not p.empty:
            partes.append(p.assign(sistema=sis))
    return pd.concat(partes, ignore_index=True) if partes else df_arcos.iloc[0:0]


def plot_arch_evolution_vs_chi(df_stats: pd.DataFrame):
    """Gráfico 15: evolución de relación de aspecto y tortuosidad vs chi."""
    if not os.path.exists(CFG.arcos_path):
        return

    print("[+] Generando Evolución (vs Chi) - Solo métricas clave...")
    df_arcos = prepare_arcos(pd.read_csv(CFG.arcos_path))

    matrix_systems = ['Disco(G)-Disco(C)', 'Tri(G)-Disco(C)', 'Cuad(G)-Disco(C)', 'Pent(G)-Disco(C)', 'Oct(G)-Disco(C)']
    df_filtered = filter_standard_size(
        df_arcos[df_arcos['sistema'].isin(matrix_systems)],
        0.15,
        0.1,
        system_col='sistema',
    )
    
    # Sistemas puros (chi = 0) como primer punto de cada mezcla. En el csv de
    # arcos los puros se llaman 'Disco(G)', 'Tri(G)', ... y tienen
    # area_ratio_nominal = 1.0, asi que filter_standard_size los descarta.
    df_filtered = pd.concat([df_filtered, puros_como_mezcla(df_arcos, df_filtered, matrix_systems)],
                            ignore_index=True)

    # Filtrar para conservar solo los valores de chi de 0 a 0.4
    df_filtered = df_filtered[(df_filtered['chi'] >= 0) & (df_filtered['chi'] <= 0.4)]

    metrics = {
        'Relacion_Aspecto': r'Arch Aspect Ratio $L_{\rm arch}/2h_{\rm arch}$',
        'Tortuosidad': 'Arch Tortuosity',
    }

    for metric, ylabel in metrics.items():
        plt.figure()
        for sys_name in matrix_systems:
            group = df_filtered[df_filtered['sistema'] == sys_name]
            if group.empty:
                continue

            if metric == 'Relacion_Aspecto':
                # Cociente de promedios <L>/(2<H>), no promedio de L/2H: un arco
                # casi plano (H ~ 0) da L/2H enorme y dominaba media y error.
                summary = (group.groupby('chi')
                           .apply(lambda g: pd.Series(aspecto_medio(g)), include_groups=False)
                           .reset_index())
            else:
                summary = group.groupby('chi')[metric].agg(['mean', 'std', 'count']).reset_index()
                summary['sem'] = summary['std'] / np.sqrt(summary['count'])
            summary = summary[summary['count'] >= 20]
            if summary.empty:
                continue

            color, marker = CFG.system_style.get(sys_name, ('black', 'o'))
            plt.errorbar(
                summary['chi'], summary['mean'], yerr=summary['sem'],
                fmt=f'{marker}-', color=color,
                label=CFG.english_labels.get(sys_name, sys_name), **ERRORBAR_KW,
            )

        # Fijar explícitamente los límites del eje x de 0 a 0.4
        plt.xlim(-0.01, 0.45)  # margen para que el punto chi = 0 no quede cortado
        
        plt.xlabel(r"Concentration $\chi$")
        plt.ylabel(ylabel)
        plt.legend(loc='best')
        # Las dos figuras 15 ahora incluyen chi = 0: sufijo para no pisar las viejas
        save_plot(f"15-evolution-{metric.lower()}{SUFIJO_NUEVO}")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("[+] Starting integrated PDF pipeline...")

    aplicar_estilo(verbose=True)

    os.makedirs(CFG.out_dir, exist_ok=True)

    df = pd.read_csv(CFG.data_path)
    df = prepare_flow_stats(df)

    # Limpieza: sistema fuera del estudio actual.
    df = df[df['sistema'] != 'Cuad(G)-Cuad(C)']

    # Orden de los sistemas: Triangle, Square, Pentagon, Octagon y Disk-Disk al final.
    # (los que no tienen datos se saltean)
    sistemas_meta = [
        ('Tri(G)-Disco(C)', 'Triangle'),
        ('Cuad(G)-Disco(C)', 'Square'),
        ('Pent(G)-Disco(C)', 'Pentagon'),
        ('Oct(G)-Disco(C)', 'Octagon'),
        ('Disco(G)-Disco(C)', 'Disk'),
    ]

    # Flujo / avalanchas
    # plot_puros_vs_orificio(df)
    plot_disco_vs_area_ratio(df)
    plot_sweeps_by_system(df, sistemas_meta)
    plot_sweeps_by_system(df, sistemas_meta, zoom=True, xlim=(0.05, 0.15))
    plot_geometric_effects(df)
    plot_geometric_effects(df, total=True)
    plot_geometric_effects(df, masa=True)
    plot_geometric_effects(df, total=True, masa=True)
    plot_geometric_effects_raw(df)
    plot_geometric_effects_linear(df)

    # Arcos
    plot_arch_analysis(df)
    plot_arch_geometry_deep_dive(df)
    plot_arch_evolution_vs_chi(df)

    print("\n[+] Pipeline finished.")
    print(f"[+] All plots generated successfully in: {CFG.out_dir}")


if __name__ == "__main__":
    main()