import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib

matplotlib.use('Agg')

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # rutas.py vive en analysis/
import rutas
MEGA_DATASET_PATH = os.path.join(rutas.AVALANCHE_DIR, "resultados_completos_arcos.csv")
OUTPUT_DIR = rutas.ARCH_HIST_DIR


def ensure_nominal_columns(df):
    if 'area_ratio_nominal' not in df.columns and 'area_ratio' in df.columns:
        df['area_ratio_nominal'] = df['area_ratio'].round(3)
    if 'r_ratio_nominal' not in df.columns and 'r_ratio' in df.columns:
        df['r_ratio_nominal'] = df['r_ratio'].round(2)
    return df


def plot_arch_histograms(dataset_path=MEGA_DATASET_PATH, output_dir=OUTPUT_DIR):
    print("\n[+] Generando análisis avanzado de composición de arcos...")

    if not os.path.exists(dataset_path):
        print(f"  -> [!] Error: No se encontró el dataset en {dataset_path}")
        return

    df = pd.read_csv(dataset_path)
    df = ensure_nominal_columns(df)

    required_cols = [
        'system_type',
        'chi',
        'area_ratio_nominal',
        'Particulas_Tipo_1',
        'Total_Particulas_Arco',
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"  -> [!] Error: faltan columnas: {missing}")
        return

    df = df[df['Total_Particulas_Arco'] > 0].copy()
    df['Fraccion_Tipo_1'] = df['Particulas_Tipo_1'] / df['Total_Particulas_Arco']

    os.makedirs(output_dir, exist_ok=True)

    grouped = df.groupby(['system_type', 'chi', 'area_ratio_nominal'], dropna=False)
    graficos_generados = 0

    for (sistema, chi, area_ratio_nominal), df_plot in grouped:
        sis_safe = str(sistema).replace(" ", "_").replace("(", "").replace(")", "")
        chi_safe = str(chi).replace(".", "p")
        ar_safe = str(area_ratio_nominal).replace(".", "p")

        subfolder = os.path.join(output_dir, sis_safe)
        os.makedirs(subfolder, exist_ok=True)

        plt.figure(figsize=(10, 6))
        bins = np.arange(0, df_plot['Total_Particulas_Arco'].max() + 2) - 0.5

        plt.hist(df_plot['Particulas_Tipo_0'], bins=bins, alpha=0.5, label='Grandes (Matriz)', color='#1f77b4')
        plt.hist(df_plot['Particulas_Tipo_1'], bins=bins, alpha=0.7, label='Chicas (Aditivo)', color='#ff7f0e')

        plt.title(f"Composición de Arcos: {sistema}\n$\\chi$={chi} | AR={area_ratio_nominal}", fontsize=13)
        plt.xlabel("Cantidad de partículas en el arco")
        plt.ylabel("Frecuencia")
        plt.xticks(np.arange(0, df_plot['Total_Particulas_Arco'].max() + 1))
        plt.legend()
        plt.grid(axis='y', linestyle='--', alpha=0.4)
        plt.savefig(os.path.join(subfolder, f"hist_chi{chi_safe}_AR{ar_safe}.png"), dpi=200)
        plt.close()

        plt.figure(figsize=(9, 5))
        sns.boxplot(data=df_plot, x='Total_Particulas_Arco', y='Fraccion_Tipo_1', color='lightgrey', showfliers=False)
        sns.stripplot(data=df_plot, x='Total_Particulas_Arco', y='Fraccion_Tipo_1', alpha=0.3, jitter=True)
        plt.axhline(chi, color='red', linestyle='--', label=fr'$\chi$={chi}')
        plt.title(f"Fracción de Aditivos por Tamaño de Arco\n{sistema} | AR={area_ratio_nominal}")
        plt.ylim(-0.05, 1.05)
        plt.legend()
        plt.savefig(os.path.join(subfolder, f"fraccion_chi{chi_safe}_AR{ar_safe}.png"), dpi=200)
        plt.close()

        graficos_generados += 2

    print("  -> Generando comparativas de barrido de área...")
    for (sistema, chi), df_trend in df.groupby(['system_type', 'chi'], dropna=False):
        stats = df_trend.groupby('area_ratio_nominal')['Fraccion_Tipo_1'].agg(['mean', 'std']).reset_index()
        stats = stats.sort_values('area_ratio_nominal')

        plt.figure(figsize=(10, 6))
        plt.errorbar(stats['area_ratio_nominal'], stats['mean'], yerr=stats['std'], fmt='o-', capsize=5, label='Composición media')
        plt.axhline(chi, color='red', linestyle=':', label=fr'Mezcla ($\chi$={chi})')
        plt.title(f"Composición de Arcos vs Area Ratio\n{sistema} | $\chi$ fijo = {chi}")
        plt.xlabel("Area Ratio")
        plt.ylabel("Fracción de aditivos en arcos")
        plt.grid(True, alpha=0.3)
        plt.legend()

        sys_safe = str(sistema).replace(" ", "_").replace("-", "_")
        chi_safe = str(chi).replace('.', 'p')
        plt.savefig(os.path.join(output_dir, f"tendencia_area_{sys_safe}_chi{chi_safe}.png"), dpi=200)
        plt.close()

    print(f"  -> Éxito: se generaron {graficos_generados} gráficos individuales.")


if __name__ == "__main__":
    plot_arch_histograms()
