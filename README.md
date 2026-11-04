# Atascos en silos 2D con mezclas binarias de formas

Simulación y análisis del flujo y atascamiento de medios granulares bidispersos (discos y polígonos) en un silo 2D vibrado. Se estudian avalanchas, arcos, fuerzas de contacto y el tensor de tensiones / presión en función de la fracción de área de partículas chicas (χ), la razón de tamaños y la forma.

- **Simulador:** C++ con **Box2D 2.4.2** (`simulation/`, ejecutable `eo-silo`).
- **Análisis:** Python (numpy, pandas, networkx, matplotlib, pyarrow) y algunos scripts de gnuplot.

## Estructura

```text
sim-binary-shapes-2d/
├── simulation/                     # Código C++ del simulador (Box2D)
│   ├── Makefile                    # make -> bin/eo-silo (bin/ y obj/ no se versionan)
│   └── src/
├── run_scripts/                    # Lanzar simulaciones
│   ├── params_template.in          # Plantilla de parámetros
│   ├── run_slurm_cluster_utn.sh    # Barridos en el cluster (SLURM)
│   ├── sweep_*.sh                  # Barridos locales (altura, ancho, orificio, mezcla)
│   ├── pack_simulations.sh, renderizar_all.sh
│   └── viejo.tar.gz                # Versión anterior del paquete UTN
├── analysis/
│   ├── 1_avalanchas_y_arcos/       # Pipeline principal (main_pipeline.py)
│   ├── 2_fuerzas/                  # Fuerza normal en arco / ROI por tipo de contacto
│   ├── 3_presion_y_mapas/          # Tensor de tensiones, presión y mapas de calor
│   ├── 4_correlaciones_gnuplot/    # Gráficos cruzados (avalanchas vs tortuosidad, Fn, ...)
│   └── utils/                      # Diagnósticos puntuales
└── processed_data/                 # Resultados (figuras); los CSV/npz/parquet quedan solo en disco
    ├── 5_graficos_finales/
    ├── fuerzas/   (pdf_v4/, pdf_en_paper/, datos_v4/*)
    ├── presion/
    └── mapas_calor/ (heatmap_*.png, frames_v4/*)
```
`*`: solo en disco local, no se versiona.

Cada carpeta de `analysis/` y `processed_data/` puede tener un `viejo.tar.gz` con las versiones anteriores de los scripts o figuras (para ver qué hay adentro: `tar tzf viejo.tar.gz`; para extraerlo: `tar xzf viejo.tar.gz`).

## Flujo de trabajo

1. **Compilar y simular**
   ```bash
   cd simulation && make            # requiere Box2D en ../../box2d y GSL
   # copiar bin/eo-silo junto a run_scripts/ y lanzar un barrido
   ```
2. **Avalanchas y arcos** (indexa simulaciones, detecta avalanchas y arcos, consolida y grafica)
   ```bash
   python analysis/1_avalanchas_y_arcos/main_pipeline.py
   python analysis/1_avalanchas_y_arcos/graficos_finales.py
   ```
   Genera `master_sim_index.csv` y `processed_data/1_avalanche_data/`.
3. **Fuerzas**
   ```bash
   python analysis/2_fuerzas/force_arch_v4.py --idx master_sim_index.csv --out processed_data/fuerzas/datos_v4
   bash   analysis/2_fuerzas/graficar_todas_configuraciones.sh   # -> processed_data/fuerzas/pdf_v4
   ```
4. **Presión y mapas de calor**
   ```bash
   python analysis/3_presion_y_mapas/presion_tensor_arch.py ...      # -> presion_por_particula.csv
   python analysis/3_presion_y_mapas/dividir_particulas_por_forma.py --csv presion_por_particula.csv --out_dir por_forma/
   python analysis/3_presion_y_mapas/hacer_mapas.py                  # todos los mapas de una vez
   python analysis/3_presion_y_mapas/comparar_presion_por_chi.py --csv resumen_presion_por_atasco.csv
   ```

Los scripts dentro de una misma carpeta de `analysis/` se importan entre sí: hay que mantenerlos juntos.

## Rutas

Todas las rutas salen de `analysis/rutas.py`, que las calcula a partir de dónde está clonado el repo: los scripts funcionan igual en `/home/juli/work/...`, en `/mnt/e/work/...` (WSL) o en el cluster.

- **Datos crudos fuera del repo:** `export SIM_DATA_DIR=/ruta/a/simulation_data/R_2.0` antes de correr el paso 1.
- **Otro contenedor de Apptainer** (scripts `.sh`): `export CONTENEDOR=/ruta/al/contenedor.sif`.
- **gnuplot** (`analysis/4_correlaciones_gnuplot/`): correrlo desde cualquier carpeta dentro del repo; la raíz la toma de `git rev-parse --show-toplevel`.

## Configuración en una PC nueva (WSL + datos en Windows / Google Drive)

El repo se clona en WSL; los datos pesados viven en Windows (y Google Drive los sube a la nube) y se ven desde el repo con enlaces simbólicos:

```bash
git clone https://github.com/montero-j/sim-binary-shapes-2d.git ~/work/lab/project/sim-binary-shapes-2d
cd ~/work/lab/project/sim-binary-shapes-2d
export PROYECTO_WIN="/mnt/c/Users/<usuario>/OneDrive/work/lab/project/sim-binary-shapes-2d"   # carpeta del proyecto en Windows
bash setup/enlazar_datos.sh "$PROYECTO_WIN"          # enlaza datos_simulaciones, csv, parquet...
mkdir -p ~/bin && cp setup/git-pushsync ~/bin/       # 'git pushsync' = push + copia del código a Windows
```

`PROYECTO_WIN` y `export PATH="$HOME/bin:$PATH"` van en `~/.bashrc`. Si la carpeta de Windows cambia de lugar, `setup/reubicar_enlaces.sh <viejo> <nuevo>` actualiza los enlaces. Después de cada commit, `git pushsync` en lugar de `git push`.

## Notas

- `setup/inventario_disco.sh`: inventario de solo lectura de todo lo de este trabajo que hay en el disco (home de WSL, OneDrive, Descargas, otros discos): clones de git, carpetas de datos, archivos grandes, duplicados, enlaces rotos y caches. Escribe `~/inventario_disco_<fecha>.txt`; no mueve ni borra nada.

- `migrar_datos_locales.sh`: correrlo una vez después del pull de la reorganización para recuperar los archivos pesados y acomodar los datos locales.
- `AUDITORIA_REORGANIZACION.md`: detalle de qué se movió, a dónde y por qué.
