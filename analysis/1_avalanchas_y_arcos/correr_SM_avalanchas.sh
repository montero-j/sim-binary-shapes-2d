#!/bin/bash
# Analisis de avalanchas del Supplemental:
#   SII     correlacion entre avalanchas consecutivas (correlacion_avalanchas_formas.py)
#   SIII.C  distribucion de tamaños, CCDF y ajuste exponencial (ccdf_avalanchas.py)
#
# Lee processed_data/1_avalanche_data/resultados_completos_esqueleto.csv y
# master_sim_index.csv (los genera main_pipeline.py: bash analysis/correr_todo.sh avalanchas).
# No pisa nada: si las salidas ya existen, el script de python corta y avisa.
#
# Uso (desde cualquier carpeta):
#   bash analysis/1_avalanchas_y_arcos/correr_SM_avalanchas.sh            # los dos
#   bash analysis/1_avalanchas_y_arcos/correr_SM_avalanchas.sh correlacion
#   bash analysis/1_avalanchas_y_arcos/correr_SM_avalanchas.sh ccdf
# Variables: CONTENEDOR (default ~/work/lab/containers/workstation_full.sif;
#            CONTENEDOR=no para usar el python3 del sistema), RESULTADOS_DIR (ver rutas.py)
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTENEDOR="${CONTENEDOR:-$HOME/work/lab/containers/workstation_full.sif}"
if [ "$CONTENEDOR" = no ]; then py() { python3 "$@"; }
else py() { apptainer exec "$CONTENEDOR" python3 "$@"; }; fi

PASOS=("$@")
[ ${#PASOS[@]} -eq 0 ] && PASOS=(correlacion ccdf)
for p in "${PASOS[@]}"; do
    case "$p" in
        correlacion) echo "===== SII: correlacion entre avalanchas consecutivas ====="
                     py "$DIR/correlacion_avalanchas_formas.py" ;;
        ccdf)        echo "===== SIII.C: CCDF de tamaños de avalancha ====="
                     py "$DIR/ccdf_avalanchas.py" ;;
        *) echo "Paso desconocido: $p (correlacion | ccdf)"; exit 1 ;;
    esac
done
