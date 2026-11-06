#!/bin/env bash

set -euo pipefail

# ==============================================================================
# CONFIGURACIÓN DINÁMICA Y PORTABILIDAD
# ==============================================================================
# Detecta la ruta donde está guardado este script
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
EJECUTABLE_NAME="eo-silo"
EJECUTABLE="${BASE_DIR}/${EJECUTABLE_NAME}"
TEMPLATE_FILE="params_template.in"

# Ejecución local (sin SLURM)
MAX_JOBS=2          # simulaciones en paralelo (cada una usa ~1 núcleo)
LIMITE="3d"         # tiempo máximo por simulación (como el --time de SLURM)

# ==============================================================================
# MATRIZ DE SIMULACIÓN
# ==============================================================================
SHAPES_BIG=("TRIANGLE")
SHAPE_SMALL="DISK"

CHIS=(0.35 0.40)
RATIOS=(0.1)
WIDTHS=(2.0)

N_REF=2000
PHYSICS_SMALL="1.27324 0.5 0.1"

# ------------------------------------------------------------------------------
# CÁLCULO DINÁMICO DEL TOTAL
# ------------------------------------------------------------------------------
TOTAL_SIMS=$(( ${#SHAPES_BIG[@]} * ${#CHIS[@]} * ${#RATIOS[@]} * ${#WIDTHS[@]} ))

echo "----------------------------------------------------------------"
echo ">>> Preparando ${TOTAL_SIMS} simulaciones (ejecución local, sin SLURM)"
echo ">>> Hasta ${MAX_JOBS} en paralelo, límite de tiempo: ${LIMITE}."
echo "----------------------------------------------------------------"

# ------------------------------------------------------------------------------
# FUNCIONES Y PREPARACIÓN
# ------------------------------------------------------------------------------
get_sides() {
    case "$1" in
        "DISK")     echo "1" ;;
        "TRIANGLE") echo "3" ;;
        "SQUARE")   echo "4" ;;
        "PENTAGON") echo "5" ;;
        *)          echo "1" ;;
    esac
}

[ ! -f "${BASE_DIR}/${TEMPLATE_FILE}" ] && echo "ERROR: No hay template" && exit 1
[ ! -x "${EJECUTABLE}" ] && echo "ERROR: No hay ejecutable" && exit 1

RAW_LINE=$(grep -E "^[0-9]+[[:space:]]+[0-9.]+" "${BASE_DIR}/${TEMPLATE_FILE}" | head -n 1)
R_BIG=$(echo "$RAW_LINE" | awk '{print $2}')
PHYSICS_BIG=$(echo "$RAW_LINE" | awk '{for(i=4;i<=NF;i++) printf $i (i==NF?"":" ")}')
SEARCH_PATTERN=$(echo "$RAW_LINE" | sed 's/\./\\./g')

LANZADAS=0

# ==============================================================================
# BUCLES DE BARRIDO
# ==============================================================================

for S_BIG in "${SHAPES_BIG[@]}"; do
    SIDES_BIG_VAL=$(get_sides "$S_BIG")
    SIDES_SMALL_VAL=$(get_sides "$SHAPE_SMALL")
    PARENT_DIR="${BASE_DIR}/MEZCLA_${S_BIG}_${SHAPE_SMALL}"
    mkdir -p "${PARENT_DIR}"

    for RATIO in "${RATIOS[@]}"; do
        for CHI in "${CHIS[@]}"; do

            read N_BIG N_SMALL R_SMALL_CALC <<< $(awk -v nref="$N_REF" -v chi="$CHI" -v rb="$R_BIG" -v ratio_target="$RATIO" -v sb="$SIDES_BIG_VAL" -v ss="$SIDES_SMALL_VAL" '
            function get_area(sides, r) {
                if (sides == 1) return 3.14159265359 * r * r;
                else return sides * r * r * (sin(3.14159265359/sides) / cos(3.14159265359/sides));
            }
            function get_r_from_area(sides, area) {
                if (sides == 1) return sqrt(area / 3.14159265359);
                else return sqrt(area / (sides * (sin(3.14159265359/sides) / cos(3.14159265359/sides))));
            }
            BEGIN {
                area_big = get_area(sb, rb);
                area_small_target = area_big * ratio_target;
                rs_calc = get_r_from_area(ss, area_small_target);
                ratio_area = area_big / area_small_target;
                if (chi <= 0.0001) { n_big = nref; n_small = 0; }
                else if (chi >= 0.9999) { n_big = 0; n_small = int(nref * ratio_area); }
                else { n_big = int((1.0 - chi) * nref); n_small = int(chi * nref * ratio_area); }
                print n_big, n_small, rs_calc
            }')

            for W in "${WIDTHS[@]}"; do
                CHI_STR="${CHI/./p}"; RATIO_STR="${RATIO/./p}"; W_STR="${W/./p}"
                RUNID="${S_BIG}_R${RATIO_STR}_X${CHI_STR}"
                RUN_DIR="${PARENT_DIR}/${RUNID}"

                if [ -d "${RUN_DIR}" ] && ls "${RUN_DIR}"/frame* 1> /dev/null 2>&1; then
                    echo "   [saltada] ${RUNID}: ya tiene frames"
                    continue
                fi

                mkdir -p "${RUN_DIR}"
                SEED=$(od -vAn -N4 -tu4 < /dev/urandom | tr -d '[:space:]')
                cp "${BASE_DIR}/${TEMPLATE_FILE}" "${RUN_DIR}/params.in"

                sed -i -E -e "s/@RHOLE@/${W}/g" -e "s/@RUNID@/${RUNID}/g" -e "s/@SEED@/${SEED}/g" "${RUN_DIR}/params.in"

                LINEA_BIG="${N_BIG} ${R_BIG} ${SIDES_BIG_VAL} ${PHYSICS_BIG}"
                LINEA_SMALL="${N_SMALL} ${R_SMALL_CALC} ${SIDES_SMALL_VAL} ${PHYSICS_SMALL}"

                if [ "$N_SMALL" -gt "0" ] && [ "$N_BIG" -gt "0" ]; then
                    sed -i -E "s/noTipoGranos:[[:space:]]*1/noTipoGranos: 2/" "${RUN_DIR}/params.in"
                    sed -i -E "s/${SEARCH_PATTERN}/${LINEA_BIG}\n${LINEA_SMALL}/" "${RUN_DIR}/params.in"
                elif [ "$N_SMALL" -gt "0" ]; then
                    sed -i -E "s/${SEARCH_PATTERN}/${LINEA_SMALL}/" "${RUN_DIR}/params.in"
                else
                    sed -i -E "s/${SEARCH_PATTERN}/${LINEA_BIG}/" "${RUN_DIR}/params.in"
                fi

                # Si ya hay MAX_JOBS corriendo, esperar a que termine alguna
                while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
                    wait -n || true
                done

                # Lanzamiento local en segundo plano
                ( cd "${RUN_DIR}" && nohup timeout "${LIMITE}" "${EJECUTABLE}" params.in > simulacion.log 2>&1 ) &
                LANZADAS=$((LANZADAS + 1))
                echo "   [lanzada] ${RUNID} (N_big=${N_BIG}, N_small=${N_SMALL}, seed=${SEED}) -> ${RUN_DIR}/simulacion.log"
            done
        done
    done
done

echo "----------------------------------------------------------------"
echo ">>> ${LANZADAS} simulaciones lanzadas. Esperando a que terminen..."
echo ">>> (Ctrl+C solo corta este script si lo corrés en primer plano;"
echo ">>>  para dejarlo corriendo al cerrar la terminal, lanzalo con nohup)"
echo "----------------------------------------------------------------"
wait
echo ">>> Todas las simulaciones terminaron."
