#!/usr/bin/env bash
# Una simulación por forma de partícula grande (DISK, TRIANGLE, SQUARE,
# PENTAGON; las chicas son discos) a chi = 0.2, con frames muy seguidos
# (para video) y cortando después de unos 4-5 atascos.
#
# Uso (desde cualquier lado):
#   bash run_scripts/run_formas_chi0p2_video.sh [carpeta_salida]
#
# Variables (con su valor por defecto):
#   CHI=0.2             fracción de área de chicas
#   RATIO=0.1           área chica / área grande
#   W=2.0               Radio_salida
#   FRAMES_CADA=20      pasos entre frames (20 x 0.001 s = 0.02 s -> 50 fps en tiempo real)
#   AVALANCHAS=5        maxAvalanchas: la simulación corta al empezar el 5.º
#                       atasco -> quedan 4-5 atascos en los frames (4 con atasco_*.dat)
#   T_MAX=600           tMax (s), por seguridad
#   MAX_JOBS=4          simulaciones en paralelo
#   SEED                semilla fija para todas (default: una al azar por forma)
#   SHAPES="DISK TRIANGLE SQUARE PENTAGON"
#   EXE                 ejecutable (default: simulation/bin/eo-silo, o run_scripts/eo-silo)
#   RENDER=1            al terminar, renderiza y arma el video (granular_render + dovideo.sh)
#   RENDER_EXE, DO_VIDEO  (default: los de renderizar_all.sh)
#   Y_MIN=-2 Y_MAX=20 FPS=50   ventana vertical y fps del video
#
# Ojo con el disco: con ~5600 granos cada frame pesa ~0.3-0.4 MB; a 0.02 s
# son ~15-20 MB por segundo simulado (incluye los 20 s de llenado).

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${BASE_DIR}/.." && pwd)"
TEMPLATE="${BASE_DIR}/params_template.in"
OUT_DIR="$(realpath -m "${1:-${BASE_DIR}/VIDEO_FORMAS_CHI_0p2}")"

CHI=${CHI:-0.2}
RATIO=${RATIO:-0.1}
W=${W:-2.0}
FRAMES_CADA=${FRAMES_CADA:-20}
AVALANCHAS=${AVALANCHAS:-5}
T_MAX=${T_MAX:-600}
MAX_JOBS=${MAX_JOBS:-4}
read -r -a SHAPES <<< "${SHAPES:-DISK TRIANGLE SQUARE PENTAGON}"
SHAPE_SMALL="DISK"
N_REF=2000
PHYSICS_SMALL="1.27324 0.5 0.1"

RENDER=${RENDER:-0}
RENDER_EXE=${RENDER_EXE:-/home/jmontero/box2d_2.4/silo-2d/render/bin/granular_render}
DO_VIDEO=${DO_VIDEO:-/home/jmontero/box2d_2.4/silo-2d/render/scripts/dovideo.sh}
Y_MIN=${Y_MIN:--2}
Y_MAX=${Y_MAX:-20}
FPS=${FPS:-50}

if [ -z "${EXE:-}" ]; then
    if [ -x "${REPO}/simulation/bin/eo-silo" ]; then EXE="${REPO}/simulation/bin/eo-silo"
    else EXE="${BASE_DIR}/eo-silo"; fi
fi
[ -f "${TEMPLATE}" ] || { echo "ERROR: no encuentro ${TEMPLATE}"; exit 1; }
[ -x "${EXE}" ] || { echo "ERROR: no encuentro el simulador en ${EXE} (make -C simulation o pasar EXE=...)"; exit 1; }
EXE="$(realpath "${EXE}")"

get_sides() {
    case "$1" in
        DISK) echo 1 ;; TRIANGLE) echo 3 ;; SQUARE) echo 4 ;; PENTAGON) echo 5 ;;
        HEXAGON) echo 6 ;; *) echo "ERROR: forma desconocida $1" >&2; exit 1 ;;
    esac
}

RAW_LINE=$(grep -E "^[0-9]+[[:space:]]+[0-9.]+" "${TEMPLATE}" | head -n 1)
R_BIG=$(echo "$RAW_LINE" | awk '{print $2}')
PHYSICS_BIG=$(echo "$RAW_LINE" | awk '{for(i=4;i<=NF;i++) printf $i (i==NF?"":" ")}')
SEARCH_PATTERN=$(echo "$RAW_LINE" | sed 's/\./\\./g')

CHI_STR="${CHI/./p}"; RATIO_STR="${RATIO/./p}"
mkdir -p "${OUT_DIR}"
echo "----------------------------------------------------------------"
echo ">>> ${#SHAPES[@]} simulaciones: chi=${CHI}, ratio=${RATIO}, W=${W}"
echo ">>> saveFrameFreq=${FRAMES_CADA}, maxAvalanchas=${AVALANCHAS}, tMax=${T_MAX}"
echo ">>> Salida: ${OUT_DIR}"
echo "----------------------------------------------------------------"

RUN_DIRS=()
for S_BIG in "${SHAPES[@]}"; do
    SB=$(get_sides "$S_BIG"); SS=$(get_sides "$SHAPE_SMALL")

    read -r N_BIG N_SMALL R_SMALL <<< "$(awk -v nref="$N_REF" -v chi="$CHI" -v rb="$R_BIG" \
        -v ratio="$RATIO" -v sb="$SB" -v ss="$SS" '
        function area(s, r) { return (s == 1) ? 3.14159265359*r*r : s*r*r*sin(3.14159265359/s)/cos(3.14159265359/s) }
        function radio(s, a) { return (s == 1) ? sqrt(a/3.14159265359) : sqrt(a/(s*sin(3.14159265359/s)/cos(3.14159265359/s))) }
        BEGIN {
            ab = area(sb, rb); as = ab * ratio
            printf "%d %d %.6f\n", int((1.0 - chi) * nref), int(chi * nref * ab / as), radio(ss, as)
        }')"

    RUNID="${S_BIG}_R${RATIO_STR}_X${CHI_STR}"
    RUN_DIR="${OUT_DIR}/${RUNID}"
    RUN_DIRS+=("${RUN_DIR}")
    if ls "${RUN_DIR}"/frames_*/frm_* >/dev/null 2>&1; then
        echo "   [saltada] ${RUNID}: ya tiene frames (borrar la carpeta para repetir)"
        continue
    fi
    mkdir -p "${RUN_DIR}"

    SEED_RUN=${SEED:-$(od -vAn -N4 -tu4 < /dev/urandom | tr -d '[:space:]')}
    P="${RUN_DIR}/params.in"
    cp "${TEMPLATE}" "$P"
    sed -i -E -e "s/@RHOLE@/${W}/g" -e "s/@RUNID@/${RUNID}/g" -e "s/@SEED@/${SEED_RUN}/g" \
              -e "s/^saveFrameFreq:.*/saveFrameFreq: ${FRAMES_CADA}/" \
              -e "s/^maxAvalanchas:.*/maxAvalanchas: ${AVALANCHAS}/" \
              -e "s/^tMax:.*/tMax: ${T_MAX}/" "$P"
    grep -q '^maxAvalanchas:' "$P" || echo "maxAvalanchas: ${AVALANCHAS}" >> "$P"
    sed -i -E "s/noTipoGranos:[[:space:]]*1/noTipoGranos: 2/" "$P"
    sed -i -E "s/${SEARCH_PATTERN}/${N_BIG} ${R_BIG} ${SB} ${PHYSICS_BIG}\n${N_SMALL} ${R_SMALL} ${SS} ${PHYSICS_SMALL}/" "$P"

    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do wait -n || true; done
    ( cd "${RUN_DIR}" && nohup "${EXE}" params.in > simulacion.log 2>&1 ) &
    echo "   [lanzada] ${RUNID} (N_big=${N_BIG}, N_small=${N_SMALL}, r_small=${R_SMALL}, seed=${SEED_RUN})"
done

echo ">>> Esperando a que terminen (seguir con: tail -f ${OUT_DIR}/*/simulacion.log)"
wait
echo ">>> Simulaciones terminadas."

for RUN_DIR in "${RUN_DIRS[@]}"; do
    n_at=$(find "${RUN_DIR}" -path "*/frames_*/atasco_[0-9]*.xy" 2>/dev/null | wc -l)
    n_fr=$(find "${RUN_DIR}" -path "*/frames_*/frm_*" 2>/dev/null | wc -l)
    echo "   $(basename "${RUN_DIR}"): ${n_fr} frames, ${n_at} snapshots de atasco, $(du -sh "${RUN_DIR}" | cut -f1)"
done

if [ "${RENDER}" = 1 ]; then
    for RUN_DIR in "${RUN_DIRS[@]}"; do
        FR=$(find "${RUN_DIR}" -maxdepth 1 -type d -name 'frames_*' | head -n 1)
        [ -n "$FR" ] || continue
        lim=$(awk '/^Radio_silo:/ {print $2 + 1}' "${RUN_DIR}/params.in")
        echo ">>> Render de $(basename "${RUN_DIR}")"
        "${RENDER_EXE}" -d "$FR" -x -"$lim" "$lim" "${Y_MIN}" "${Y_MAX}" < /dev/null
        "${DO_VIDEO}" "$FR" frm anim "${FPS}" < /dev/null
    done
fi
