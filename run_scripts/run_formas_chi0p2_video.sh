#!/usr/bin/env bash
# Una simulación por forma de partícula grande (DISK, TRIANGLE, SQUARE,
# PENTAGON; las chicas son discos) a chi = 0.2, con frames muy seguidos
# (para video) y cortando después de unos 4-5 atascos.
#
# Pensado para correr en la PC local (WSL), sin SLURM.
#
# Uso (desde cualquier lado):
#   bash run_scripts/run_formas_chi0p2_video.sh [carpeta_salida]
#   RENDER=1 bash run_scripts/run_formas_chi0p2_video.sh      # ademas arma los videos
# Para cerrar la terminal sin cortarlo:
#   nohup bash run_scripts/run_formas_chi0p2_video.sh > videos_chi0p2.out 2>&1 &
#
# carpeta_salida: default ~/sims/videos_formas_chi0p2 (disco de WSL: es mucho
# mas rapido que /mnt/c y no lo sube OneDrive/Google Drive). No conviene
# datos_simulaciones/ del repo si es un enlace a la carpeta de Windows.
#
# Variables (con su valor por defecto):
#   CHI=0.2             fracción de área de chicas
#   RATIO=0.1           área chica / área grande
#   W=2.0               Radio_salida
#   FRAMES_CADA=20      pasos entre frames (20 x 0.001 s = 0.02 s -> 50 fps en tiempo real)
#   AVALANCHAS=5        maxAvalanchas: la simulación corta al empezar el 5.º
#                       atasco -> quedan 4-5 atascos en los frames (4 con atasco_*.dat)
#   T_MAX=600           tMax (s), por seguridad
#   MAX_JOBS            simulaciones en paralelo (default: min(4, nucleos - 1))
#   SEED                semilla fija para todas (default: una al azar por forma)
#   SHAPES="DISK TRIANGLE SQUARE PENTAGON"
#   EXE                 ejecutable (default: simulation/bin/eo-silo, o run_scripts/eo-silo)
#   CONTENEDOR          si se define (ruta a un .sif), corre eo-silo con
#                       "apptainer exec $CONTENEDOR". Si no se define y al ejecutable le
#                       faltan bibliotecas (se compilo dentro del contenedor), usa solo
#                       ~/work/lab/containers/box2d_2.4.sif
#   RENDER=1            al terminar, renderiza y arma el video (granular_render + dovideo.sh)
#   RENDER_EXE, DO_VIDEO  default: render/bin/granular_render y render/scripts/dovideo.sh
#                       del repo (carpeta render/), si no los de renderizar_all.sh
#   SIN_LLENADO=1       antes de renderizar mueve los frames del llenado (t < tFill) a
#                       <corrida>/frames_llenado/ (no se borran); 0 = video con llenado
#   Y_MIN=-2 Y_MAX=20 FPS=50   ventana vertical y fps del video (Y_MAX ~ 70: silo entero)
#
# Ojo con el disco: con ~5600 granos cada frame pesa ~0.3-0.4 MB; a 0.02 s
# son ~15-20 MB por segundo simulado (incluye los 20 s de llenado).

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${BASE_DIR}/.." && pwd)"
TEMPLATE="${BASE_DIR}/params_template.in"
OUT_DIR="$(realpath -m "${1:-${HOME}/sims/videos_formas_chi0p2}")"

CHI=${CHI:-0.2}
RATIO=${RATIO:-0.1}
W=${W:-2.0}
FRAMES_CADA=${FRAMES_CADA:-20}
AVALANCHAS=${AVALANCHAS:-5}
T_MAX=${T_MAX:-600}
NUCLEOS=$(nproc 2>/dev/null || echo 2)
MAX_JOBS=${MAX_JOBS:-$(( NUCLEOS > 4 ? 4 : (NUCLEOS > 1 ? NUCLEOS - 1 : 1) ))}
read -r -a SHAPES <<< "${SHAPES:-DISK TRIANGLE SQUARE PENTAGON}"
SHAPE_SMALL="DISK"
N_REF=2000
PHYSICS_SMALL="1.27324 0.5 0.1"

RENDER=${RENDER:-0}
SIN_LLENADO=${SIN_LLENADO:-1}
if [ -z "${RENDER_EXE:-}" ]; then
    if [ -x "${REPO}/render/bin/granular_render" ]; then RENDER_EXE="${REPO}/render/bin/granular_render"
    else RENDER_EXE=/home/jmontero/box2d_2.4/silo-2d/render/bin/granular_render; fi
fi
if [ -z "${DO_VIDEO:-}" ]; then
    if [ -f "${REPO}/render/scripts/dovideo.sh" ]; then DO_VIDEO="${REPO}/render/scripts/dovideo.sh"
    else DO_VIDEO=/home/jmontero/box2d_2.4/silo-2d/render/scripts/dovideo.sh; fi
fi
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
RUN_PREFIX=()
SIF_DEFAULT="${HOME}/work/lab/containers/box2d_2.4.sif"
if [ -z "${CONTENEDOR:-}" ] && ldd "${EXE}" 2>/dev/null | grep -q "not found"; then
    echo ">>> A ${EXE} le faltan bibliotecas fuera del contenedor:"
    ldd "${EXE}" | grep "not found" | sed 's/^/      /'
    if [ -f "${SIF_DEFAULT}" ]; then CONTENEDOR="${SIF_DEFAULT}"
    else echo "ERROR: pasar CONTENEDOR=/ruta/al/box2d_2.4.sif"; exit 1; fi
fi
if [ -n "${CONTENEDOR:-}" ]; then
    [ -f "${CONTENEDOR}" ] || { echo "ERROR: no encuentro el contenedor ${CONTENEDOR}"; exit 1; }
    RUN_PREFIX=(apptainer exec "$(realpath "${CONTENEDOR}")")
    echo ">>> eo-silo se corre dentro de ${CONTENEDOR}"
fi

if [ "${RENDER}" = 1 ]; then
    [ -x "${RENDER_EXE}" ] || { echo "ERROR: no encuentro el render en ${RENDER_EXE} (pasar RENDER_EXE=...)"; exit 1; }
    [ -f "${DO_VIDEO}" ] || { echo "ERROR: no encuentro ${DO_VIDEO} (pasar DO_VIDEO=...)"; exit 1; }
fi

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
echo ">>> Hasta ${MAX_JOBS} en paralelo (${NUCLEOS} nucleos)"
LIBRE_GB=$(df -Pk "${OUT_DIR}" | awk 'NR==2 {printf "%d", $4 / 1048576}')
echo ">>> Libre en disco: ${LIBRE_GB} GB (calculo grueso: unos 2-5 GB por simulacion con FRAMES_CADA=20)"
[ "${LIBRE_GB}" -lt $(( 3 * ${#SHAPES[@]} )) ] && echo ">>> OJO: puede no alcanzar el disco (subir FRAMES_CADA o liberar espacio)"
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
    ( cd "${RUN_DIR}" && nohup "${RUN_PREFIX[@]}" "${EXE}" params.in > simulacion.log 2>&1 ) &
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
        FR=$(find "${RUN_DIR}" -maxdepth 1 -type d -name 'frames_*' ! -name 'frames_llenado' | head -n 1)
        [ -n "$FR" ] || continue
        if [ "${SIN_LLENADO}" = 1 ]; then
            t_fill=$(awk '/^tFill:/ {print $2}' "${RUN_DIR}/params.in")
            dt=$(awk '/^timeStep:/ {print $2}' "${RUN_DIR}/params.in")
            paso_fill=$(awk -v a="$t_fill" -v b="$dt" 'BEGIN { printf "%d", a / b }')
            mkdir -p "${RUN_DIR}/frames_llenado"
            n=0
            for f in "$FR"/frm_*.xy; do
                [ -e "$f" ] || continue
                paso=$(basename "$f" .xy); paso=$((10#${paso#frm_}))
                if [ "$paso" -lt "$paso_fill" ]; then mv "$f" "${RUN_DIR}/frames_llenado/"; n=$((n + 1)); fi
            done
            echo "   $(basename "${RUN_DIR}"): ${n} frames del llenado movidos a frames_llenado/"
        fi
        lim=$(awk '/^Radio_silo:/ {print $2 + 1}' "${RUN_DIR}/params.in")
        echo ">>> Render de $(basename "${RUN_DIR}")"
        "${RENDER_EXE}" -d "$FR" -x -"$lim" "$lim" "${Y_MIN}" "${Y_MAX}" < /dev/null
        "${DO_VIDEO}" "$FR" frm anim "${FPS}" < /dev/null
    done
fi
