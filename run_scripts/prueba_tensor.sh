#!/usr/bin/env bash
# Prueba corta del tensor de estrés calculado en los atascos.
#
# Corre eo-silo con run_scripts/params_prueba_tensor.in (silo chico, 250
# cuadrados y 120 discos chicos, unos 3-4 minutos) y compara el tensor de
# cada grano calculado en la simulación con el que reconstruye
# analysis/3_presion_y_mapas/presion_tensor_arch.py, con la normal de Box2D
# y con la normal inferida de la geometría.
# Ver simulation/README_tensor_atasco.md.
#
# Uso:
#   bash run_scripts/prueba_tensor.sh [carpeta_de_trabajo]
#
# Variables:
#   EO_SILO  ejecutable (default: simulation/bin/eo-silo; compilar con
#            "cd simulation && make", o "make BOX2D=/ruta/a/box2d")
#   PYTHON   intérprete con numpy, pandas y networkx (default: python3)

set -euo pipefail

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EO_SILO="${EO_SILO:-$RAIZ/simulation/bin/eo-silo}"
PYTHON="${PYTHON:-python3}"
TRABAJO="${1:-$RAIZ/prueba_tensor}"

if [[ ! -x "$EO_SILO" ]]; then
    echo "No se encuentra el ejecutable $EO_SILO (compilar con 'cd simulation && make')." >&2
    exit 1
fi
EO_SILO="$(cd "$(dirname "$EO_SILO")" && pwd)/$(basename "$EO_SILO")"

mkdir -p "$TRABAJO"
cd "$TRABAJO"
cp "$RAIZ/run_scripts/params_prueba_tensor.in" params.in
echo "# Simulación en $TRABAJO (salida en salida.log)..."
"$EO_SILO" params.in > salida.log 2>&1
echo "# Atascos guardados: $(ls frames_prueba/atasco_tensor_*.dat 2>/dev/null | wc -l)"

R_SALIDA="$(awk '/^Radio_salida:/ {print $2}' params.in)"
for normales in auto inferidas; do
    echo
    echo "########## Análisis con --normales $normales"
    "$PYTHON" "$RAIZ/analysis/utils/comparar_tensor_simulacion.py" \
        --frames frames_prueba --r_salida "$R_SALIDA" --normales "$normales" \
        --csv "comparacion_$normales.csv"
done
