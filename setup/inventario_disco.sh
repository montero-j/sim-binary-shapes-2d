#!/bin/bash
# Uso: bash setup/inventario_disco.sh [carpeta_extra ...]
#
# Inventario de SOLO LECTURA de todo lo que hay en el disco relacionado con
# este trabajo, para despues decidir que se queda, que va comprimido a una
# papelera etiquetada y que se borra. No mueve, no borra, no abre archivos
# (solo mira nombres, tamanos y fechas, asi que no baja de la nube los
# archivos de OneDrive que estan "solo en linea").
#
# Recorre:
#   - el home de WSL ($HOME)
#   - C:\Users\<usuario>\{OneDrive,Downloads,Desktop,Documents} y work\
#   - los otros discos de Windows montados (/mnt/d, /mnt/e, ...)
#   - las carpetas extra que se pasen como argumento
#
# Escribe el reporte en ~/inventario_disco_<fecha>.txt.
set -uo pipefail

SALIDA="$HOME/inventario_disco_$(date +%Y-%m-%d_%H%M).txt"
GRANDE_MB="${GRANDE_MB:-100}"   # umbral para "archivos grandes"

# Nombres que delatan cosas de este trabajo (carpetas o archivos)
PATRONES=(
    'sim-binary*' '*silo*' 'box2d*' 'R_2.0' 'BARRIDO_*' 'Chi_*_W_*' 'frames_*'
    'master_sim_index*' 'processed_data*' 'datos_simulaciones*' 'simulation_data*'
    'resultados_*' '*octogono*' '*Octogono*' '*.sif' 'lista_atascos*'
    'presion_por_particula*' 'heatmap*' '*avalanch*' '*_R0p*_X0p*_W*'
    '0_backup*' '*BACKUP*' 'viejo*' '*_old*' '*.tar.gz' '*.zip'
)

# Carpetas que no vale la pena recorrer
PODAR=(-name AppData -o -name 'Application Data' -o -name node_modules
       -o -name .git -o -name __pycache__ -o -name '$Recycle.Bin'
       -o -name 'System Volume Information' -o -name WindowsApps)

RAICES=("$HOME")
for u in /mnt/c/Users/*/; do
    case "$(basename "$u")" in Public|Default|'Default User'|'All Users') continue;; esac
    for d in OneDrive Downloads Desktop Documents work Escritorio Descargas Documentos; do
        [ -d "$u$d" ] && RAICES+=("${u%/}/$d")
    done
done
for m in /mnt/[d-z]; do [ -d "$m" ] && [ "$(ls -A "$m" 2>/dev/null)" ] && RAICES+=("$m"); done
RAICES+=("$@")

seccion() { printf '\n\n==================== %s ====================\n' "$*"; }
# tamano legible + fecha de ultima modificacion + ruta, de una carpeta
resumen_dir() {
    local d="$1" t f
    t=$(du -sh "$d" 2>/dev/null | cut -f1)
    f=$(find "$d" -type f -printf '%TY-%Tm-%Td\n' 2>/dev/null | sort | tail -1)
    printf '%8s  %10s  %s\n' "${t:-?}" "${f:----}" "$d"
}

{
echo "Inventario generado el $(date '+%Y-%m-%d %H:%M') en $(hostname)"
echo "Carpetas recorridas:"; printf '  %s\n' "${RAICES[@]}"
echo "Columnas: tamano  ultima_modificacion  ruta"

seccion "1. Home de WSL: primer nivel (incluye ocultos)"
for e in "$HOME"/* "$HOME"/.[!.]*; do [ -e "$e" ] && resumen_dir "$e"; done | sort -k3

seccion "2. Home de WSL: segundo nivel de las carpetas de trabajo"
for e in "$HOME"/*/; do
    [ -d "$e" ] || continue
    for s in "$e"* ; do [ -d "$s" ] && resumen_dir "$s"; done
done | sort -k3

seccion "3. Caches y entornos que se pueden regenerar"
for c in "$HOME/.cache" "$HOME/.cache/pip" "$HOME/.apptainer/cache" "$HOME/.singularity/cache" \
         "$HOME/.conda/pkgs" "$HOME/miniconda3/pkgs" "$HOME/anaconda3/pkgs" "$HOME/.local/lib" \
         "$HOME/.vscode-server" "$HOME/.npm" "$HOME/.local/share/Trash"; do
    [ -d "$c" ] && resumen_dir "$c"
done
find "${RAICES[@]}" \( "${PODAR[@]}" \) -prune -o -type d \( -name venv -o -name .venv -o -name 'env' \) \
    -print -prune 2>/dev/null | while IFS= read -r d; do resumen_dir "$d"; done

seccion "4. Clones de git (para ver duplicados)"
find "${RAICES[@]}" \( -name AppData -o -name node_modules \) -prune -o -type d -name .git -print -prune 2>/dev/null |
while IFS= read -r g; do
    r=$(dirname "$g")
    url=$(git -C "$r" config --get remote.origin.url 2>/dev/null)
    rama=$(git -C "$r" branch --show-current 2>/dev/null)
    cambios=$(git -C "$r" status --porcelain 2>/dev/null | wc -l)
    ultimo=$(git -C "$r" log -1 --format='%cd %h' --date=short 2>/dev/null)
    printf '%s\n    remoto: %s | rama: %s | archivos modificados: %s | ultimo commit: %s\n' \
        "$r" "${url:-(sin remoto)}" "${rama:-?}" "$cambios" "${ultimo:-?}"
done

seccion "5. Carpetas y archivos con nombres de este trabajo"
args=(); for p in "${PATRONES[@]}"; do args+=(-o -iname "$p"); done
find "${RAICES[@]}" \( "${PODAR[@]}" \) -prune -o \( -false "${args[@]}" \) -print -prune 2>/dev/null |
    sort -u | while IFS= read -r e; do
        if [ -d "$e" ]; then resumen_dir "$e"
        else printf '%8s  %10s  %s\n' "$(du -sh "$e" 2>/dev/null | cut -f1)" \
                    "$(date -r "$e" +%Y-%m-%d 2>/dev/null)" "$e"; fi
    done

seccion "6. Archivos de mas de ${GRANDE_MB} MB"
find "${RAICES[@]}" \( "${PODAR[@]}" \) -prune -o -type f -size +"${GRANDE_MB}"M \
    -printf '%10s  %TY-%Tm-%Td  %p\n' 2>/dev/null | sort -k3 |
    awk '{ $1 = sprintf("%7.1fG", $1/2^30); print }'

seccion "7. Posibles duplicados (mismo nombre y tamano, > ${GRANDE_MB} MB)"
find "${RAICES[@]}" \( "${PODAR[@]}" \) -prune -o -type f -size +"${GRANDE_MB}"M \
    -printf '%f|%s\t%p\n' 2>/dev/null | sort |
    awk -F'\t' '{ n[$1]++; r[$1] = r[$1] "\n    " $2 } END { for (k in n) if (n[k] > 1) print k r[k] }'

seccion "8. Enlaces simbolicos rotos"
find "${RAICES[@]}" \( "${PODAR[@]}" \) -prune -o -xtype l -printf '%p -> %l\n' 2>/dev/null

seccion "9. Basura tipica (cantidad y tamano total)"
for t in '__pycache__' '.ipynb_checkpoints' 'core' 'core.*' 'slurm-*.out' 'nohup.out' '*.o' '*~' '*.swp'; do
    lista=$(find "${RAICES[@]}" \( -name AppData -o -name node_modules -o -name .git \) -prune -o -name "$t" -print 2>/dev/null)
    [ -z "$lista" ] && continue
    n=$(printf '%s\n' "$lista" | wc -l)
    tam=$(printf '%s\n' "$lista" | tr '\n' '\0' | du -shc --files0-from=- 2>/dev/null | tail -1 | cut -f1)
    printf '%-20s %5d  %s\n' "$t" "$n" "$tam"
done

seccion "10. Espacio en disco"
df -h "$HOME" /mnt/c /mnt/[d-z] 2>/dev/null | awk '!v[$0]++'
} > "$SALIDA"

echo "Listo: $SALIDA ($(wc -l < "$SALIDA") lineas)"
echo "Revisalo por arriba (puede tener rutas personales) y pasamelo."
