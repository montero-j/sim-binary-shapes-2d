# Tensor de estrés: discos y polígonos regulares

`estres.tex` documenta todo el cálculo del tensor de estrés en los atascos:

- las fuerzas de contacto de Box2D y lo que se guarda de ellas;
- la definición del tensor de cada grano y sus propiedades para discos y polígonos;
- la normal del contacto, de Box2D o inferida;
- las implementaciones en la simulación (`simulation/src/jam_output.cpp`) y en el análisis (`analysis/3_presion_y_mapas/presion_tensor_arch.py`);
- las regiones y los promedios;
- la verificación.

Los datos de la figura de comparación están en `datos/`. Los genera `analysis/utils/comparar_tensor_simulacion.py --csv`.

Para compilarlo (el PDF no se guarda en el repositorio). Funciona con `pdflatex`, `xelatex` o `lualatex`:

```bash
cd docs/estres && latexmk -pdf estres.tex && latexmk -c
```
