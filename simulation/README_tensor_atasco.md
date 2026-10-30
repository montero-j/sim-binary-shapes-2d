# Tensor de estrés calculado en la simulación

En cada atasco, antes de romper el arco, `eo-silo` escribe en `frames_<dirID>/`:

| Archivo | Contenido |
|---|---|
| `atasco_<jam>.xy` | Granos activos y paredes (`saveAtascoFrame`, sin cambios). |
| `atasco_contacts_<jam>.dat` | Contactos en el formato de siempre: `ID_A ID_B PosX PosY NormalImpulse TangentImpulse`, una fila por punto de manifold, impulsos con signo. Es la entrada de `analysis/3_presion_y_mapas/presion_tensor_arch.py`. |
| `atasco_tensor_<jam>.dat` | **Nuevo.** Tensor de estrés de cada grano, calculado con la normal real de Box2D. |
| `atasco_contacts_full_<jam>.dat` | **Nuevo.** Los mismos contactos, en el mismo orden, con la normal de Box2D, las fuerzas y los centros. |

Las tres salidas corresponden al mismo instante: posiciones e impulsos del último `Step`.

El código está en `src/jam_output.cpp`. El escritor de `atasco_contacts` se reconstruyó a partir del binario de producción (`run_scripts/viejo.tar.gz`), cuyo fuente no estaba en el repositorio. Se reprodujeron la cabecera, el orden de las filas, el formato numérico (`%.6e`), la inclusión de los contactos con paredes y de los puntos con impulso nulo.

## Convenciones

Box2D 2.4 (`b2ContactSolver`):

- la normal `n` del manifold apunta del cuerpo A al cuerpo B;
- la tangente es `t = b2Cross(n, 1) = (n.y, -n.x)`;
- la fuerza media durante el paso sobre B es `F_B = (normalImpulse n + tangentImpulse t) / dt`, y sobre A es `F_A = -F_B`.

Tensor de cada grano (fórmula de Love, **compresión < 0**):

```
s_ij = (1/A) sum_c F_i l_j,     l = punto de contacto - centro de masa
```

- `A` es el área que usa Box2D para la masa (`ComputeMass`). Para los polígonos es el área de los vértices, sin la piel `b2_polygonRadius = 0.01`.
- `sn` es la parte debida solo a las fuerzas normales; la tangencial es `s - sn`.
- Con polígonos, `sn` no es simétrico y la fricción aporta a la traza: el brazo `l` no es paralelo a la normal.

Columnas de `atasco_tensor_<jam>.dat`:

```
gID tipo nLados x y A m sxx sxy syx syy snxx snxy snyx snyy Fx Fy tau z z_w z_f vx vy w
```

- `Fx`, `Fy` y `tau`: fuerza y torque netos de contacto. En reposo, `F + m g = 0` y `tau = 0`.
- `z`: contactos con fuerza normal total positiva. Cuenta manifolds, no puntos: un contacto cara–cara tiene dos puntos y cuenta una vez.
- `z_w`: de esos contactos, los que son con paredes.
- `z_f`: de esos contactos, los que tienen los dos puntos activos (cara–cara).
- `vx`, `vy`, `w`: velocidades lineal y angular. Sirven para descartar granos que no están en reposo.

Relación con `presion_tensor_arch.py`, que usa `S = sum r F` con el brazo primero:

```
A sxx = S00,  A sxy = S10,  A syx = S01,  A syy = S11,   P = -(sxx + syy)/2
```

## Uso en el análisis

`analysis/3_presion_y_mapas/presion_tensor_arch.py` usa la normal de `atasco_contacts_full_<jam>.dat` cuando el archivo existe y coincide fila por fila con `atasco_contacts_<jam>.dat`. Si no, la infiere de la geometría, igual que antes, así que los datos viejos siguen funcionando. La columna `normales` del resumen por atasco indica cuál se usó: `box2d` o `inferidas (...)`.

## Comparación

```bash
python3 analysis/utils/comparar_tensor_simulacion.py --frames /ruta/frames_<dirID> --r_salida <Radio_salida>
```

El script hace cuatro comparaciones:

1. Compara, contacto por contacto, la normal que infiere el análisis con la de Box2D.
2. Compara, grano por grano, P, `sxx`, `syy` y las partes simétrica y antisimétrica de `sxy`.
3. Controla el tensor del simulador:
   - equilibrio de fuerzas con la gravedad;
   - torque nulo;
   - `A (syx - sxy) = tau`;
   - virial global contra las paredes.
4. Con `--r_salida`, calcula ⟨P⟩ pesado por área en el arco y en la ROI con los dos tensores.

Con `--normales auto` (el valor por defecto), el análisis usa las mismas normales que `presion_tensor_arch.py`. Con `--normales inferidas` se fuerza la inferencia geométrica, como con los datos viejos.

## Validación (2026-10-21)

Prueba local: silo de semiancho 5 y alto 40, con 250 cuadrados (radio 0.5) y 120 discos (radio 0.158). `Radio_salida: 1.1` y `tEsperaAtasco: 5.0`, como en producción. Se analizaron 3 atascos.

Para repetirla (unos 3–4 minutos):

```bash
cd simulation && make && cd ..
bash run_scripts/prueba_tensor.sh
```

- Los parámetros están en `run_scripts/params_prueba_tensor.in`, con la semilla fija.
- El script corre la simulación en `prueba_tensor/`, que está en el `.gitignore`.
- Después corre la comparación dos veces: con `--normales auto` y con `--normales inferidas`.
- Con `EO_SILO=/ruta/a/eo-silo` se usa otro ejecutable. Con `PYTHON=...` se usa otro intérprete, que necesita numpy, pandas y networkx.
- Los números pueden variar un poco con otro compilador o con otra compilación de Box2D, porque cambian las trayectorias. Las conclusiones no cambian.

- El formato de `atasco_contacts` coincide con el del binario de producción (cabecera, columnas, `%.6e`). Las trayectorias no son idénticas, porque cambian el compilador y la compilación de Box2D.
- Tensor del simulador:
  - equilibrio de fuerzas y torque en cada grano a ~1e-6 de m g;
  - velocidades ≈ 0;
  - `A (syx - sxy) = tau` a 4e-8;
  - virial global a 1e-10.
- `presion_tensor_arch.py` coincide grano por grano a ~3e-6 en la mediana, que es la precisión de los archivos.
- Un 2–4 % de los granos tiene errores de más del 10 %, todos por la normal inferida:
  - en los contactos polígono–polígono esquina contra esquina, la normal queda girada 90°;
  - el punto de contacto de Box2D queda a ~0.01 (la piel `b2_polygonRadius`) de un vértice de cada polígono, y desde la posición del punto no se puede saber cuál es la cara de referencia;
  - discos y paredes salen exactos.
- Las medias ⟨P⟩ por región (arco y ROI) cambiaron ≤ 0.03 %.
- Con la normal de Box2D (`--normales auto`), el análisis coincide con la simulación en **todos** los granos: error máximo 6e-5 y ninguno por encima del 1 %.
- Con `tEsperaAtasco: 1.0`, un atasco guardado 1 s después de un unjam que no reanudó el flujo no estaba en reposo (|F + m g| ~ m g). Conviene filtrar los atascos con las columnas de velocidad o con el residuo de equilibrio.
