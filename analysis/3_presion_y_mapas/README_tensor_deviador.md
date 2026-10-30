# Tensor de esfuerzos, presión y deviador — `presion_tensor_arch.py`

Guía línea por línea del cálculo de presión `P` y esfuerzo de corte (deviador) `τ` por partícula en cada atasco. Los números de línea corresponden a `presion_tensor_arch.py` en este mismo directorio. Si el script cambia, actualizar esta tabla.

El documento `presion_tensor_arch.md` explica la teoría general. Este README sigue el código. La comparación con el material de clase (Born–Huang, círculo de Mohr, simetría por equilibrio de rotación) está en `README_tensor_vs_clases.md`.

---

## 1. Qué calcula, en una línea

Para cada partícula `p` de un atasco, a partir de las fuerzas de contacto:

$$
\sigma^{(p)}_{ij} \;=\; -\frac{1}{A_p}\sum_{c\,\in\,p} r^{c}_i\, f^{c}_j
\qquad
P = \tfrac12(\sigma_{xx}+\sigma_{yy})
\qquad
\tau = \sqrt{\left(\tfrac{\sigma_{xx}-\sigma_{yy}}{2}\right)^2+\sigma_{xy}^2}
$$

- $r^c$ es el vector desde el centro de la partícula al punto de contacto.
- $f^c$ es la fuerza que recibe la partícula en ese contacto.
- $A_p$ es el área de la partícula.

El signo menos hace que la compresión sea positiva. `τ` es el radio del círculo de Mohr, $(\sigma_1-\sigma_2)/2$: es la magnitud del deviador en 2D y no depende de la orientación de los ejes.

---

## 1b. Cálculo completo, paso a paso

Resumen de lo que hacen `presion_tensor_arch.py` (por grano y por atasco) y `comparar_presion_por_chi.py` (promedio por configuración y barras de error). Los detalles de cada función están en las secciones 5 y 6.

### Fuerza de cada contacto

Box2D da, para cada contacto entre los cuerpos $A$ y $B$, el punto de contacto $\vec p$ y los impulsos normal $f_n$ y tangencial $f_t$.

- **Normal** $\hat n$, de $A$ hacia $B$:
  - discos: dirección radial;
  - polígonos: normal a la cara de contacto;
  - paredes: normal al segmento.
  Si los dos lados dan una normal, se usa la del cuerpo en el que el punto de contacto está más adentro, como hace Box2D.
- **Tangente:** $\hat t = (-n_y,\; n_x)$.
- **Fuerza sobre $A$:** el impulso se pasa a fuerza dividiendo por el paso de tiempo $\Delta t$:

$$
\vec F_A = \frac{-f_n\,\hat n + f_t\,\hat t}{\Delta t},
\qquad
\vec F_B = -\vec F_A \quad \text{(tercera ley de Newton)}
$$

Los contactos contra pared se conservan: aportan al tensor del grano, pero la pared no tiene tensor.

### Momento de las fuerzas de cada grano

Para cada grano $i$, sumando sobre sus contactos $c$, con $\vec r^{\,c} = \vec p^{\,c} - \vec x_i$ (del centro del grano al punto de contacto):

$$
S_{\alpha\beta} = \sum_{c} r^{c}_{\alpha}\, F^{c}_{\beta}
\qquad (\alpha,\beta = x,y)
$$

### Tensor de esfuerzos

Se divide por el área $A_i$ del grano y se cambia el signo para que la compresión sea positiva. La parte cruzada se simetriza:

$$
\sigma_{xx} = -\frac{S_{xx}}{A_i},
\qquad
\sigma_{yy} = -\frac{S_{yy}}{A_i},
\qquad
\sigma_{xy} = -\frac{S_{xy} + S_{yx}}{2A_i}
$$

Con equilibrio de torques, $S_{xy} = S_{yx}$ y la simetrización no cambia nada; con la normal corregida la diferencia es despreciable (sección 8).

### Presión y deviador

**Presión**, el centro del círculo de Mohr:

$$
P_i = \frac{\sigma_{xx} + \sigma_{yy}}{2}
$$

**Deviador**, el tensor menos su parte isótropa:

$$
s = \sigma - P_i\,\mathbb{1} =
\begin{pmatrix}
\dfrac{\sigma_{xx}-\sigma_{yy}}{2} & \sigma_{xy} \\[2mm]
\sigma_{xy} & -\dfrac{\sigma_{xx}-\sigma_{yy}}{2}
\end{pmatrix}
$$

**Magnitud del deviador**, el radio del círculo de Mohr o esfuerzo de corte máximo:

$$
\tau_i = \sqrt{\left(\frac{\sigma_{xx}-\sigma_{yy}}{2}\right)^{2} + \sigma_{xy}^{2}}
= \frac{\sigma_1 - \sigma_2}{2},
\qquad
\sigma_{1,2} = P_i \pm \tau_i
$$

En unidades reducidas, divididas por la gravedad:

$$
\tilde P_i = \frac{P_i}{g},
\qquad
\tilde\tau_i = \frac{\tau_i}{g},
\qquad
g = 9.8
$$

### Regiones y especies

- **Arco:** los granos del arco que bloquea el orificio.
- **ROI:** los granos que no son del arco y tienen el centro en la franja $y^{\rm arco}_{\min} \le y \le y^{\rm arco}_{\max}$, de pared a pared.
- **Bulk:** el resto.

Cada grano es además grande (B) o chico (S), según la columna `tipo` del `.xy`.

### Promedio por región

Para cada atasco $j$ y región $\mathcal R$ (y especie, si corresponde) se guardan en `resumen_presion_por_atasco.csv` las sumas pesadas por área. Son las columnas `SumP_reducida_<region>`, `SumTau_reducida_<region>` y `Area_<region>`:

$$
\Sigma^{P}_{j} = \sum_{i\in\mathcal R} \tilde P_i A_i,
\qquad
\Sigma^{\tau}_{j} = \sum_{i\in\mathcal R} \tilde\tau_i A_i,
\qquad
\mathcal A_j = \sum_{i\in\mathcal R} A_i
$$

Para cada configuración (forma, $\chi$), el valor medio es el cociente de las sumas sobre todos los atascos, es decir, un promedio pesado por área:

$$
\langle P\rangle = \frac{\sum_j \Sigma^{P}_{j}}{\sum_j \mathcal A_j},
\qquad
\langle \tau\rangle = \frac{\sum_j \Sigma^{\tau}_{j}}{\sum_j \mathcal A_j}
$$

### Barras de error

Cada atasco es una muestra, porque los granos de un mismo atasco no son independientes. Con la media de cada atasco $\langle P\rangle_j = \Sigma^{P}_{j}/\mathcal A_j$ y pesos $w_j = \mathcal A_j$:

$$
\text{SEM} = \frac{\text{std}_w\!\left(\langle P\rangle_j\right)}{\sqrt{n_{\rm eff}}},
\qquad
n_{\rm eff} = \frac{\left(\sum_j w_j\right)^2}{\sum_j w_j^2}
$$

Para $\langle\tau\rangle$ es lo mismo.

### Observaciones

- $\langle\tau\rangle$ es el promedio de la magnitud del deviador de cada grano, no la magnitud del deviador del tensor promedio. Esta última es siempre menor o igual, porque los ejes principales de cada grano apuntan en direcciones distintas. `estres_arco_vs_roi.py` imprime las dos (`tau` y `tau_m`).
- $P$ no depende de los ejes, porque es la traza. $\sigma_{xy}$ entra solo en $\tau$, $\sigma_1$ y $\sigma_2$.
- El promedio de $\sigma_{xy}$ sobre una región simétrica da aproximadamente cero (izquierda y derecha del orificio se cancelan), por eso se grafica $\tau$ y no $\sigma_{xy}$.

---

## 2. Flujo general

```
main() ─┬─ cargar_master()              índice de simulaciones (radios, lados, chi)
        ├─ listar_atascos()/leer_lista()  pares (atasco_XXXXX.xy, atasco_contacts_XXXXX.dat)
        ├─ generar_tareas()              un atasco = una tarea
        └─ Pool → process_one_jam_task() ──┬─ parse_xy()                     partículas
                                           ├─ parse_paredes()                segmentos de paredes/piso
                                           ├─ lee contactos (pd.read_csv)
                                           ├─ detectar_arco_primario()       región "arco"
                                           └─ calcular_tensor_y_presion_vectorizado()   ← tensor, P, τ
                                                 └─ fuerzas_contacto() → normales_contacto()
           Registro.guardar_lote()   →  partes/*.parquet + resumen_presion_por_atasco.csv
           combinar_partes()         →  presion_por_particula_<Forma>.parquet
           resumen_por_configuracion() → resumen_presion_por_configuracion.csv
```

---

## 3. Entradas

### 3.1 Partículas — `parse_xy()` (líneas 35–56)

Lee `atasco_XXXXX.xy` (snapshot del atasco, escrito por el simulador):

| Formato de línea | Significado |
|---|---|
| `id 1 x y r tipo` | disco: centro `(x, y)`, radio `r` |
| `id n x1 y1 … xn yn tipo` | polígono de `n` lados: vértices absolutos |
| `id < 0` | paredes/piso: no son partículas (línea 42); las lee `parse_paredes()` |

| Líneas | Qué hace |
|---|---|
| 45–48 | Disco: `area = π r²`, `r_eff = r`, `coords = None`. |
| 50–55 | Polígono: centro = **promedio de los vértices** (línea 51); `r_eff` = distancia máxima centro–vértice (circunradio, línea 53); `area` con la fórmula del zapato (`poly_area`, líneas 20–23). |

El promedio de vértices coincide con el centro de masa solo para polígonos regulares, que es el caso de este proyecto.

### 3.1b Paredes — `parse_paredes()` (líneas 58–74)

Las líneas con `id < 0` del mismo `.xy` son el silo (`ChainShape`, `id -100`) y bordes (`EdgeShape`): `id n x1 y1 … xn yn tipo`. Se convierten en los `n−1` segmentos consecutivos de la cadena (abierta, no se cierra). Se usan para la normal de los contactos polígono–pared (sección 5.2).

### 3.2 Contactos (líneas 583–585)

Lee `atasco_contacts_XXXXX.dat` con 6 columnas:

| Columna | Significado (verificado, sección 8) |
|---|---|
| `id_a`, `id_b` | cuerpos en contacto (A = fixture A de Box2D; `< 0` = pared) |
| `px`, `py` | punto de contacto (world manifold) |
| `f_norm` | **impulso** normal del paso (se divide por `dt`, línea 434) |
| `f_tan` | **impulso** tangencial, **con signo** (convención de Box2D) |

Estos archivos los escribió el binario del cluster (su código fuente no está en el repo); el formato se verificó con el equilibrio estático de los granos.

### 3.3 Parámetros por atasco — `process_one_jam_task()` (líneas 550–640)

| Líneas | Qué hace |
|---|---|
| 556–563 | `r_nominal = r_eff · cos(π/n)`: para polígonos pasa del circunradio al **inradio**; para discos queda igual. Se usa solo para separar grandes de chicos. |
| 565–580 | Umbral grande/chico: `umbral_r = (r_large + r_small)/2`, con radios del índice maestro o, si faltan, de los radios presentes en el atasco. Con χ = 0 o monodisperso, `umbral_r = 0` (todo es "B"). |
| 587–588 | `r_salida` y `g` de `params.in` (`parse_radio_salida`, `parse_g`, líneas 77–109). `g` se toma en valor absoluto. |
| 589–591 | Detecta el arco; si no hay arco, el atasco se descarta (`"No arco detectado"`). |
| 593–594 | Llama al cálculo del tensor (sección 5), pasándole las paredes de `parse_paredes(xy_file)`. |
| 596–638 | Agrega metadatos (formas, χ, `jam_id`, `sim_path`) al resumen y a las filas por partícula. |

`dt` llega por argumento (`--dt`, default `0.001`, línea 927) y tiene que coincidir con `timeStep` de la simulación.

---

## 4. Detección del arco — `detectar_arco_primario()` (líneas 112–212)

Define qué partículas son la región `arco`.

| Líneas | Qué hace |
|---|---|
| 134–143 | Un par por cada dos partículas en contacto, sin paredes. Si un par tiene 2 puntos de contacto, queda una sola arista. |
| 145–157 | Grafo con costo de arista `10^{y_max}` del par: favorece mucho los caminos bajos. |
| 191–194 | Candidatos a extremos: partículas apoyadas en el piso (`y − r < 0.2`) a la izquierda (`x < 0`) y a la derecha (`x > 0`). |
| 202–212 | Un Dijkstra por cada extremo izquierdo; se queda con el camino de costo mínimo hacia cualquier extremo derecho. |
| 161–176, 178–189 | `purgar_cimientos`: recorta de cada punta los granos que solo son cimiento (apoyados en el piso fuera del orificio o con ≥ 2 apoyos inferiores ajenos al arco). |

---

## 5. Fuerzas de contacto y tensor

### 5.1 Normal del contacto — `normales_contacto()` (líneas 353–416)

> **Normal de Box2D cuando está disponible.** Desde la rama `tensor-en-simulacion`, el simulador escribe también `atasco_contacts_full_<jam>.dat`, con la normal real del manifold de cada fila de `atasco_contacts_<jam>.dat` (ver `simulation/README_tensor_atasco.md`).
>
> `process_one_jam_task` la lee con `leer_normales_box2d` y la pasa a `fuerzas_contacto(..., normales=...)`. La inferencia de esta sección solo se usa en dos casos:
>
> - con datos viejos, sin ese archivo;
> - si el archivo no coincide fila por fila con los contactos.
>
> La columna `normales` del resumen por atasco dice cuál se usó.
>
> Motivo: la inferencia falla en los contactos **esquina contra esquina** entre polígonos. Ahí el punto de Box2D queda a ~0.01 (la piel `b2_polygonRadius`) de un vértice de cada polígono y la normal sale girada 90°. Afecta al 2–4 % de los granos, con errores de más del 10 % en su tensor.

Devuelve, para cada contacto, la normal unitaria **de A hacia B** (convención de Box2D). Con `metodo='nuevo'` (el que se usa):

| Caso | Normal | Líneas |
|---|---|---|
| A o B es **disco** | dirección radial del disco: $(\mathbf p_c-\mathbf c_A)/\lvert\cdot\rvert$ si A es disco, $(\mathbf c_B-\mathbf p_c)/\lvert\cdot\rvert$ si B es disco. **Exacta**: en Box2D el punto de contacto con un círculo está sobre la recta normal que pasa por su centro. | 377–390 |
| polígono–polígono, polígono–pared | se calcula la normal de la superficie de cada lado (`normales_poligonos`, líneas 272–318; `normales_paredes`, líneas 321–350) junto con su **interioridad**: distancia de la proyección del punto de contacto al vértice más cercano de esa cara. Se usa el lado con mayor interioridad (línea 414): la **cara de referencia** de Box2D, donde el contacto cae en el medio de la cara y no sobre un vértice. | 391–416 |

`metodo='viejo'` (líneas 409–413) reproduce la versión anterior: siempre el lado más cercano del polígono A (o de B si A es pared). Falla cuando el contacto es un **vértice** del polígono usado: la normal sale de uno de los dos lados adyacentes en vez de la cara (o centro de disco) del otro cuerpo. Solo queda para comparar en el diagnóstico.

Casos de prueba (geometría con normal conocida):

| Caso | Normal vieja | Normal nueva |
|---|---|---|
| disco–disco | ✅ | ✅ |
| vértice de cuadrado A contra disco B | ❌ (1, 0) en vez de (0.71, 0.71) | ✅ |
| cara de A contra vértice de B | ✅ | ✅ |
| vértice de A contra cara de B | ❌ (0.71, −0.71) en vez de (0, −1) | ✅ |
| piso contra cara de cuadrado | ✅ | ✅ |
| piso contra vértice de cuadrado | ❌ (−0.71, 0.71) en vez de (0, 1) | ✅ |

### 5.2 Fuerza de cada contacto — `fuerzas_contacto()` (líneas 419–439)

```python
427–429  merge de contactos con posición/forma/vértices de A (_A) y B (_B); how='left': las paredes no se pierden
430      nx, ny = normales_contacto(df, paredes, metodo)       # de A hacia B
431      tx, ty = -ny, nx                                      # t = n rotada +90°
435      Fx_A = (-f_norm*nx + f_tan*tx) / dt
436      Fy_A = (-f_norm*ny + f_tan*ty) / dt
437–438  a_valida / b_valida: el lado es un grano (no pared)
```

- **Parte normal:** `-f_norm·n̂` es la fuerza que B ejerce sobre A. Empuja a A hacia su centro.
- **Parte tangencial:** `+f_tan·t̂` con `t̂ = (−n_y, n_x)`. En Box2D la tangente es `b2Cross(n, 1) = (n_y, −n_x) = −t̂`, y el impulso sobre B es `fn·n̂ + ft·(n_y, −n_x)`. Sobre A va el opuesto: `−fn·n̂ − ft·(n_y, −n_x) = −fn·n̂ + ft·t̂`.
- **`/ dt`:** convierte el impulso del paso en fuerza.
- La fuerza sobre B es siempre **−F_A** (3ª ley).
- `signo_tan` y `dividir_dt` (argumentos) existen solo para el diagnóstico; el script usa los valores por defecto.

> **Corregido (versión anterior: `− f_tan·t̂`).** El signo tangencial estaba invertido: con el signo actual el residuo de equilibrio es 0.0002, con el anterior era 15 (sección 8). Afectaba a `τ` siempre y a `P` solo en polígonos. **Los resultados de τ, y de P en polígonos, generados antes de esta corrección hay que recalcularlos.**

### 5.3 Tensor de cada partícula — `calcular_tensor_y_presion_vectorizado()` (líneas 442–545)

| Líneas | Qué hace |
|---|---|
| 448 | Fuerzas de contacto (sección 5.2), con las paredes del atasco. |
| 449–455 | Vectores rama $\mathbf r = \mathbf p_c - \mathbf c$ para A y para B. |
| 457–460 | $S_{ij}=\sum_c r_i F_j$: `S00 = r_x F_x`, `S11 = r_y F_y`, `S01 = r_x F_y`, `S10 = r_y F_x`. Sobre A con `F_A`; sobre B con `−F_A`. Si un lado es pared, aporta 0. |
| 462–466 | Suma las contribuciones de todos los contactos donde la partícula es A y donde es B. |
| 470 | Une con la tabla de partículas. Las partículas sin contactos quedan con `S = 0`. |

### 5.4 Presión (líneas 471–475)

```python
471  P = -0.5 * (S00 + S11) / area          # = (σxx + σyy)/2, compresión > 0
475  P_reducida = P / g
```

### 5.5 Tensor y deviador (líneas 481–492)

```python
483  σxx = -S00 / area
484  σyy = -S11 / area
485  σxy = -0.5 * (S01 + S10) / area         # se simetriza
488  τ   = sqrt( ((σxx − σyy)/2)**2 + σxy**2 )
489  τ_reducida = τ / g
490–492  σxx/g, σyy/g, σxy/g                 # se guardan para recalcular invariantes después
```

| Línea | Detalle |
|---|---|
| 481–482 | Si `area = 0` (no debería pasar), el tensor vale 0 en vez de dividir por cero. |
| 485 | **Simetrización:** con torques balanceados (equilibrio estático) `S01 ≈ S10`. La parte antisimétrica se descarta. Verificado con datos: el torque neto relativo es ≈ 0 en discos y en polígonos (sección 8). |
| 488 | Radio del círculo de Mohr = $(\sigma_1 - \sigma_2)/2$ = magnitud del deviador 2D. Es invariante ante rotación de ejes. **No** es lo mismo que $\sigma_{xy}$ solo, que depende de los ejes. |

Relaciones útiles:
- Esfuerzos principales: $\sigma_{1,2} = P \pm \tau$.
- Relación de esfuerzos (fricción movilizada): $\mu^* = \tau/P$.

Algunos trabajos llaman deviador a $q = \sigma_1-\sigma_2 = 2\tau$; conviene aclarar en el paper cuál se usa.

### 5.6 Regiones y especies (líneas 494–502)

| Líneas | Qué hace |
|---|---|
| 503–509 | ROI: granos que no son del arco con el centro en la franja `y` entre el grano más bajo y el más alto del arco, de pared a pared, más `--roi_extra` por encima (default 0; antes era 10, ~10 diámetros de grano grande). Es la misma franja que usa `force_arch_v4.py` para los contactos. Sin arco: `[0, 20]` (esos atascos se descartan igual). |
| 498–499 | `region = 'arco'` si está en el arco; si no, `'roi'` si está en la franja; si no, `'bulk'`. |
| 501–502 | `tipo = 'B'` (grande) si `r_nominal ≥ umbral_r`; si no, `'S'` (chico). |

### 5.7 Resumen por atasco (líneas 503–529)

Para cada región (`arco`, `roi`, `bulk`) y por especie (`_B`, `_S`), guarda **sumas** y no promedios, así después se pueden juntar varios atascos:

| Columna | Valor |
|---|---|
| `SumP_<reg>[_tipo]` | $\sum_p P_p A_p$ |
| `SumTau_<reg>[_tipo]` | $\sum_p \tau_p A_p$ |
| `Area_<reg>[_tipo]` | $\sum_p A_p$ |
| (+ versiones `_reducida_`, divididas por g) | |

### 5.8 Salida por partícula (líneas 533–542)

`pos_x, pos_y, shape, tipo, region, area, P, P_reducida, tau, tau_reducida, sigma_xx_reducida, sigma_yy_reducida, sigma_xy_reducida, coords`. Los vértices van serializados a JSON (línea 540) para poder dibujar los polígonos.

---

## 6. Agregación por configuración — `resumen_por_configuracion()` (líneas 881–910)

| Líneas | Qué hace |
|---|---|
| 886–888 | Usa solo atascos de lotes completos (`lotes_ok.txt`). |
| 889–898 | Agrupa por forma grande/chica, χ, r_chica y carpeta de simulación; suma las columnas `SumP_*`, `SumTau_*`, `Area_*`. |
| 899–907 | `mean_P_<reg> = ΣSumP / ΣArea` y `mean_tau_<reg> = ΣSumTau / ΣArea`: promedios **pesados por área de partícula**, sobre todas las partículas de todos los atascos de la configuración. |

> `mean_tau` es el promedio de $\lvert\tau\rvert$ de cada partícula. **No** es el τ del tensor promedio de la región. Como cada partícula tiene su propia orientación principal, $\langle\tau_p\rangle \ge \tau(\langle\sigma\rangle)$. Para el deviador macroscópico (por ejemplo $\mu^*$ de la región), primero promediar $\sigma_{xx}, \sigma_{yy}, \sigma_{xy}$ pesados por área (columnas de la salida por partícula) y después aplicar la fórmula de la línea 488.

> El área usada es la del **sólido** (suma de áreas de partículas), no la de la región. El resultado es el esfuerzo medio dentro de los granos; el esfuerzo de continuo de la región es ≈ φ veces ese valor (φ = fracción de empaquetamiento).

---

## 7. Salidas en disco

| Archivo | Contenido |
|---|---|
| `<out>/presion_por_particula_<Forma>.parquet` | una fila por partícula y atasco (sección 5.8) |
| `<out>/resumen_presion_por_atasco.csv` | sumas por región/especie, una fila por atasco (5.7) |
| `<out>/resumen_presion_por_configuracion.csv` | promedios pesados por área (6) |
| `<out>/lotes_ok.txt`, `procesados.txt`, `partes/` | control de lotes para `--reanudar` (`Registro`, líneas 782–848) |

---

## 8. Verificaciones

Todo se verifica con `analysis/utils/verificar_contactos.py`, que usa **las mismas funciones** (`fuerzas_contacto`, `normales_contacto`, `parse_paredes`) que este script:

```bash
python analysis/utils/verificar_contactos.py --sim_path <carpeta_de_una_simulacion> --jam_id 1 2 3 4 5
```

Mide, por tipo de grano (discos / polígonos / que tocan pared): el residuo del equilibrio de fuerzas |ΣF − mg|/(mg), el torque neto relativo |S01 − S10|/|S| (simetría del tensor), el % de granos con P < 0 y la relación τ/P.

### Confirmadas con datos (2026-10-02)

| Pregunta | Resultado |
|---|---|
| ¿Signo de la tangencial? | Convención de Box2D. Discos (`BARRIDO_2D_DISK_DISK/Chi_0p2_W_2p0`): residuo **0.0002** con el signo actual, **15** con el anterior. Cuadrado–disco (`BARRIDO_WIDTHS_SQUARE_DISK/Ratio_0p1_Chi_0p2_W_2p0`): 0.0008 vs 14.7. |
| ¿`f_tan` tiene signo? | Sí: ~50 % de valores negativos. |
| ¿Impulsos o fuerzas? | Impulsos: dividir por `dt` es correcto (sin dividir, el residuo es ≈ 1). |
| ¿Tensor simétrico en discos? | Sí: torque neto relativo 0.0000 (mediana y percentil 90); P < 0 en 0 % de los granos. |
| ¿Tensor simétrico en cuadrados? | **Sí, con la normal nueva.** Cuadrado–disco (`BARRIDO_WIDTHS_SQUARE_DISK/Ratio_0p1_Chi_0p2_W_2p0`, atascos 1–5, ~1500 cuadrados por atasco): torque neto relativo mediana 0.0000 y percentil 90 **0.0004**; residuo de fuerzas 0.0008; P < 0 en 0.07 % de los granos. Con la normal **vieja**: percentil 90 **0.68** y residuo **8.1** en los polígonos (0.66 de percentil 90 en los que tocan pared y 0.056 en los discos que tocan polígonos). Origen: la normal en contactos por vértice (sección 5.1). |

Sobre τ > P en ~5 % de los granos (discos): no es un error. Un contacto con fricción ya genera en el tensor del grano un esfuerzo principal levemente negativo (∝ ½[fn ± √(fn² + ft²)]), así que granos con pocos contactos y mucha fricción movilizada pueden tener τ algo mayor que P.

### Verificación de σxy y τ — `analysis/utils/verificar_tensor.py`

Prueba independiente de la parte fuera de la diagonal y del deviador, con la misma función `calcular_tensor_y_presion_vectorizado()`:

```bash
python analysis/utils/verificar_tensor.py --sintetico                          # casos con resultado conocido
python analysis/utils/verificar_tensor.py --sim_path <simulacion> --jam_id 1 2 # rota un atasco real
python analysis/utils/verificar_tensor.py --parquet <out>/presion_por_particula_Disco.parquet --incluir 'Chi_0p2_W_2p0'
```

| Prueba | Qué confirma |
|---|---|
| Disco comprimido a lo largo de u | $\sigma = (2F/A)\,u u^T$: compresión a lo largo de y = x da σxy > 0; P = τ = F/A. |
| Corte puro con tangenciales (disco y cuadrado) | σxy = −2T/A y τ = 2T/A: la tangencial entra con el signo correcto. |
| Dos discos en contacto | La 3ª ley (fuerza −F_A sobre B) da el mismo tensor en los dos granos. |
| Rotación (al azar y atasco real) | P y τ no cambian; σ rota como $R\sigma R^T$. |
| Espejo y paredes (datos) | σxy(−x) = −σxy(x); junto a las paredes, σxy > 0 a la izquierda y < 0 a la derecha (las paredes sostienen peso por fricción, efecto Janssen), con \|σxy\|/σxx ≤ fricción de la pared. |

Las pruebas sintéticas dan error relativo < 1e-6 (2026-10-07).

### Pendientes

1. ~~**Clasificación grande/chico por umbral de radio.**~~ **Resuelto:** ahora `process_one_jam_task()` usa la columna `tipo` del `.xy` (la especie exacta de cada grano) con `analysis/especies.py`: la especie de mayor área mediana es `B`, la otra `S`; con χ = 0, una sola especie o áreas iguales, todo `B`. El umbral de radio (inradio vs. radios de igual área) clasificaba mal, por ejemplo, triángulos grandes en Tri–Disco con ratio 0.4. Cambia solo la separación `_B`/`_S` (P y τ por grano y los totales por región no cambian): hay que recalcular para usar las columnas por especie.
2. **Recalcular** `presion_tensor_arch.py` para las mezclas con polígonos: la normal nueva quedó confirmada con datos (cuadrado–disco, sección 8), así que P y τ calculados con la versión anterior están mal en esas mezclas (en polígonos y también en los discos que los tocan). Disco–Disco no cambia. Opcional: repetir el diagnóstico en una mezcla con triángulos, que tienen los vértices más agudos.
3. **`simulation/src/save_tensors.cpp`** tiene el mismo signo tangencial invertido. Hoy no se llama desde ningún lado del simulador, pero si se vuelve a usar hay que corregirlo.
