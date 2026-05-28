# Simulación de Polarización Óptica

Este repositorio contiene una simulación numérica de polarización óptica usando el formalismo de vectores de Jones, matrices de polarización, retardadores de fase y parámetros de Stokes. Este incluye dos versiones del código: una en Python y otra en MATLAB.

## Descripción general

El objetivo del código es modelar la propagación y análisis de un campo óptico polarizado. Se parte de un campo complejo definido sobre una malla polar y se le asigna un estado de polarización inicial. Posteriormente, se simula el paso del campo por polarizadores lineales, placas retardadoras y analizadores de polarización.

El código permite estudiar fenómenos como:

- Ley de Malus.
- Polarización lineal horizontal y vertical.
- Polarización diagonal y antidiagonal.
- Polarización circular derecha e izquierda.
- Efecto de una placa de cuarto de onda.
- Efecto de una placa de media onda.
- Cálculo de parámetros de Stokes.
- Representación en la esfera de Poincaré.
- Comparación entre resultados simulados y datos experimentales.

## Archivos incluidos

```text
polarizacionSimulacion.py
polarizacionSimulacion.m
```

### `polarizacionSimulacion.py`

Versión en Python de la simulación. Este archivo genera figuras en formato PDF dentro de una carpeta llamada `imagenes`. Incluye comparación con datos experimentales, cálculo de incertidumbres, RMSE, parámetros de Stokes y gráficas listas para usarse en reportes científicos.

### `polarizacionSimulacion.m`

Versión en MATLAB de la simulación. Permite visualizar el campo complejo, sus componentes de polarización, la intensidad transmitida por un polarizador y el cálculo de parámetros de Stokes a partir de diferentes analizadores.

## Requisitos

### Python

Para ejecutar la versión en Python se recomienda usar Python 3.10 o superior.

Dependencias principales:

```bash
pip install numpy matplotlib scienceplots
```

El código utiliza:

```python
numpy
matplotlib
scienceplots
pathlib
```

Además, el script intenta cargar un estilo externo:

```python
aip.mplstyle
```

Por lo tanto, si no se tiene este archivo en la misma carpeta del script, se puede comentar o modificar la línea:


### MATLAB

Para ejecutar la versión de MATLAB solo se requiere una instalación estándar de MATLAB con soporte para operaciones matriciales y gráficas básicas.

## Ejecución

### Ejecutar en Python

Desde la terminal, dentro de la carpeta del repositorio:

```bash
python polarizacionSimulacion.py
```

Al ejecutarse, el código crea automáticamente la carpeta:

```text
imagenes/
```

y guarda dentro de ella las figuras generadas.

### Ejecutar en MATLAB

Abrir MATLAB, colocarse en la carpeta del repositorio y ejecutar:

```matlab
polarizacionSimulacion
```

## Figuras generadas en Python

La versión de Python genera varias gráficas relevantes para el análisis experimental:

```text
imagenes/sim_leyMalus.pdf
imagenes/sim_rmse_malus.pdf
imagenes/sim_poincare.pdf
imagenes/sim_stokes.pdf
imagenes/sim_QWP_circular.pdf
```

Estas figuras corresponden a:

- Comparación entre la Ley de Malus simulada y datos experimentales.
- RMSE acumulado entre el modelo teórico y los datos normalizados.
- Representación de estados de polarización en la esfera de Poincaré.
- Comparación de parámetros de Stokes simulados y experimentales.
- Comparación entre un polarizador sin retardador y con una placa de cuarto de onda.

## Modelo físico

El campo óptico se define inicialmente como un campo complejo de la forma:

```math
E(r,\phi) =
\left(\frac{r}{w_0}\right)^m
e^{-r^2/w_0^2}
e^{im\phi}
```

donde:

- `r` es la coordenada radial.
- `phi` es la coordenada angular.
- `w0` es el ancho característico del haz.
- `m` representa la carga topológica.

Después, al campo se le asigna una polarización mediante vectores de Jones. Los estados considerados son:

```text
H   -> Horizontal
V   -> Vertical
R   -> Circular derecha
L   -> Circular izquierda
D   -> Diagonal
LD  -> Antidiagonal
```

## Polarizadores

El paso por un polarizador lineal se modela usando una matriz de Jones dependiente del ángulo:

```math
M(\theta)=
\begin{pmatrix}
\cos^2\theta & \cos\theta\sin\theta \\
\sin\theta\cos\theta & \sin^2\theta
\end{pmatrix}
```

La intensidad transmitida se calcula como:

```math
I = \sum \left(|E_x|^2 + |E_y|^2\right)
```

Esto permite reproducir la Ley de Malus:

```math
I(\theta) = I_0 \cos^2\theta
```

## Retardadores

El código también incluye la simulación de placas retardadoras, especialmente una placa de cuarto de onda orientada a 45 grados. Esta permite transformar una polarización lineal en una polarización circular ideal.

También se incluye el análisis de una placa de media onda, la cual cambia la orientación del estado de polarización.

## Parámetros de Stokes

A partir de proyecciones sobre diferentes analizadores, se calculan los parámetros de Stokes:

```math
S_0 = I_H + I_V
```

```math
S_1 = I_H - I_V
```

```math
S_2 = I_D - I_A
```

```math
S_3 = I_R - I_L
```

donde:

- `IH` es la intensidad horizontal.
- `IV` es la intensidad vertical.
- `ID` es la intensidad diagonal.
- `IA` es la intensidad antidiagonal.
- `IR` es la intensidad circular derecha.
- `IL` es la intensidad circular izquierda.

El grado de polarización se calcula como:

```math
P =
\frac{\sqrt{S_1^2+S_2^2+S_3^2}}{S_0}
```

## Datos experimentales

La versión en Python incluye datos experimentales de intensidad en función del ángulo del polarizador. Estos datos se normalizan respecto a la potencia máxima medida y se comparan con el modelo teórico de la Ley de Malus.

También se consideran incertidumbres experimentales asociadas a:

- Incertidumbre angular.
- Incertidumbre en potencia medida.

Los valores usados en el código son:

```python
ANGLE_UNC_DEG = 1.0
POWER_UNC_MW = 0.01
```


## Notas

- La versión en Python está más orientada a generar figuras listas para un reporte científico.
- La versión en MATLAB es más directa y útil para visualizar paso a paso el comportamiento del campo y sus componentes.
- Ambos códigos usan la misma idea física: representar la polarización con vectores de Jones y analizar las intensidades transmitidas después de elementos ópticos.
- Para reproducibilidad, se recomienda mantener los datos experimentales y las constantes físicas directamente documentadas dentro del código.

## Autor

Santiago Mejía
