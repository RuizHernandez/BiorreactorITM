
# Cinetica Biorreactor GUI

Este proyecto permite calcular parámetros cinéticos microbianos a partir de datos experimentales obtenidos en cultivos tipo batch, incluyendo:

- Velocidad específica de crecimiento (μ)
- Tiempo de duplicación (td)
- Rendimiento Yxs
- Tasa de consumo específico (qs)
- Tasa máxima de crecimiento
- Tasa máxima de consumo de sustrato
- Tasa máxima de consumo de oxígeno
- Estimación del coeficiente de transferencia de oxígeno (kLa)

## Contenido del paquete

- `cinetica_biorreactor.py`: Script principal con interfaz gráfica para seleccionar archivos.
- `datos_experimentales.csv`: Archivo de ejemplo con columnas requeridas.

## Requisitos

Instala las siguientes bibliotecas si no las tienes:

```bash
pip install numpy pandas matplotlib scipy
```

## Cómo usar

1. Ejecuta el archivo `cinetica_biorreactor.py`.
2. Se abrirá una ventana para seleccionar un archivo `.csv` con tus datos.
3. El archivo debe tener las siguientes columnas:

```
tiempo, biomasa, sustrato, oxigeno
```

4. Al cargar el archivo, el script generará gráficos de:
   - Crecimiento (biomasa)
   - Consumo de sustrato
   - Consumo de oxígeno

5. Finalmente, mostrará una tabla con los resultados calculados.
