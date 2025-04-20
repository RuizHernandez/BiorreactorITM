
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.interpolate import make_interp_spline
import math
import pandas as pd
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
ruta_csv = filedialog.askopenfilename(
    title="Selecciona tu archivo CSV de datos experimentales",
    filetypes=[("CSV files", "*.csv")]
)

if not ruta_csv:
    raise FileNotFoundError("No se seleccionó ningún archivo.")

df = pd.read_csv(ruta_csv)
tiempo = df["tiempo"].values
biomasa = df["biomasa"].values
sustrato = df["sustrato"].values
oxigeno = df["oxigeno"].values

fase_exp = (tiempo >= 4) & (tiempo <= 12)
ln_biomasa = np.log(biomasa[fase_exp])
mu, intercept, *_ = linregress(tiempo[fase_exp], ln_biomasa)
td = math.log(2) / mu
delta_x = biomasa[6] - biomasa[0]
delta_s = sustrato[0] - sustrato[6]
Yxs = delta_x / delta_s if delta_s != 0 else 0
qs = mu / Yxs if Yxs != 0 else 0

pendientes_biomasa = np.diff(biomasa) / np.diff(tiempo)
max_growth_rate = np.max(pendientes_biomasa)
pendientes_sustrato = np.diff(sustrato) / np.diff(tiempo)
max_subs_rate = -np.min(pendientes_sustrato)
pendientes_o2 = np.diff(oxigeno) / np.diff(tiempo)
max_o2_rate = np.min(pendientes_o2)

rpm = 450
N = rpm / 60
D = 0.03
V = 0.0015
vvm = 3.0
Q = vvm * V / 60

k = 0.002
alpha = 0.7
beta = 0.4
gamma = 0.5

Pg_V = (N**3 * D**5) / V
vs = Q / (math.pi * (0.06)**2 / 4)
kla_estimado = k * (Pg_V)**alpha * vs**beta * N**gamma

tiempo_suave = np.linspace(tiempo.min(), tiempo.max(), 300)
bio_spline = make_interp_spline(tiempo, biomasa, k=3)(tiempo_suave)
sus_spline = make_interp_spline(tiempo, sustrato, k=3)(tiempo_suave)
o2_spline = make_interp_spline(tiempo, oxigeno, k=3)(tiempo_suave)

plt.figure(figsize=(18, 5))

plt.subplot(1, 3, 1)
plt.plot(tiempo, biomasa, 'o', label='Biomasa')
plt.plot(tiempo_suave, bio_spline, '-', label='Curva Suave')
plt.title("Crecimiento (Biomasa)")
plt.xlabel("Tiempo (h)")
plt.ylabel("g/L")
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(tiempo, sustrato, 's', label='Sustrato', color='green')
plt.plot(tiempo_suave, sus_spline, '-', color='green')
plt.title("Consumo de sustrato")
plt.xlabel("Tiempo (h)")
plt.ylabel("g/L")
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(tiempo, oxigeno, 'd', label='O2 disuelto (%)", color='red')
plt.plot(tiempo_suave, o2_spline, '-', color='red')
plt.title("Oxígeno disuelto")
plt.xlabel("Tiempo (h)")
plt.ylabel("%")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

resultados = {
    "Velocidad específica (μ) [h⁻¹]": mu,
    "Tiempo de duplicación (td) [h]": td,
    "Yxs [g/g]": Yxs,
    "qs [g/g·h]": qs,
    "Tasa máx. crecimiento [g/L·h]": max_growth_rate,
    "Tasa máx. consumo sustrato [g/L·h]": max_subs_rate,
    "Tasa máx. consumo O2 [%/h]": max_o2_rate,
    "kLa estimado [h⁻¹]": kla_estimado
}

df_resultados = pd.DataFrame(resultados, index=["Valor"])
print("\nResultados calculados:\n")
print(df_resultados)
