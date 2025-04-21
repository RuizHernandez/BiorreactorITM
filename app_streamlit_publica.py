# ----------------------------------------------------------------------------
# ⚠️ Versión funcional de uso público – Propiedad de Itan Ruiz Hernández
# Esta aplicación se proporciona únicamente para fines académicos y demostrativos.
# El código completo y comentado se encuentra protegido con contraseña.
# No se permite su uso comercial, redistribución ni modificación sin autorización.
# Solicita acceso al código completo: itan.ruizh@itm.edu.mx
# Licencia: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/)
# ----------------------------------------------------------------------------

import streamlit as st
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

# Mensaje visible de protección
st.markdown("🔒 **Versión protegida para demostración científica. Código completo disponible bajo solicitud académica.**")

st.set_page_config(page_title="BiorreactorITM: Simulación Cinética", layout="centered")
st.title("🧪 BiorreactorITM: Simulación Cinética de *Pseudomonas*")

def parse_input(text):
    return np.array([float(x.strip()) for x in text.replace("\n", ",").split(",") if x.strip()])

# Entrada de datos
st.subheader("📥 Ingrese datos experimentales")
tiempo_input = st.text_area("⏱ Tiempo (h)", "0\n4\n8\n12\n16\n20\n24\n48\n72")
biomasa_input = st.text_area("🦠 Biomasa (g/L)", "0.26567\n2.3\n5.5333\n5.6\n5.733\n5.2667\n5.4467\n5.5677\n3.43")
sustrato_input = st.text_area("🍬 Sustrato (g/L)", "10\n9.21\n8.7\n2.98\n0\n0\n0\n0\n0")
oxigeno_input = st.text_area("🌬️ Oxígeno disuelto (%)", "100\n71.1\n8.0\n2.8\n4.8\n56.3\n62.2\n100\n100")

try:
    tiempo = parse_input(tiempo_input)
    biomasa = parse_input(biomasa_input)
    sustrato = parse_input(sustrato_input)
    oxigeno = parse_input(oxigeno_input)

    if len({len(tiempo), len(biomasa), len(sustrato), len(oxigeno)}) != 1:
        st.warning("Todos los vectores deben tener la misma longitud.")
    else:
        inicio_exp = st.select_slider("📍 Inicio fase exponencial (h)", options=tiempo.tolist(), value=tiempo[0])
        fin_exp = st.select_slider("📍 Fin fase exponencial (h)", options=tiempo.tolist(), value=tiempo[min(2, len(tiempo)-1)])

        if inicio_exp < fin_exp:
            idx_inicio = np.where(tiempo == inicio_exp)[0][0]
            idx_fin = np.where(tiempo == fin_exp)[0][0] + 1

            x_exp = tiempo[idx_inicio:idx_fin]
            y_exp = np.log(biomasa[idx_inicio:idx_fin])
            mu, intercept, r_value, *_ = linregress(x_exp, y_exp)

            td = math.log(2) / mu
            delta_x = biomasa[6] - biomasa[0]
            delta_s = sustrato[0] - sustrato[6]
            Yxs = delta_x / delta_s if delta_s != 0 else 0
            qs = mu / Yxs if Yxs != 0 else 0

            resultados = {
                "μ (h⁻¹)": mu,
                "td (h)": td,
                "Yxs (g/g)": Yxs,
                "qs (g/g·h)": qs
            }

            st.subheader("📊 Resultados cinéticos")
            st.dataframe(pd.DataFrame(resultados.items(), columns=["Parámetro", "Valor"]))

            st.subheader("📈 Gráfica de biomasa")
            fig1, ax1 = plt.subplots()
            ax1.plot(tiempo, biomasa, 'o-', color='blue')
            ax1.axvspan(inicio_exp, fin_exp, color='orange', alpha=0.3)
            ax1.set_xlabel("Tiempo (h)")
            ax1.set_ylabel("Biomasa (g/L)")
            ax1.set_title("Crecimiento microbiano")
            ax1.grid(True)
            st.pyplot(fig1)

except Exception as e:
    st.error(f"Error: {e}")
