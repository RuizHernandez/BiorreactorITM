import streamlit as st
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
import plotly.express as px
from scipy.integrate import odeint
import math

# FUNCION PARA PARSEAR DATOS

def parse_input(text):
    return np.array([float(x.strip()) for x in text.replace("\n", ",").split(",") if x.strip()])

# CONFIGURACION DE PAGINA
st.set_page_config(page_title="ReactorPro: Cinética Microbiana por Itan Ruiz (ITM)", layout="centered")
st.title("ReactorPro: Cinética Microbiana por Itan Ruiz (ITM)")
st.markdown("*Este programa está diseñado exclusivamente para estudios con **Pseudomonas reptilivora B-6bs***")

# INSTRUCCIONES
st.markdown("""
Ingrese sus datos de forma **vertical o separada por comas** para tiempo, biomasa, sustrato y oxígeno disuelto.
Puede ajustar manualmente la fase exponencial con los sliders. Los resultados, gráficas e interpretación se actualizan automáticamente.

La ecuación utilizada para estimar el **coeficiente de transferencia de oxígeno (kLa)** está basada en condiciones para un **biorreactor de 3 L marca Applikon Systems** con **dos propelas tipo Rushton**.
""")

# PARÁMETROS DEL BIORREACTOR
st.subheader("Parámetros del biorreactor")
rpm = st.number_input("Revoluciones por minuto (rpm)", min_value=100, max_value=1000, value=450, step=10)
vvm = st.number_input("Aireación (vvm)", min_value=0.1, max_value=5.0, value=3.0, step=0.1)
D_cm = st.number_input("Diámetro del impulsor (D) en cm", min_value=1.0, max_value=10.0, value=3.0, step=0.1)
D = D_cm / 100
V_litros = st.number_input("Volumen del biorreactor (L)", min_value=0.1, max_value=10.0, value=1.5, step=0.1)
V = V_litros / 1000

# ECUACIONES UTILIZADAS
st.subheader("📘 Ecuaciones utilizadas en el análisis cinético")

st.markdown("**1. Tasa específica de crecimiento (μ):**")
st.latex(r"\mu = \frac{\ln X_2 - \ln X_1}{t_2 - t_1}")
st.markdown("Se usa para calcular la velocidad de crecimiento celular durante la fase exponencial, a partir de dos puntos consecutivos de biomasa.")

st.markdown("**2. Tiempo de duplicación (td):**")
st.latex(r"t_d = \frac{\ln 2}{\mu}")
st.markdown("Indica el tiempo que tarda la población microbiana en duplicarse bajo condiciones exponenciales.")

st.markdown("**3. Rendimiento biomasa/sustrato (Yxs):**")
st.latex(r"Y_{X/S} = \frac{\Delta X}{\Delta S}")
st.markdown("Cuantifica la eficiencia de conversión de sustrato en biomasa.")

st.markdown("**4. Tasa específica de consumo de sustrato (qS):**")
st.latex(r"q_S = \frac{-\Delta S}{\Delta X \cdot \Delta t}")
st.markdown("Describe la velocidad con la que una unidad de biomasa consume el sustrato.")

st.markdown("**5. Coeficiente de transferencia de oxígeno (kLa):**")
st.latex(r"k_La = k \left( \frac{P_g}{V} \right)^\alpha v_s^\beta N^\gamma")
st.markdown("Permite estimar la capacidad del sistema para transferir oxígeno del gas al medio líquido, importante en procesos aerobios.")

# EXPLICACION
with st.expander("Desglose de términos en la fórmula de kLa"):
    st.markdown("""
    - **k**: constante empírica del sistema
    - **Pg/V**: potencia específica = (N³ D⁵) / V
    - **vs**: velocidad superficial de gas = Q/A
    - **N**: rpm / 60
    - **D**: impulsor en metros
    - **V**: volumen en m³
    """)

# ENTRADA DE DATOS
st.subheader("Ingrese datos experimentales")
tiempo_input = st.text_area("Tiempo (h)", "0\n4\n8\n12\n16\n20\n24\n48\n72")
biomasa_input = st.text_area("Biomasa (g/L)", "0.26567\n2.3\n5.5333\n5.6\n5.733\n5.2667\n5.4467\n5.5677\n3.43")
sustrato_input = st.text_area("Sustrato (g/L)", "10\n9.21\n8.7\n2.98\n0\n0\n0\n0\n0")
oxigeno_input = st.text_area("Oxígeno disuelto (%)", "100\n71.1\n8.0\n2.8\n4.8\n56.3\n62.2\n100\n100")

try:
    tiempo = parse_input(tiempo_input)
    biomasa = parse_input(biomasa_input)
    sustrato = parse_input(sustrato_input)
    oxigeno = parse_input(oxigeno_input)

    if len({len(tiempo), len(biomasa), len(sustrato), len(oxigeno)}) != 1:
        st.warning("Todos los vectores deben tener la misma longitud.")
    else:
        st.subheader("Selección de fase exponencial")
        inicio_exp = st.select_slider("Inicio fase exponencial (h)", options=tiempo.tolist(), value=tiempo[0])
        fin_exp = st.select_slider("Fin fase exponencial (h)", options=tiempo.tolist(), value=tiempo[min(2, len(tiempo)-1)])

        if inicio_exp >= fin_exp:
            st.warning("El tiempo de inicio debe ser menor al tiempo final.")
        else:
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

            N = rpm / 60
            Q = vvm * V / 60
            A = math.pi * (0.06 ** 2) / 4
            vs = Q / A
            Pg_V = (N ** 3) * (D ** 5) / V
            k = 0.002
            alpha = 0.7
            beta = 0.4
            gamma = 0.5
            kla_estimado = k * (Pg_V) ** alpha * vs ** beta * N ** gamma

            resultados = {
                "μ (h⁻¹)": mu,
                "Tiempo duplicación td (h)": td,
                "Yxs (g/g)": Yxs,
                "qs (g/g·h)": qs,
                "kLa estimado (h⁻¹)": kla_estimado,
                "R² fase exp.": r_value ** 2,
                "Inicio fase exp. (h)": inicio_exp,
                "Fin fase exp. (h)": fin_exp,
            }

            st.subheader("Resultados cinéticos")
            st.dataframe(pd.DataFrame(resultados.items(), columns=["Parámetro", "Valor"]))
            
            st.subheader("Interpretación de resultados")
            interp = []

            if mu > 0.4:
                interp.append("Alta velocidad especifica de crecimiento (μ > 0.4 h⁻¹)")
            elif mu > 0.1:
                interp.append("Velocidad especifica de crecimiento moderada (μ entre 0.1 y 0.4 h⁻¹)")
            else:
                interp.append("Baja velocidad especifica de crecimiento (μ ≤ 0.1 h⁻¹)")

            if Yxs > 0.6:
                interp.append("Alta eficiencia en conversión de sustrato a biomasa (Yxs > 0.6)")
            elif Yxs > 0.3:
                interp.append("Conversión moderada del sustrato (Yxs entre 0.3 y 0.6)")
            else:
                interp.append("Baja eficiencia de conversión de sustrato (Yxs < 0.3), posible estrés celular")

            if qs < 0.1:
                interp.append("Consumo de sustrato lento o eficiente en relación con la biomasa")
            elif qs < 1.0:
                interp.append("Consumo proporcional de sustrato (qS moderado)")
            else:
                interp.append("Alto consumo de sustrato, revisar posible desperdicio o metabolismo no balanceado")

            if kla_estimado < 0.01:
                interp.append("Bajo coeficiente de transferencia de oxígeno (kLa < 0.01 h⁻¹), posible limitación de oxígeno")
            else:
                interp.append("Transferencia de oxígeno adecuada (kLa > 0.01 h⁻¹)")

            duracion_exp = fin_exp - inicio_exp
            interp.append(f"La fase exponencial tuvo una duración de {duracion_exp} h (desde {inicio_exp} h hasta {fin_exp} h)")

            for i in interp:
                if any(p in i.lower() for p in ["limitación", "estrés", "alto consumo", "desperdicio", "baja eficiencia", "crecimiento bajo"]):
                    st.markdown(f"<span style='color:red;'>⚠️ {i}</span>", unsafe_allow_html=True)
                else:
                    st.markdown(f"<span style='color:green;'>✅ {i}</span>", unsafe_allow_html=True)


            # Versión con gráficas originales (matplotlib) con actualización en tiempo real
# Reemplaza la sección de "Animación de variables (Plotly)" por estas gráficas:

            st.subheader("Gráfica de biomasa")
            fig1, ax1 = plt.subplots()
            ax1.plot(tiempo, biomasa, 'o-', label='Biomasa', color='blue')
            ax1.axvspan(inicio_exp, fin_exp, color='orange', alpha=0.3, label='Fase exponencial')
            ax1.set_xlabel("Tiempo (h)")
            ax1.set_ylabel("Biomasa (g/L)")
            ax1.set_title("Crecimiento y fase exponencial detectada")
            ax1.grid(True)
            ax1.legend()
            st.pyplot(fig1)

            st.subheader("Gráfica de sustrato")
            fig2, ax2 = plt.subplots()
            ax2.plot(tiempo, sustrato, 's-', label='Sustrato', color='green')
            ax2.set_xlabel("Tiempo (h)")
            ax2.set_ylabel("Sustrato (g/L)")
            ax2.set_title("Consumo de sustrato")
            ax2.grid(True)
            ax2.legend()
            st.pyplot(fig2)

            st.subheader("Gráfica de oxígeno disuelto")
            fig3, ax3 = plt.subplots()
            ax3.plot(tiempo, oxigeno, 'd-', label='Oxígeno disuelto (%)', color='red')
            ax3.set_xlabel("Tiempo (h)")
            ax3.set_ylabel("Oxígeno (%)")
            ax3.set_title("Perfil de oxígeno disuelto")
            ax3.grid(True)
            ax3.legend()
            st.pyplot(fig3)


            st.subheader("Modelado en batch (odeint)")

            # PARÁMETROS PARA MODELADO
            # ------------------------------
            # CONFIGURACIÓN DE LA INTERFAZ
            # ------------------------------
            st.title("Simulación de Biorreactor Batch con *Pseudomonas*")
            st.write("Este modelo simula el crecimiento bacteriano, consumo de sustrato, producción de antibiótico y oxígeno disuelto usando un modelo tipo Monod-Luedeking-Piret.")

            # ------------------------------
            # MOSTRAR LAS ECUACIONES DEL MODELO
            # ------------------------------
            st.markdown("### 📘 Modelo matemático utilizado para simular el cultivo batch")
            
            st.markdown("**1. Tasa específica de crecimiento con limitación doble (Monod):**")
            st.latex(r"""
            \mu = \mu_{max} \cdot \frac{S}{K_S + S} \cdot \frac{O_2}{K_O + O_2}
            """)
            st.markdown("""
Este modelo combina la **limitación por sustrato** y por **oxígeno disuelto** para calcular la velocidad de crecimiento bacteriano de forma más realista en cultivos aerobios.

- $\\mu_{\\text{max}}$: tasa máxima de crecimiento
- $K_S$: constante de saturación del sustrato
- $K_O $: constante de saturación del oxígeno disuelto
- $S$, $O_2$: concentraciones actuales de sustrato y oxígenno disuelto

Permite simular el comportamiento de *Pseudomonas* bajo condiciones limitantes dobles (nutrientes + oxígeno).
""")
            
            st.markdown("**2. Variación de biomasa:**")
            st.latex(r"""
            \frac{dX}{dt} = \mu \cdot X
            """)
            st.markdown("Describe el aumento de la biomasa celular a lo largo del tiempo como función de la tasa de crecimiento.")
            
            st.markdown("**3. Consumo de sustrato:**")
            st.latex(r"""
            \frac{dS}{dt} = -\frac{1}{Y_{X/S}} \cdot \mu \cdot X
            """)
            st.markdown("Cuantifica la disminución del sustrato disponible, proporcional al crecimiento celular y al rendimiento biomasa/sustrato.")
            
            st.markdown("**4. Producción de metabolito secundario (antibiótico):**")
            st.latex(r"""
            \frac{dP}{dt} = \alpha \cdot \mu \cdot X + \beta \cdot X
            """)
            st.markdown("Modelo de **Luedeking-Piret**: considera producción tanto asociada al crecimiento celular (fase exponencial, controlada por **alpha**, es decir que mide cuando del producto se genera mientras la celula crece activamente. Por otra mano, la parte no asociada (fase estacionaria, controlada por **beta**. Es decir que mide cuanto del producto se genera independientemente del crecimiento celular especialmente en fase estacionaria")
            
            st.markdown("**5. Transferencia y consumo de oxígeno disuelto:**")
            st.latex(r"""
            \frac{dO_2}{dt} = k_La \cdot (O_2^* - O_2) - q_{O_2} \cdot X
            """)
            st.markdown("Esta ecuación represemnta el balance dinámico de óxigeno disuelto en medio liquído. Usada dentro de modelos de simulación batch para calcular cómo varía el oxígeno a lo largo del tiempo. Para usarse, requiere un valor de kLa, la cual puede venir de la ecuación anterior.")
            

# ------------------------------
# SLIDERS PARA PARÁMETROS
# ------------------------------
            st.sidebar.header("🔧 Parámetros del modelo")

            mu_max = st.sidebar.slider("μmax (h⁻¹)", 0.1, 1.0, 0.4)         
            Ks = st.sidebar.slider("Ks (g/L)", 0.1, 5.0, 0.5)
            Ko = st.sidebar.slider("Ko (mg/L)", 0.01, 1.0, 0.1)
            Yxs = st.sidebar.slider("Yxs (gX/gS)", 0.1, 1.0, 0.5)
            alpha = st.sidebar.slider("α (mg producto/gX)", 0.0, 1.0, 0.1)
            beta = st.sidebar.slider("β (mg producto/gX/h)", 0.0, 0.1, 0.02)
            kla = st.sidebar.slider("kLa (h⁻¹)", 0, 1, 10)
            O2_sat = st.sidebar.slider("O₂ saturado (mg/L)", 5.0, 10.0, 8.0)
            qO2 = st.sidebar.slider("qO₂ (mg O₂/gX/h)", 0.1, 2.0, 0.5)

# ------------------------------
# CONDICIONES INICIALES
# ------------------------------
            st.sidebar.header("📌 Condiciones iniciales")

            X0 = st.sidebar.number_input("X₀: Biomasa inicial (g/L)", 0.01, 10.0, 0.1)
            S0 = st.sidebar.number_input("S₀: Sustrato inicial (g/L)", 0.1, 50.0, 10.0)
            P0 = st.sidebar.number_input("P₀: Producto inicial (mg/L)", 0.0, 100.0, 0.0)
            O20 = st.sidebar.number_input("O₂₀: Oxígeno disuelto inicial (mg/L)", 0.0, 10.0, 7.5)
            tiempo_final = st.sidebar.slider("⏱ Tiempo final (h)", 12, 200, 72)

# ------------------------------
# SISTEMA DE ECUACIONES
# ------------------------------
            def modelo(y, t):
                X, S, P, O2 = y
                mu = mu_max * (S / (Ks + S)) * (O2 / (Ko + O2))
                dXdt = mu * X
                dSdt = -(1 / Yxs) * mu * X
                dPdt = alpha * mu * X + beta * X
                dO2dt = kla * (O2_sat - O2) - qO2 * X
                return [dXdt, dSdt, dPdt, dO2dt]

            # Tiempo de simulación
            t = np.linspace(0, tiempo_final, 500)
            cond_iniciales = [X0, S0, P0, O20]

            # Resolver el sistema
            sol = odeint(modelo, cond_iniciales, t)
            X, S, P, O2 = sol.T

            # ------------------------------
            # GRÁFICAS
            # ------------------------------
            st.subheader("📈 Resultados de la simulación")

            fig, ax = plt.subplots()
            ax.plot(t, X, label='Biomasa (X)', linewidth=2)
            ax.plot(t, S, label='Sustrato (S)', linewidth=2)
            ax.plot(t, P, label='Producto (P)', linewidth=2)
            ax.plot(t, O2, label='Oxígeno disuelto (O₂)', linewidth=2)
            ax.set_xlabel("Tiempo (h)")
            ax.set_ylabel("Concentración")
            ax.set_title("Simulación de biorreactor batch")
            ax.legend()
            ax.grid(True)
            st.pyplot(fig)
            
            
            

    st.markdown("""**Referencias:**

- Doran, P. M. (2012). *Bioprocess engineering principles* (2nd ed.). Academic Press.

- García-Ochoa, F., & Gómez, E. (2009). Bioreactor scale-up and oxygen transfer rate in microbial processes: An overview. *Biotechnology Advances, 27*(2), 153–176. https://doi.org/10.1016/j.biotechadv.2008.10.006

- Luedeking, R., & Piret, E. L. (1959). A kinetic study of the lactic acid fermentation. Batch process at controlled pH. *Journal of Biochemical and Microbiological Technology and Engineering, 1*(4), 393–412. https://doi.org/10.1002/jbmte.390010406

- Okpokwasili, G. C., & Okorie, B. B. (1988). Biodegradation potentials of bacteria isolated from a refinery effluent. *Journal of Aquatic Sciences, 3*(2), 93–98.

- Pirt, S. J. (1975). *Principles of microbe and cell cultivation*. Wiley.

- Shuler, M. L., & Kargi, F. (2017). *Bioprocess engineering: Basic concepts* (3rd ed.). Prentice Hall.
""")

except Exception as e:
    st.error(f"Error: {e}")




