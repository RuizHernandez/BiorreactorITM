
# BiorreactorITM: Simulador cinético para *Pseudomonas reptilivora B-6bs*

**Autor:** Itan Homero Ruiz Hernández  
**Desarrollado con:** Python + Streamlit
**Agradecimientos:** Dariana Berenice, Ana Karen, Karitza Barrios, Alexis Soria por su participación en prueba cerrada 

---

## 🔐 Código fuente protegido

Por motivos de protección intelectual, el código completo y comentado de esta aplicación está disponible únicamente bajo solicitud.

- El repositorio actual contiene una versión funcional para uso en Streamlit.
- El archivo `.py` original se encuentra protegido con contraseña y alojado de forma privada.

Si deseas acceso para fines académicos, envía un correo a:

📩 **D12120039@morelia.tecnm.mx**

Incluye en el mensaje:
- Tu nombre completo
- Institución
- Finalidad de uso

La contraseña será compartida únicamente a personas autorizadas.

---

## ⚠️ Licencia

Este proyecto se encuentra bajo licencia **CC BY-NC-ND 4.0**  
(No comercial / Sin obras derivadas)

🔗 [Ver términos](https://creativecommons.org/licenses/by-nc-nd/4.0/)


##  Descripción

Esta aplicación permite realizar simulaciones cinéticas de crecimiento microbiano, consumo de sustrato, dinámica del oxígeno disuelto y producción de metabolitos secundarios (como antibióticos), usando datos experimentales y modelos matemáticos clásicos aplicados a *Pseudomonas reptilivora B-6bs* en un biorreactor tipo batch.

Incluye interpretación automática, ajuste de parámetros, visualización interactiva y referencias científicas en formato APA.

---

## Acceso rápido

App en línea (Streamlit Cloud):  
https://ruizhernandez.streamlit.app

Repositorio:  
https://github.com/RuizHernandez/BiorreactorITM

---

## Características principales

- Cálculo de parámetros cinéticos:
  - Tasa específica de crecimiento (μ)
  - Tiempo de duplicación (td)
  - Rendimiento biomasa/sustrato (Yxs)
  - Tasa específica de consumo de sustrato (qS)
  - Coeficiente kLa estimado

- Simulación matemática:
  - Modelo de Monod con limitación por sustrato y oxígeno
  - Modelo de producción de metabolitos (Luedeking-Piret)
  - Balance dinámico de oxígeno disuelto

- Visualización:
  - Gráficas interactivas
  - Interpretaciones automáticas
  - Ecuaciones en formato LaTeX
  - Referencias científicas en formato APA 7

---

## 🛠Instalación local

```bash
git clone https://github.com/[tu_usuario]/BiorreactorITM.git
cd BiorreactorITM
pip install -r requirements.txt
streamlit run cinetica_app.py

---

