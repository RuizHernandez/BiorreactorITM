
# BiorreactorITM: Simulador cinético para *Pseudomonas reptilivora B-6bs*

**Autor:** Itan Homero Ruiz Hernández  
**Desarrollado con:** Python + Streamlit
**Agradecimientos:** Dariana Berenice, Ana Karen, Karitza Barrios, Alexis Soria por su participación en prueba cerrada 

---

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

## ⚠️ Licencia

Este código está protegido bajo la licencia **CC BY-NC-ND 4.0**.

- No se permite el uso comercial.
- No se permite su modificación o redistribución.
- El código puede visualizarse únicamente con fines académicos.

🔗 [Ver términos de la licencia](https://creativecommons.org/licenses/by-nc-nd/4.0/)
