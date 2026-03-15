# =============================================================================
# SIMULACION DE BIORREACTOR APPLIKON 3.0 L  -  MODO BATCH
# Modelos acoplados: Haldane + Pirt + Luedeking-Piret
# =============================================================================
# Organismo referencia: bacteria mesofila generica en glucosa
# Vessel: Applikon 3.0 L vidrio  |  Volumen de trabajo: 2.0 L
# Python 3.10+  |  Dependencias: numpy, scipy, matplotlib
# =============================================================================

# =============================================================================
# SECCION 1 - IMPORTACIONES
# =============================================================================
import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from typing import Literal

# Forzar UTF-8 en la salida de consola (necesario en Windows con cp1252)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ─────────────────────────────────────────────────────────────────────────────
# SECCIÓN 2 · CLASE PRINCIPAL DEL MODELO
# ─────────────────────────────────────────────────────────────────────────────

class HaldanePirtPiret:
    """
    Modelo acoplado para simulación de biorreactor batch (Applikon 3.0 L).

    Integra tres modelos clásicos de bioprocesos:
    ─────────────────────────────────────────────
    1. Haldane         — cinética de crecimiento con inhibición por sustrato
    2. Pirt            — consumo de sustrato con término de mantenimiento
    3. Luedeking-Piret — formación de producto (primario, secundario o mixto)

    Parámetros cinéticos
    ────────────────────
    mu_max  : float  — tasa máxima de crecimiento específico (h⁻¹)
    Ks      : float  — constante de semisaturación de Monod (g/L)
    Ki      : float  — constante de inhibición por sustrato (g/L)
    kd      : float  — tasa de muerte celular endógena (h⁻¹)
    Yg      : float  — rendimiento verdadero de crecimiento (g biomasa / g sustrato)
    ms      : float  — coeficiente de mantenimiento (g sustrato / g biomasa / h)
    alpha   : float  — coeficiente L-P de crecimiento asociado (g producto / g biomasa)
    beta    : float  — coeficiente L-P de mantenimiento (g producto / g biomasa / h)

    Condiciones iniciales
    ─────────────────────
    X0 : float — concentración inicial de biomasa (g/L)
    S0 : float — concentración inicial de sustrato (g/L)
    P0 : float — concentración inicial de producto (g/L)

    metabolite_type : {'primary', 'secondary', 'mixed'}
        'primary'   — producción solo acoplada al crecimiento  (β efectivo = 0)
        'secondary' — producción desacoplada del crecimiento   (α efectivo = 0)
        'mixed'     — ambos términos activos simultáneamente
    """

    # ── Parámetros default: bacterias mesófilas en glucosa ────────────────────
    def __init__(
        self,
        mu_max: float = 0.40,   # h⁻¹          — tasa máxima de crecimiento
        Ks: float     = 0.15,   # g/L           — constante de semisaturación
        Ki: float     = 8.0,    # g/L           — constante de inhibición
        kd: float     = 0.01,   # h⁻¹           — tasa de muerte celular
        Yg: float     = 0.50,   # g X / g S     — rendimiento verdadero
        ms: float     = 0.02,   # g S / g X / h — coeficiente de mantenimiento
        alpha: float  = 0.30,   # g P / g X     — coeficiente Luedeking-Piret α
        beta: float   = 0.05,   # g P / g X / h — coeficiente Luedeking-Piret β
        X0: float     = 0.05,   # g/L           — biomasa inicial
        S0: float     = 10.0,   # g/L           — sustrato inicial
        P0: float     = 0.0,    # g/L           — producto inicial
        metabolite_type: Literal['primary', 'secondary', 'mixed'] = 'mixed',
        volume: float = 2.0,    # L             — volumen de trabajo (vessel 3.0 L)
    ):
        # ── Cinética de Haldane ───────────────────────────────────────────────
        self.mu_max = mu_max
        self.Ks     = Ks
        self.Ki     = Ki
        self.kd     = kd

        # ── Modelo de Pirt ────────────────────────────────────────────────────
        self.Yg = Yg
        self.ms = ms

        # ── Modelo de Luedeking-Piret ─────────────────────────────────────────
        self.alpha = alpha
        self.beta  = beta

        # ── Condiciones iniciales ─────────────────────────────────────────────
        self.X0 = X0
        self.S0 = S0
        self.P0 = P0

        # ── Configuración del biorreactor ─────────────────────────────────────
        self.volume         = volume
        self.metabolite_type = metabolite_type

        # ── Validación básica ─────────────────────────────────────────────────
        allowed = {'primary', 'secondary', 'mixed'}
        if metabolite_type not in allowed:
            raise ValueError(
                f"metabolite_type='{metabolite_type}' no válido. "
                f"Opciones: {allowed}"
            )

    # ─────────────────────────────────────────────────────────────────────────
    def haldane(self, S: float | np.ndarray) -> float | np.ndarray:
        """
        Modelo de Haldane: tasa específica de crecimiento con inhibición
        por sustrato a alta concentración.

            μ(S) = (μmax · S) / (Ks + S + S²/Ki)

        Parámetros
        ----------
        S : float o ndarray — concentración de sustrato (g/L)

        Retorna
        -------
        mu : float o ndarray — tasa específica de crecimiento (h⁻¹)

        Notas biológicas
        ----------------
        - S → 0       : μ → 0           (limitación por sustrato)
        - S = S_opt   : μ = μ_max_real  (concentración óptima)
        - S → ∞       : μ → 0           (inhibición total por exceso)
        - S_opt = √(Ks · Ki)  [g/L]
        """
        S_arr = np.asarray(S, dtype=float)
        S_safe = np.maximum(S_arr, 0.0)   # protección numérica: S no negativo
        mu = (self.mu_max * S_safe) / (self.Ks + S_safe + S_safe**2 / self.Ki)
        # Devolver escalar si la entrada fue escalar
        return float(mu) if np.ndim(S) == 0 else mu

    # ─────────────────────────────────────────────────────────────────────────
    def s_opt(self) -> float:
        """
        Concentración de sustrato que maximiza μ(S).

            S_opt = √(Ks · Ki)   [g/L]

        Obtenida igualando dμ/dS = 0.  En este punto μ alcanza su valor
        máximo real (ligeramente menor que μmax por la curvatura del denominador).

        Retorna
        -------
        s_opt : float — concentración óptima de sustrato (g/L)
        """
        return float(np.sqrt(self.Ks * self.Ki))

    # ─────────────────────────────────────────────────────────────────────────
    def odes(self, t: float, y: list) -> list:
        """
        Sistema de ecuaciones diferenciales ordinarias (ODEs) acopladas.

        Estado del sistema:  y = [X, S, P]

        ─── ODE 1: Balance de biomasa (Haldane + muerte) ────────────────────
            dX/dt = (μ(S) − kd) · X
                    ↑ crecimiento neto por sustrato disponible
                                   ↑ pérdida por muerte celular endógena

        ─── ODE 2: Balance de sustrato (Pirt) ───────────────────────────────
            dS/dt = −(μ(S)/Yg + ms) · X
                         ↑              ↑
                    sustrato para    mantenimiento celular
                    producir biomasa (independiente de μ)

        ─── ODE 3: Formación de producto (Luedeking-Piret) ──────────────────
            dP/dt = α · (dX/dt)   +   β · X
                    ↑ acoplado         ↑ desacoplado
                    al crecimiento     (mantenimiento metabólico)

        Parámetros
        ----------
        t : float — tiempo actual (h)
        y : list  — [X (g/L), S (g/L), P (g/L)]

        Retorna
        -------
        [dX/dt, dS/dt, dP/dt]  con unidades g/L/h
        """
        X, S, P = y

        # ── Protección numérica: concentraciones no negativas ─────────────────
        X = max(X, 0.0)
        S = max(S, 0.0)

        # ── Tasa específica de crecimiento (Haldane) ──────────────────────────
        mu = self.haldane(S)     # h⁻¹

        # ── ODE 1: Biomasa ─────────────────────────────────────────────────────
        dXdt = (mu - self.kd) * X
        #        ↑ crecimiento neto     ↑ muertes celulares
        #        mu·X                  kd·X

        # ── ODE 2: Sustrato (modelo de Pirt) ──────────────────────────────────
        dSdt = -(mu / self.Yg + self.ms) * X
        #         ↑ sustrato → biomasa     ↑ mantenimiento (Pirt)
        #         requiere 1/Yg g S        ms g S por g X por hora
        #         por g X producido        independiente del crecimiento

        # ── ODE 3: Producto (Luedeking-Piret) ─────────────────────────────────
        if self.metabolite_type == 'primary':
            # Metabolito primario: producción directamente acoplada al crecimiento
            # Ejemplos: etanol (S. cerevisiae), ácido láctico (Lactobacillus)
            dPdt = self.alpha * dXdt

        elif self.metabolite_type == 'secondary':
            # Metabolito secundario: producción desacoplada del crecimiento
            # Ejemplos: penicilina (P. chrysogenum), estreptomicina
            dPdt = self.beta * X

        else:  # 'mixed'
            # Metabolito mixto: ambos mecanismos contribuyen simultáneamente
            # Ejemplos: ácido acético, ácido cítrico, 2,3-butanodiol
            dPdt = self.alpha * dXdt + self.beta * X

        return [dXdt, dSdt, dPdt]

    # ─────────────────────────────────────────────────────────────────────────
    def substrate_depletion_event(self, t: float, y: list) -> float:
        """
        Evento de parada para solve_ivp: agotamiento de sustrato.

        La integración se detiene cuando S cruza el umbral de 0.01 g/L
        en sentido descendente.  Se usa para evitar continuar la simulación
        cuando el sustrato está prácticamente agotado (S ≤ 0.01 g/L).

        Parámetros
        ----------
        t : float — tiempo actual (h)
        y : list  — estado del sistema [X, S, P]

        Retorna
        -------
        float — S − 0.01  (cero cuando S = umbral de agotamiento)
        """
        return y[1] - 0.01   # cruza cero cuando S = 0.01 g/L

    # ─────────────────────────────────────────────────────────────────────────
    def run(
        self,
        t_span: tuple = (0, 72),
        dt_out: float = 0.5,
        dense_output: bool = False,
    ):
        """
        Integra el sistema de ODEs con scipy.integrate.solve_ivp (RK45).

        Parámetros
        ----------
        t_span       : tuple — (t_inicio, t_fin) en horas
        dt_out       : float — paso de tiempo para salida (h)
        dense_output : bool  — si True, activa interpolación densa de scipy

        Retorna
        -------
        sol : OdeResult — objeto scipy con atributos .t, .y, .status, .message
                          sol.y tiene forma (3, n_pasos): filas = [X, S, P]
        """
        # Configurar evento de parada con atributos requeridos por solve_ivp.
        # Se define como función local porque los métodos bound no admiten
        # asignación de atributos arbitrarios en Python.
        def event(t, y):
            """S - 0.01: cruza cero cuando el sustrato se agota."""
            return y[1] - 0.01
        event.terminal  = True   # detiene integración al detectar el cruce
        event.direction = -1     # solo cruces descendentes (S bajando a 0.01)

        # Puntos de evaluación solicitados
        t_eval = np.arange(t_span[0], t_span[1] + dt_out, dt_out)
        t_eval = t_eval[t_eval <= t_span[1]]

        sol = solve_ivp(
            fun          = self.odes,
            t_span       = t_span,
            y0           = [self.X0, self.S0, self.P0],
            method       = 'RK45',
            t_eval       = t_eval,
            events       = event,
            dense_output = dense_output,
            rtol         = 1e-6,   # tolerancia relativa adecuada para bioprocesos
            atol         = 1e-8,   # tolerancia absoluta (g/L, muy pequeño)
        )
        return sol


# ─────────────────────────────────────────────────────────────────────────────
# SECCIÓN 3 · VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────────────────────

def plot_simulation(
    model: HaldanePirtPiret,
    sol,
    fig_path: str = 'biorreactor_simulacion.png',
    dpi: int = 300,
) -> plt.Figure:
    """
    Genera figura de 4 subplots (2×2) con los resultados de la simulación.

    ┌────────────────────────┬────────────────────────┐
    │ [0,0] Curva Haldane    │ [0,1] X y S vs tiempo  │
    │       μ(S) vs S        │       (doble eje Y)     │
    ├────────────────────────┼────────────────────────┤
    │ [1,0] P vs tiempo      │ [1,1] dS/dt y dP/dt    │
    │       (L-P type)       │       tasas instánt.   │
    └────────────────────────┴────────────────────────┘

    Parámetros
    ----------
    model    : HaldanePirtPiret — instancia con parámetros del modelo
    sol      : OdeResult        — resultado de solve_ivp (.t, .y)
    fig_path : str              — ruta de salida de la figura PNG
    dpi      : int              — resolución de exportación (puntos/pulgada)

    Retorna
    -------
    fig : matplotlib.Figure — objeto de figura para uso posterior
    """
    # ── Extraer series de tiempo ──────────────────────────────────────────────
    t = sol.t
    X = sol.y[0]
    S = sol.y[1]
    P = sol.y[2]

    # ── Recalcular tasas instantáneas en cada punto de salida ─────────────────
    mu_t = model.haldane(S)
    dXdt = (mu_t - model.kd) * X
    dSdt = -(mu_t / model.Yg + model.ms) * X

    if model.metabolite_type == 'primary':
        dPdt = model.alpha * dXdt
    elif model.metabolite_type == 'secondary':
        dPdt = model.beta * X
    else:
        dPdt = model.alpha * dXdt + model.beta * X

    # ── Paleta accesible de 7 colores (Wong 2011, Nature Methods) ────────────
    #    Sin rojo/verde puros como únicos diferenciadores
    C = {
        'azul':      '#0072B2',   # biomasa / Haldane
        'naranja':   '#E69F00',   # sustrato
        'cian':      '#56B4E9',   # dS/dt
        'vermilion': '#D55E00',   # marcadores / dP/dt
        'magenta':   '#CC79A7',   # producto P
        'amarillo':  '#F0E442',   # auxiliar
        'negro':     '#000000',   # referencias
    }

    # ── Estilo global de figura ───────────────────────────────────────────────
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 10,
        'axes.titlesize': 11,
        'axes.titleweight': 'bold',
        'axes.labelsize': 10,
        'legend.fontsize': 9,
    })

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    nombres_tipo = {
        'primary':   'Primario',
        'secondary': 'Secundario',
        'mixed':     'Mixto',
    }
    fig.suptitle(
        f'Simulación · Biorreactor Applikon 3.0 L (batch, 2.0 L de trabajo)\n'
        f'Haldane + Pirt + Luedeking-Piret  |  Tipo de metabolito: '
        f'{nombres_tipo[model.metabolite_type]}',
        fontsize=12, fontweight='bold', y=1.01
    )

    def _estilo_ax(ax):
        """Elimina marcos superior y derecho; agrega cuadrícula suave."""
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, linestyle='--', alpha=0.4, linewidth=0.7)

    # ─────────────────────────────────────────────────────────────────────────
    # SUBPLOT [0,0] · Curva de Haldane μ(S) vs S
    # ─────────────────────────────────────────────────────────────────────────
    ax00 = axes[0, 0]
    S_rng  = np.linspace(0, model.S0 * 1.25, 600)
    mu_rng = model.haldane(S_rng)

    s_opt_val  = model.s_opt()
    mu_at_sopt = model.haldane(s_opt_val)
    mu_at_ks   = model.haldane(model.Ks)

    # Curva principal
    ax00.plot(S_rng, mu_rng, color=C['azul'], lw=2.2, label='μ(S)  Haldane')

    # Línea horizontal μmax
    ax00.axhline(
        model.mu_max, color=C['negro'], lw=0.9, ls=':', alpha=0.55,
        label=f'μmax = {model.mu_max:.2f} h⁻¹'
    )

    # Marcador en S_opt (máximo real de μ)
    ax00.scatter(
        [s_opt_val], [mu_at_sopt],
        s=90, color=C['vermilion'], zorder=6,
        label=f'S_opt = {s_opt_val:.2f} g/L  (μ={mu_at_sopt:.3f} h⁻¹)'
    )
    ax00.annotate(
        f'  S_opt={s_opt_val:.2f} g/L\n  μ={mu_at_sopt:.3f} h⁻¹',
        xy=(s_opt_val, mu_at_sopt),
        xytext=(s_opt_val + model.S0 * 0.07, mu_at_sopt - 0.035),
        fontsize=8, color=C['vermilion'],
        arrowprops=dict(arrowstyle='->', color=C['vermilion'], lw=0.8)
    )

    # Marcador en Ks (referencia de semisaturación)
    ax00.scatter(
        [model.Ks], [mu_at_ks],
        s=65, marker='s', color=C['naranja'], zorder=6,
        label=f'Ks = {model.Ks} g/L  (μ={mu_at_ks:.3f} h⁻¹)'
    )

    # Zona de inhibición sombreada
    ax00.axvspan(
        s_opt_val, S_rng[-1], alpha=0.07, color=C['vermilion'],
        label='Zona de inhibición'
    )

    ax00.set_xlabel('Concentración de sustrato  S  (g/L)')
    ax00.set_ylabel('Tasa específica de crecimiento  μ  (h⁻¹)')
    ax00.set_title('Cinética de Haldane — Inhibición por sustrato')
    ax00.legend(loc='upper right', fontsize=8)
    ax00.set_xlim(left=0)
    ax00.set_ylim(bottom=0)
    _estilo_ax(ax00)

    # ─────────────────────────────────────────────────────────────────────────
    # SUBPLOT [0,1] · Biomasa X y Sustrato S vs tiempo (doble eje Y)
    # ─────────────────────────────────────────────────────────────────────────
    ax01      = axes[0, 1]
    ax01_twin = ax01.twinx()

    lX, = ax01.plot(t, X, color=C['azul'], lw=2.2,
                    label='Biomasa X (g/L)')
    lS, = ax01_twin.plot(t, S, color=C['naranja'], lw=2.2, ls='--',
                          label='Sustrato S (g/L)')

    ax01.set_xlabel('Tiempo (h)')
    ax01.set_ylabel('Biomasa X (g/L)', color=C['azul'], fontweight='bold')
    ax01_twin.set_ylabel('Sustrato S (g/L)', color=C['naranja'], fontweight='bold')
    ax01.tick_params(axis='y', colors=C['azul'])
    ax01_twin.tick_params(axis='y', colors=C['naranja'])
    ax01.set_title('Biomasa y Sustrato vs Tiempo')
    ax01.set_xlim(left=0)
    ax01.set_ylim(bottom=0)
    ax01_twin.set_ylim(bottom=0)

    # Leyenda combinada
    ax01.legend(
        [lX, lS],
        [lX.get_label(), lS.get_label()],
        loc='center right', fontsize=9
    )

    # Limpiar marcos (el twin tiene su propio spine derecho)
    ax01.spines['top'].set_visible(False)
    ax01_twin.spines['top'].set_visible(False)
    ax01.grid(True, linestyle='--', alpha=0.4, linewidth=0.7)

    # ─────────────────────────────────────────────────────────────────────────
    # SUBPLOT [1,0] · Producto P vs tiempo
    # ─────────────────────────────────────────────────────────────────────────
    ax10 = axes[1, 0]

    etiquetas_lp = {
        'primary':   'Primario — dP/dt = α · dX/dt',
        'secondary': 'Secundario — dP/dt = β · X',
        'mixed':     'Mixto — dP/dt = α · dX/dt + β · X',
    }
    ax10.plot(t, P, color=C['magenta'], lw=2.2,
              label=etiquetas_lp[model.metabolite_type])
    ax10.fill_between(t, P, alpha=0.12, color=C['magenta'])

    ax10.set_xlabel('Tiempo (h)')
    ax10.set_ylabel('Concentración de producto P (g/L)')
    ax10.set_title('Formación de Producto — Luedeking-Piret')
    ax10.legend(loc='upper left', fontsize=9)
    ax10.set_xlim(left=0)
    ax10.set_ylim(bottom=0)
    _estilo_ax(ax10)

    # ─────────────────────────────────────────────────────────────────────────
    # SUBPLOT [1,1] · Tasas instantáneas dS/dt y dP/dt vs tiempo
    # ─────────────────────────────────────────────────────────────────────────
    ax11 = axes[1, 1]

    ax11.plot(t, dSdt, color=C['cian'],      lw=2.2, label='dS/dt  (g/L/h)')
    ax11.plot(t, dPdt, color=C['vermilion'], lw=2.2, ls='-.', label='dP/dt  (g/L/h)')
    ax11.axhline(0, color=C['negro'], lw=0.8, ls=':', alpha=0.6)

    ax11.set_xlabel('Tiempo (h)')
    ax11.set_ylabel('Tasa  (g/L/h)')
    ax11.set_title('Tasas Instantáneas de Consumo y Producción')
    ax11.legend(loc='lower right', fontsize=9)
    ax11.set_xlim(left=0)
    _estilo_ax(ax11)

    # ── Exportar ──────────────────────────────────────────────────────────────
    plt.tight_layout()
    fig.savefig(fig_path, dpi=dpi, bbox_inches='tight')
    print(f"    → Figura guardada: {fig_path}  ({dpi} dpi)")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# SECCIÓN 4 · TABLA RESUMEN
# ─────────────────────────────────────────────────────────────────────────────

def print_summary_table(results: list[dict]) -> None:
    """
    Imprime tabla de texto con métricas finales de cada simulación.

    Columnas
    --------
    Tipo        : tipo de metabolito (primary / mixed / secondary)
    t_final (h) : tiempo al que finalizó la simulación (agotamiento o t_max)
    X_max (g/L) : concentración máxima de biomasa alcanzada
    S_final     : sustrato residual al final
    P_final     : producto acumulado al final
    Y_obs       : rendimiento experimental  Y_obs = ΔX / ΔS  [g X / g S]
    """
    cabecera = (
        f"{'Tipo':<12}  {'t_final (h)':>12}  {'X_max (g/L)':>12}  "
        f"{'S_final (g/L)':>14}  {'P_final (g/L)':>14}  {'Y_obs':>8}"
    )
    sep = '═' * len(cabecera)
    print(f'\n{sep}')
    print('  RESUMEN DE SIMULACIONES  —  Biorreactor Applikon 3.0 L (batch)')
    print(sep)
    print(cabecera)
    print('─' * len(cabecera))
    for r in results:
        print(
            f"{r['tipo']:<12}  {r['t_final']:>12.1f}  {r['X_max']:>12.4f}  "
            f"{r['S_final']:>14.4f}  {r['P_final']:>14.4f}  {r['Y_obs']:>8.4f}"
        )
    print(f'{sep}\n')


# ─────────────────────────────────────────────────────────────────────────────
# SECCIÓN 5 · BLOQUE PRINCIPAL  (__main__)
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    """
    Ejecuta tres simulaciones secuenciales (primario → mixto → secundario),
    genera una figura PNG para cada una y finaliza con una tabla resumen.
    """

    print('\n' + '=' * 60)
    print('  BiorreactorITM — Simulación Haldane · Pirt · Luedeking-Piret')
    print('  Applikon 3.0 L  |  Batch  |  2.0 L volumen de trabajo')
    print('=' * 60)

    tipos_metabolito = ['primary', 'mixed', 'secondary']
    resultados = []

    for tipo in tipos_metabolito:
        print(f'\n{"─"*60}')
        print(f'  Tipo de metabolito: {tipo.upper()}')
        print(f'{"─"*60}')

        # ── Instanciar y ejecutar el modelo ───────────────────────────────────
        model = HaldanePirtPiret(metabolite_type=tipo)
        sol   = model.run(t_span=(0, 72), dt_out=0.5)

        # ── Métricas de resultado ─────────────────────────────────────────────
        t_final = sol.t[-1]
        X_final = max(sol.y[0, -1], 0.0)
        S_final = max(sol.y[1, -1], 0.0)
        P_final = max(sol.y[2, -1], 0.0)
        X_max   = float(sol.y[0].max())

        # Rendimiento experimental observado (evitar división por cero)
        delta_S = model.S0 - S_final
        Y_obs   = (X_final - model.X0) / delta_S if delta_S > 1e-6 else 0.0

        print(f'  Estado solver          : {sol.message}')
        print(f'  Tiempo final           : {t_final:.1f} h')
        print(f'  Biomasa máxima X_max   : {X_max:.4f} g/L')
        print(f'  Sustrato final S       : {S_final:.4f} g/L')
        print(f'  Producto final P       : {P_final:.4f} g/L')
        print(f'  Rendimiento Y_obs      : {Y_obs:.4f} g X / g S')

        resultados.append({
            'tipo':    tipo,
            't_final': t_final,
            'X_max':   X_max,
            'S_final': S_final,
            'P_final': P_final,
            'Y_obs':   Y_obs,
        })

        # ── Visualización individual ──────────────────────────────────────────
        fig = plot_simulation(
            model    = model,
            sol      = sol,
            fig_path = f'biorreactor_simulacion_{tipo}.png',
            dpi      = 300,
        )
        plt.show()
        plt.close(fig)

    # ── Tabla comparativa final ───────────────────────────────────────────────
    print_summary_table(resultados)
    print('  Simulaciones completadas exitosamente.\n')


# =============================================================================
# NOTAS PARA EXTENSIÓN A FED-BATCH
# =============================================================================
#
# Para adaptar este script a operación fed-batch (alimentación continua
# con dilución D variable) se requieren los siguientes cambios:
#
# ── 1. NUEVO ESTADO: Volumen V(t) ─────────────────────────────────────────
#    Agregar V como cuarta variable de estado:
#        y = [X, S, P, V]
#    ODE adicional:
#        dV/dt = F_in(t) - F_out(t)
#    Para fed-batch sin extracción: F_out = 0 → dV/dt = F_in(t)
#
# ── 2. TASA DE DILUCIÓN VARIABLE ─────────────────────────────────────────
#        D(t) = F_in(t) / V(t)     [h⁻¹]
#    D cambia en el tiempo a medida que V aumenta con la alimentación.
#
# ── 3. ODEs MODIFICADAS (términos de dilución) ────────────────────────────
#        dX/dt = (μ(S) - kd - D) · X
#        dS/dt = -(μ/Yg + ms)·X + D·(S_feed - S)
#        dP/dt = (α·dXdt + β·X) - D·P
#    donde S_feed es la concentración de glucosa en la corriente de entrada.
#
# ── 4. PARÁMETROS ADICIONALES EN __init__ ────────────────────────────────
#        self.S_feed = 200.0   # g/L  concentración de glucosa en entrada
#        self.F_in_func = ...  # función F_in(t) → L/h (ver estrategias abajo)
#
# ── 5. ESTRATEGIAS DE ALIMENTACIÓN F_in(t) ───────────────────────────────
#    a) Constante      : F_in(t) = cte  →  dilución decreciente al crecer V
#    b) Exponencial    : F_in(t) = F0·exp(μ_set·t)  →  mantiene μ = μ_set
#    c) Por señal DO   : F_in proporcional al porcentaje de O₂ disuelto
#                        (requiere modelo de transferencia de oxígeno kLa)
#    d) Óptimo         : resolver problema de control óptimo (Pontryagin)
#                        minimizando tiempo o maximizando P_final
#
# ── 6. RESTRICCIÓN DE VOLUMEN DE TRABAJO ─────────────────────────────────
#    Agregar evento de parada cuando V(t) ≥ V_max (ej. 2.8 L en vessel 3.0 L)
#    para evitar desbordamiento del biorreactor Applikon.
#
# ── 7. EJEMPLO MÍNIMO DE PERFIL DE ALIMENTACIÓN ──────────────────────────
#    def F_in_exponencial(t, mu_set=0.15, F0=0.05):
#        return F0 * np.exp(mu_set * t)   # L/h
#
# =============================================================================
