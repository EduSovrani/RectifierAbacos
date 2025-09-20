#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Curvas do retificador monofásico com filtro capacitivo (Cap. 10)
Gera plots "simples" (sem Streamlit) das Fig. 10.9, 10.10, 10.27, 10.28, 10.29.

Equações usadas:
- Relação ωRC ↔ x=Vcmin/Vpk via Eq. (10.69).
- Integrais A1..B2 (10.158–10.161) e coeficientes normalizados (10.169–10.171).
- THD = sqrt(sum_{n>=2} Cn^2) / C1   (por padrão usa SOMENTE ÍMPARES; selecione com --harmonics)
- φ1 = atan2(B1, A1)  (graus)
- FP = cos(φ1) / sqrt(1 + THD^2)
"""

import math
import argparse

# --- Lazy imports para acelerar startup ---
import importlib
import types

class _LazyMod:
    """Carrega o módulo só no primeiro acesso a um atributo."""
    def __init__(self, modname, on_before_import=None):
        self.__dict__["_modname"] = modname
        self.__dict__["_loaded"] = None
        self.__dict__["_hook"] = on_before_import

    def _load(self):
        m = self.__dict__["_loaded"]
        if m is None:
            hook = self.__dict__["_hook"]
            if hook:
                hook()
            m = importlib.import_module(self.__dict__["_modname"])
            self.__dict__["_loaded"] = m
        return m

    def __getattr__(self, name):
        return getattr(self._load(), name)

    def __setattr__(self, name, value):
        setattr(self._load(), name, value)

    def __call__(self, *a, **kw):
        return self._load()(*a, **kw)

# Backend “Agg” por padrão (renderização offscreen) — troca para TkAgg quando você chama plt.show()
def _ensure_matplotlib_base():
    import importlib
    m = importlib.import_module("matplotlib")
    # Se ainda não configurou backend explicitamente, usa Agg (mais leve)
    try:
        # Só define se nenhum backend foi selecionado
        if not getattr(m, "rcParams", None) or m.get_backend().lower() == "agg":
            m.use("Agg")
    except Exception:
        pass

# Quando o usuário chamar plt.show(), trocamos para TkAgg (se disponível) e fazemos o import
def _before_pyplot():
    _ensure_matplotlib_base()

np = _LazyMod("numpy")                 # acessa np.* só quando necessário
pd = _LazyMod("pandas")                # acessa pd.* só quando necessário
plt = _LazyMod("matplotlib.pyplot", _before_pyplot)  # carrega pyplot na 1ª chamada/uso


# -------------------------- Núcleo matemático --------------------------

def alpha_from_x(x: float) -> float:
    # α = π/2 − asin(x)
    return math.pi/2 - math.asin(x)

def beta_from_w(w: float) -> float:
    # β = π/2 + atan(−w)
    return math.pi/2 + math.atan(-w)

def f_equation(w: float, x: float) -> float:
    """Equação implícita f(w,x)=0 que liga ωRC a x (Eq. 10.69)."""
    a = alpha_from_x(x)
    b = beta_from_w(w)
    if w <= 0.0:
        return -1e9
    # ωRC*(1 − cos α) − (β cos β)/2 − ωRC*cos β*(1 − e^{(α+β−π)/ωRC}) = 0
    return w*(1 - math.cos(a)) - 0.5*b*math.cos(b) - w*math.cos(b)*(1 - math.exp((a + b - math.pi)/w))

def solve_w_for_x(x: float, w_min: float = 1e-3, w_max: float = 2e2, max_iter: int = 200) -> float:
    """Resolve f(w,x)=0 por bisseção com varredura log para achar o bracket."""
    a, b = w_min, w_max
    fa, fb = f_equation(a, x), f_equation(b, x)
    if fa*fb > 0:
        grid = np.logspace(math.log10(w_min), math.log10(w_max), 600)
        vals = [f_equation(w, x) for w in grid]
        for i in range(len(grid)-1):
            if vals[i]*vals[i+1] < 0:
                a, b = grid[i], grid[i+1]
                fa, fb = vals[i], vals[i+1]
                break
        else:
            # Sem troca de sinal: devolve o ponto de |f| mínimo (aproximação)
            return float(grid[int(np.argmin(np.abs(vals)))])
    for _ in range(max_iter):
        m = 0.5*(a+b)
        fm = f_equation(m, x)
        if abs(fm) < 1e-10 or (b-a) < 1e-8:
            return m
        if fa*fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5*(a+b)

# Integrais em [θ3, θ2] (Eqs. 10.158–10.161)
def A1(n, t3, t2):
    if n == 1:
        # ∫ cos^2θ dθ = θ/2 + sin(2θ)/4
        return (0.5*t2 + 0.25*math.sin(2*t2)) - (0.5*t3 + 0.25*math.sin(2*t3))
    return 0.5*((math.sin((n-1)*t2)-math.sin((n-1)*t3))/(n-1) +
                (math.sin((n+1)*t2)-math.sin((n+1)*t3))/(n+1))

def A2(n, t3, t2):
    if n == 1:
        # ∫ sinθ cosθ dθ = -cos(2θ)/4
        return (-math.cos(2*t2) + math.cos(2*t3))/4.0
    return 0.5*((math.cos((n-1)*t2)-math.cos((n-1)*t3))/(n-1) -
                (math.cos((n+1)*t2)-math.cos((n+1)*t3))/(n+1))

def B1(n, t3, t2):
    if n == 1:
        # ∫ cosθ sinθ dθ = -cos(2θ)/4
        return (-math.cos(2*t2) + math.cos(2*t3))/4.0
    return 0.5*((-math.cos((n+1)*t2)+math.cos((n+1)*t3))/(n+1) +
                (-math.cos((n-1)*t2)+math.cos((n-1)*t3))/(n-1))

def B2(n, t3, t2):
    if n == 1:
        # ∫ sin^2θ dθ = θ/2 − sin(2θ)/4
        return (0.5*t2 - 0.25*math.sin(2*t2)) - (0.5*t3 - 0.25*math.sin(2*t3))
    return 0.5*((math.sin((n-1)*t2)-math.sin((n-1)*t3))/(n-1) -
                (math.sin((n+1)*t2)-math.sin((n+1)*t3))/(n+1))

def coeffs_normalized(n: int, x: float, w: float):
    """Retorna (ā_n, b̄_n, c̄_n) normalizados para a harmônica n (convenção do livro: sem 1/π)."""
    theta3 = math.asin(x)             # (10.84)
    theta2 = math.pi + math.atan(-w)  # (10.83)
    y = w                             # (10.168) y = ωCR
    a_bar = y*A1(n, theta3, theta2) + A2(n, theta3, theta2)   # (10.169)
    b_bar = y*B1(n, theta3, theta2) + B2(n, theta3, theta2)   # (10.170)
    c_bar = math.hypot(a_bar, b_bar)                          # (10.171)
    return a_bar, b_bar, c_bar

def Ic_eff_over_Vpk_R(w: float, x: float) -> float:
    """(R*Ic_rms/Vpk) da Fig. 10.10 (a partir das integrais fechadas)."""
    b = beta_from_w(w)
    theta1 = math.pi + math.asin(x)   # (10.82)
    theta2 = math.pi + math.atan(-w)  # (10.83)
    theta3 = math.asin(x)             # (10.84)
    if theta2 <= theta3:
        return float('nan')
    # ∫ cos^2 θ dθ:
    F = lambda t: 0.5*t + 0.25*math.sin(2.0*t)
    term1 = (w**2 / math.pi) * (F(theta2) - F(theta3))
    # ∫ e^{2θ/w} dθ:
    delta = theta1 - theta2
    term2 = (math.cos(b)**2 / math.pi) * (-w/2.0) * (math.exp(-2.0*delta / w) - 1.0)
    val = term1 + term2
    return math.sqrt(val) if val >= 0 else float('nan')

# -------------------------- Cálculo das curvas --------------------------

def compute_curves(x_min=0.40, x_max=0.985, n_pts=150, nmax=50, harmonics_mode="odd"):
    """
    harmonics_mode:
      - "odd"  : usa somente ímpares (3,5,7,...)  [padrão, como no livro]
      - "all"  : usa pares e ímpares (2..nmax)
      - "even" : usa somente pares (2,4,6,...)
    """
    xs = np.linspace(x_min, x_max, n_pts)
    w_vals = np.array([solve_w_for_x(float(x)) for x in xs])

    # Fig. 10.10: parâmetro (R*Ic_rms/Vpk)
    y_vals = np.array([Ic_eff_over_Vpk_R(float(w), float(x)) for x, w in zip(xs, w_vals)])

    # seleção de harmônicos conforme modo
    if harmonics_mode == "all":
        harmonics = range(2, nmax + 1)
    elif harmonics_mode == "even":
        harmonics = range(2, nmax + 1, 2)
    else:  # "odd"
        harmonics = range(3, nmax + 1, 2)

    c1_list, thd_list, phi1_deg_list, pf_list = [], [], [], []
    for x, w in zip(xs, w_vals):
        a1, b1, c1 = coeffs_normalized(1, x, w)

        # soma das potências harmônicas (n conforme seleção acima)
        c_sq = 0.0
        for n in harmonics:
            _, _, cn = coeffs_normalized(n, x, w)
            c_sq += cn * cn

        thd = math.sqrt(c_sq) / c1
        phi1_deg = math.degrees(math.atan2(a1, b1))  # φ1 = atan2(B1, A1)
        pf = math.cos(math.radians(phi1_deg)) / math.sqrt(1.0 + thd*thd)

        c1_list.append(c1)
        thd_list.append(thd)
        phi1_deg_list.append(phi1_deg)
        pf_list.append(pf)

    df = pd.DataFrame({
        "Vcmin_over_Vpk": xs,
        "omegaRC": w_vals,
        "R_Ic_eff_over_Vpk": y_vals,
        "c1_bar": c1_list,
        "THD_pu": thd_list,
        "phi1_deg": phi1_deg_list,
        "PF": pf_list
    })
    return df

# -------------------------- Plots e salvamento --------------------------

def make_plots_and_save(df: pd.DataFrame, save_prefix: str = "fig"): # type: ignore
    # Fig. 10.9
    plt.figure()
    plt.plot(df["Vcmin_over_Vpk"], df["omegaRC"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel(r"$\omega RC$")
    plt.title("Fig. 10.9 — ωRC vs Vcmin/Vpk (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_9.png", dpi=200, bbox_inches="tight")

    # Fig. 10.10
    plt.figure()
    plt.plot(df["omegaRC"], df["R_Ic_eff_over_Vpk"], linewidth=2)
    plt.xlabel(r"$\omega RC$"); plt.ylabel(r"$R\cdot I_{c,ef}/V_{pk}$")
    plt.title("Fig. 10.10 — (R·Ic_ef/Vpk) vs ωRC (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_10.png", dpi=200, bbox_inches="tight")

    # Fig. 10.27
    plt.figure()
    plt.plot(df["Vcmin_over_Vpk"], df["THD_pu"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel("THD (pu)")
    plt.title("Fig. 10.27 — THD (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_27_THD.png", dpi=200, bbox_inches="tight")

    # Fig. 10.28
    plt.figure()
    plt.plot(df["Vcmin_over_Vpk"], df["phi1_deg"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel(r"$\phi_1$ (graus)")
    plt.title("Fig. 10.28 — Ângulo de deslocamento (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_28_phi1.png", dpi=200, bbox_inches="tight")

    # Fig. 10.29
    plt.figure()
    plt.plot(df["Vcmin_over_Vpk"], df["PF"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel("FP")
    plt.title("Fig. 10.29 — Fator de potência (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_29_PF.png", dpi=200, bbox_inches="tight")

    # plt.show()

# -------------------------- Main --------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Plots simples das curvas (Cap. 10) sem Streamlit.")
    p.add_argument("--xmin", type=float, default=0.40, help="x_min = Vcmin/Vpk mínimo (default 0.50)")
    p.add_argument("--xmax", type=float, default=0.985, help="x_max = Vcmin/Vpk máximo (default 0.975)")
    p.add_argument("--pts",  type=int,   default=150,  help="Número de pontos (default 1000)")
    p.add_argument("--nmax", type=int,   default=50,   help="Maior harmônico considerado para THD (default 49)")
    p.add_argument("--harmonics", choices=["odd","all","even"], default="odd",
                   help="Quais harmônicos entram no THD/FP: odd (padrão), all, even")
    p.add_argument("--csv",  type=str,   default="curvas_cap10.csv", help="Nome do CSV de saída")
    p.add_argument("--prefix", type=str, default="fig", help="Prefixo dos PNGs salvos")
    return p.parse_args()

# if __name__ == "__main__":
#     args = parse_args()
#     if args.xmax <= args.xmin:
#         raise SystemExit("Erro: --xmax deve ser > --xmin.")
#     df = compute_curves(args.xmin, args.xmax, args.pts, args.nmax, harmonics_mode=args.harmonics)
#     df.to_csv(args.csv, index=False)
#     make_plots_and_save(df, args.prefix)

    


# ==========================  BLOCO NOVO — TRIFÁSICO  ==========================
# Convenções do trifásico:
#   x ≜ Vcmin/Vp
#   w ≜ ωRC
#
# Implementações baseadas em modelo clássico do retificador trifásico a diodos com filtro capacitivo:
# - Entre os picos (separados por Δt = (π/3)/ω), o capacitor descarrega exponencialmente em R:
#       v(t) = Vp * exp( -t/(RC) )  com  1/(RC) = ω/w
#   Logo: x = Vcmin/Vp = exp( - (π/3) / w )  ⇒  w = - (π/3) / ln(x)
#   (equivalente à Eq. 10.133 rearranjada como f(w,x)= x - exp(-π/(3w)) = 0)
#
# - Corrente no capacitor entre recargas: i_C = C dv/dt = -v/R  (a corrente na carga é v/R e a fonte está em aberto)
#   Desprezando o pulso de recarga infinitesimal (ideal), a RMS no intervalo Δt é:
#       Ic_rms^2 = (1/Δt) ∫_0^{Δt} (v/R)^2 dt
#                 = (Vp^2/R^2) * (1/Δt) * [ (1/(2α)) (1 - e^{-2αΔt}) ],   α = ω/w
#       Δt = (π/3)/ω  ⇒  αΔt = (π)/(3w)  ⇒  e^{-2αΔt} = x^2
#       1/(2α) = w/(2ω),  1/Δt = 3ω/π
#       ⇒ (R*Ic_rms/Vp) = sqrt( (3w/(2π)) * (1 - x^2) )
#   (forma funcional usada para a Fig. 10.24 em função de w e x; caso sua Eq. 10.146 inclua correções
#    pelos pulsos de recarga finitos, esta expressão pode ser ajustada para coincidir exatamente.)

def alpha_from_x_tri(x: float) -> float:
    """α = π/2 - asin(r). Clipa r para evitar problemas numéricos em |r|≈1."""
    r_clip = min(max(x, -1 + 1e-12), 1 - 1e-12)
    return math.pi/2 - math.asin(r_clip)

def beta_from_wrc_tri(wrc: float) -> float:
    """β = π/2 - atan(ωRC) (equivale a β = π/2 + atan(-ωRC))."""
    return math.pi/2 + math.atan(-wrc)

def safe_exp_tri(x: float) -> float:
    """exp com proteção para over/underflow."""
    if x > 700:   # ~1e304
        return float("inf")
    if x < -700:  # ~0.0
        return 0.0
    return math.exp(x)

def F_eq_tri(wrc: float, x: float) -> float:
    """
    Equação (10.133) rearranjada como F=0:
      ωRC(1 - cosα) - (β cosβ)/2 - ωRC cosβ [1 - e^{(α+β-π/3)/ωRC}] = 0
    """
    a = alpha_from_x_tri(x)
    b = beta_from_wrc_tri(wrc)
    term1 = wrc * (1 - math.cos(a))
    term2 = (b * math.cos(b)) / 2.0
    expo  = (a + b - math.pi/3) / max(wrc, 1e-12)   # (α+β-π/3)/(ωRC) = -γ/(ωRC)
    term3 = wrc * math.cos(b) * (1 - safe_exp_tri(expo))
    return term1 - term2 - term3

def bisection_root_tri(func, x, a, b, max_iter=200, tol=1e-10):
    """Bisseção simples para achar raiz de func(wrc, r)=0 no intervalo [a,b]."""
    fa, fb = func(a, x), func(b, x)
    if not (np.isfinite(fa) and np.isfinite(fb)) or fa * fb > 0:
        return math.nan
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        fm = func(m, x)
        if not np.isfinite(fm):
            # encolhe para o lado "saudável"
            if np.isfinite(func(a, x)):
                b = m
            else:
                a = m
            continue
        if abs(fm) < tol or (b - a) / 2 < tol:
            return m
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)

def solve_w_for_x_tri(x: float,
                        w_min: float = 1e-3,
                        w_max: float = 200.0,
                        nscan: int = 8000) -> float:
    """
    Varre ωRC em [w_min, w_max], detecta mudanças de sinal de F_eq e resolve por bisseção.
    Mantém apenas raízes com γ >= 0  (i.e., α + β <= π/3).
    Se houver mais de uma, retorna a de maior ωRC (mais próxima da fronteira γ=0).
    """
    xs = np.linspace(w_min, w_max, nscan)
    vals = [F_eq_tri(wrc, x) for wrc in xs]
    roots = []
    for i in range(len(xs) - 1):
        fa, fb = vals[i], vals[i + 1]
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0 or fa * fb < 0:
            root = bisection_root_tri(F_eq_tri, x, xs[i], xs[i + 1])
            if np.isfinite(root):
                a = alpha_from_x_tri(x)
                b = beta_from_wrc_tri(root)
                # condição física: γ = π/3 - (α+β) >= 0
                if a + b <= math.pi/3 + 1e-10:
                    roots.append(root)
    if not roots:
        return math.nan
    return max(roots)

# def Ic_eff_over_Vpk_R_tri(w: float, x: float) -> float:
#     """(R*Ic_rms/Vpk) da Fig. 10.10 (a partir das integrais fechadas)."""
#     b = beta_from_wrc_tri(w)
#     theta1 = math.pi/3.0 + math.asin(x)   # (10.82)
#     theta2 = math.pi + math.atan(-w)  # (10.83)
#     theta3 = math.asin(x)             # (10.84)
#     if theta2 <= theta3:
#         return float('nan')
#     # ∫ cos^2 θ dθ:
#     F = lambda t: 0.5*t + 0.25*math.sin(2.0*t)
#     term1 = (w**2 / math.pi) * (F(theta2) - F(theta3))
#     # ∫ e^{2θ/w} dθ:
#     delta = theta1 - theta2
#     term2 = (math.cos(b)**2 / math.pi) * (-w/2.0) * (math.exp(-2.0*delta / w) - 1.0)
#     val = term1 + term2
#     return math.sqrt(val) if val >= 0 else float('nan')

# def Ic_eff_over_Vpk_R(w: float, x: float) -> float:
#     """(R*Ic_rms/Vpk) da Fig. 10.10 (a partir das integrais fechadas)."""
#     b = beta_from_w(w)
#     theta1 = math.pi + math.asin(x)   # (10.82)
#     theta2 = math.pi + math.atan(-w)  # (10.83)
#     theta3 = math.asin(x)             # (10.84)
#     if theta2 <= theta3:
#         return float('nan')
#     # ∫ cos^2 θ dθ:
#     F = lambda t: 0.5*t + 0.25*math.sin(2.0*t)
#     term1 = (w**2 / math.pi) * (F(theta2) - F(theta3))
#     # ∫ e^{2θ/w} dθ:
#     delta = theta1 - theta2
#     term2 = (math.cos(b)**2 / math.pi) * (-w/2.0) * (-math.exp(2.0*delta / w) - 1.0)
#     val = term1 + term2
#     return math.sqrt(val) if val >= 0 else float('nan')

def Ic_eff_over_Vpk_R_tri(w: float, x: float) -> float:
    b = beta_from_wrc_tri(w)
    x = max(-1.0, min(1.0, x))
    theta1 = math.pi/3.0 + math.asin(x)
    theta2 = math.pi - math.atan(w)   # = π + atan(-w)
    theta3 = math.asin(x)
    if theta2 <= theta3:
        return float('nan')

    # ∫ cos^2 θ dθ
    F = lambda t: 0.5*t + 0.25*math.sin(2.0*t)
    term1 = (3.0/math.pi) * (w**2) * (F(theta2) - F(theta3))

    # termo exponencial (sua forma, estável com expm1)
    delta = theta1 - theta2   # < 0
    term2 = (3.0/math.pi) * (math.cos(b)**2) * (-w/2.0) * math.expm1(-2.0*delta / w)

    val = term1 + term2
    return math.sqrt(val) if val >= 0.0 else float('nan')

# ----------------- helpers só pra legibilidade -----------------
def _S(x): return math.sin(x)
def _C(x): return math.cos(x)

# ================== TRECHO "a" — limites [θ3, θ2] ==================

def A1a_tri(n, t3, t2):
    """∫_{t3}^{t2} cos(θ)·cos(nθ) dθ"""
    if n == 1:
        # ∫ cos^2θ dθ = θ/2 + sin(2θ)/4
        return (0.5*t2 + 0.25*_S(2*t2)) - (0.5*t3 + 0.25*_S(2*t3))
    return 0.5 * (
        (_S((n-1)*t2) - _S((n-1)*t3)) / (n-1) +
        (_S((n+1)*t2) - _S((n+1)*t3)) / (n+1)
    )

def A2a_tri(n, t3, t2):
    """∫_{t3}^{t2} sin(θ)·cos(nθ) dθ"""
    if n == 1:
        # ∫ sinθ cosθ dθ = -cos(2θ)/4
        return (-_C(2*t2) + _C(2*t3)) / 4.0
    return 0.5 * (
        (_C((n-1)*t2) - _C((n-1)*t3)) / (n-1) -
        (_C((n+1)*t2) - _C((n+1)*t3)) / (n+1)
    )

def B1a_tri(n, t3, t2):
    """∫_{t3}^{t2} cos(θ)·sin(nθ) dθ"""
    if n == 1:
        # ∫ sinθ cosθ dθ = -cos(2θ)/4
        return (-_C(2*t2) + _C(2*t3)) / 4.0
    return 0.5 * (
        (-_C((n+1)*t2) + _C((n+1)*t3)) / (n+1) +
        (-_C((n-1)*t2) + _C((n-1)*t3)) / (n-1)
    )

def B2a_tri(n, t3, t2):
    """∫_{t3}^{t2} sin(θ)·sin(nθ) dθ"""
    if n == 1:
        # ∫ sin^2θ dθ = θ/2 − sin(2θ)/4
        return (0.5*t2 - 0.25*_S(2*t2)) - (0.5*t3 - 0.25*_S(2*t3))
    return 0.5 * (
        (_S((n-1)*t2) - _S((n-1)*t3)) / (n-1) -
        (_S((n+1)*t2) - _S((n+1)*t3)) / (n+1)
    )

# ================== TRECHO "b" — limites [θ3+π/3, θ2+π/3] ==================
# Mudança de variável φ = θ − π/3  → limites voltam a [t3, t2]
# e usa-se:
#  cos(n(φ+π/3)) = cos(nπ/3)cos(nφ) − sin(nπ/3)sin(nφ)
#  sin(n(φ+π/3)) = sin(nπ/3)cos(nφ) + cos(nπ/3)sin(nφ)

def A1b_tri(n, t3, t2):
    """∫_{t3+π/3}^{t2+π/3} cos(θ−π/3)·cos(nθ) dθ"""
    c = _C(n*math.pi/3.0)
    s = _S(n*math.pi/3.0)
    # combinação linear das "bases" em [t3, t2]
    return c*A1a_tri(n, t3, t2) - s*B1a_tri(n, t3, t2)

def A2b_tri(n, t3, t2):
    """∫_{t3+π/3}^{t2+π/3} sin(θ−π/3)·cos(nθ) dθ"""
    c = _C(n*math.pi/3.0)
    s = _S(n*math.pi/3.0)
    return c*A2a_tri(n, t3, t2) - s*B2a_tri(n, t3, t2)

def B1b_tri(n, t3, t2):
    """∫_{t3+π/3}^{t2+π/3} cos(θ−π/3)·sin(nθ) dθ"""
    c = _C(n*math.pi/3.0)
    s = _S(n*math.pi/3.0)
    return c*B1a_tri(n, t3, t2) + s*A1a_tri(n, t3, t2)

def B2b_tri(n, t3, t2):
    """∫_{t3+π/3}^{t2+π/3} sin(θ−π/3)·sin(nθ) dθ"""
    c = _C(n*math.pi/3.0)
    s = _S(n*math.pi/3.0)
    return c*B2a_tri(n, t3, t2) + s*A2a_tri(n, t3, t2)

# ================== SOMA DOS DOIS TRECHOS (definição do livro) ==================

def A1_tri(n, t3, t2): return A1a_tri(n, t3, t2) + A1b_tri(n, t3, t2)
def A2_tri(n, t3, t2): return A2a_tri(n, t3, t2) + A2b_tri(n, t3, t2)
def B1_tri(n, t3, t2): return B1a_tri(n, t3, t2) + B1b_tri(n, t3, t2)
def B2_tri(n, t3, t2): return B2a_tri(n, t3, t2) + B2b_tri(n, t3, t2)

# ================== COEFICIENTES NORMALIZADOS (trifásico) ==================

def coeffs_normalized_tri(n: int, x: float, w: float):
    """
    Retorna (ā_n, b̄_n, c̄_n) normalizados para a harmônica n (convenção do livro: sem 1/π).
    x = Vcmin/Vp   ;   w = ωCR  (Fig./Eqs. do livro)
    """
    theta3 = math.asin(x)             # (10.84)
    theta2 = math.pi + math.atan(-w)  # (10.83)
    if theta2 <= theta3:
        return float('nan'), float('nan'), float('nan')
    y = w                              # (10.168) y = ωCR

    a_bar = y*A1_tri(n, theta3, theta2) + A2_tri(n, theta3, theta2)   # (10.169)
    b_bar = y*B1_tri(n, theta3, theta2) + B2_tri(n, theta3, theta2)   # (10.170)
    c_bar = math.hypot(a_bar, b_bar)                                  # (10.179)
    return a_bar, b_bar, c_bar

def wrap_pi(x):
    while x <= -math.pi: x += 2*math.pi
    while x >  math.pi: x -= 2*math.pi
    return x

# Convenção: positivo = corrente ADIANTA a tensão de linha
def phi_lead_deg(a1, b1):
    phi = -math.atan2(b1, a1) + 2*math.pi/3   # ajuste de 120°
    return math.degrees(wrap_pi(phi))

def phi_book_deg(x, w):
    a1_bar, b1_bar, _ = coeffs_normalized_tri(1, x, w)
    return -math.degrees(math.atan2(b1_bar, a1_bar)) + 30.0

def sanity_checks():
    xs = [0.87, 0.90, 0.95, 0.98]
    w  = 0.6  # típico
    print("n múltiplos de 6 devem zerar:")
    for n in [6, 12, 18]:
        a, b, c = coeffs_normalized_tri(n, xs[0], w)
        print(n, a, b, c)
    print("ângulos (livro):")
    for x in xs:
        print(x, round(phi_book_deg(x, w), 2))

def compute_curves_tri(x_min=0.50, x_max=0.985, n_pts=150, nmax=50, harmonics_mode="odd"):
    """
    harmonics_mode:
      - "odd"  : usa somente ímpares (3,5,7,...)  [padrão, como no livro]
      - "all"  : usa pares e ímpares (2..nmax)
      - "even" : usa somente pares (2,4,6,...)
    """
    xs = np.linspace(x_min, x_max, n_pts)
    w_vals = np.array([solve_w_for_x_tri(float(x)) for x in xs])
    y_vals = np.array([Ic_eff_over_Vpk_R_tri(float(w), float(x)) for x, w in zip(xs, w_vals)])

    # seleção de harmônicos conforme modo
    if harmonics_mode == "all":
        harmonics = range(2, nmax + 1)
    elif harmonics_mode == "even":
        harmonics = range(2, nmax + 1, 2)
    else:  # "odd"
        harmonics = range(3, nmax + 1, 2)

    c1_list, thd_list, phi1_deg_list, pf_list = [], [], [], []
    for x, w in zip(xs, w_vals):
        a1, b1, c1 = coeffs_normalized_tri(1, x, w)

        # soma das potências harmônicas (n conforme seleção acima)
        c_sq = 0.0
        for n in harmonics:
            _, _, cn = coeffs_normalized_tri(n, x, w)
            c_sq += cn * cn

        thd = math.sqrt(c_sq) / c1
        phi1_deg = phi_lead_deg(a1, b1)
        pf = math.cos(math.radians(phi1_deg)) / math.sqrt(1.0 + thd*thd)

        c1_list.append(c1)
        thd_list.append(thd)
        phi1_deg_list.append(phi1_deg)
        pf_list.append(pf)

    df_tri = pd.DataFrame({
        "Vcmin_over_Vpk": xs,
        "omegaRC": w_vals,
        "R_Ic_eff_over_Vpk": y_vals,
        "c1_bar": c1_list,
        "THD_pu": thd_list,
        "phi1_deg": phi1_deg_list,
        "PF": pf_list
    })
    return df_tri


def make_plots_and_save_tri(df_tri, save_prefix: str = "fig_tri"):
    import matplotlib.pyplot as plt
    # Fig. 10.22 — ωRC vs Vcmin/Vp (trifásico)
    plt.figure()
    plt.plot(df_tri["Vcmin_over_Vpk"], df_tri["omegaRC"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{p}$"); plt.ylabel(r"$\omega RC$")
    plt.title("Fig. 10.22 — ωRC vs Vcmin/Vp (trifásico)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_22_tri.png", dpi=200, bbox_inches="tight")

    # Fig. 10.24 — (R·Ic_ef/Vp) vs ωRC (trifásico)
    plt.figure()
    plt.plot(df_tri["omegaRC"], df_tri["R_Ic_eff_over_Vpk"], linewidth=2)
    plt.xlabel(r"$\omega RC$"); plt.ylabel(r"$R\cdot I_{c,ef}/V_{p}$")
    plt.title("Fig. 10.24 — (R·Ic_ef/Vp) vs ωRC (trifásico)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_24_tri.png", dpi=200, bbox_inches="tight")

    # Fig. 10.27
    plt.figure()
    plt.plot(df_tri["Vcmin_over_Vpk"], df_tri["THD_pu"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel("THD (pu)")
    plt.title("Fig. 10.27 — THD (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_27_THD.png", dpi=200, bbox_inches="tight")

    # Fig. 10.28
    plt.figure()
    plt.plot(df_tri["Vcmin_over_Vpk"], df_tri["phi1_deg"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel(r"$\phi_1$ (graus)")
    plt.title("Fig. 10.28 — Ângulo de deslocamento (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_28_phi1.png", dpi=200, bbox_inches="tight")

    # Fig. 10.29
    plt.figure()
    plt.plot(df_tri["Vcmin_over_Vpk"], df_tri["PF"], linewidth=2)
    plt.xlabel(r"$V_{cmin}/V_{pk}$"); plt.ylabel("FP")
    plt.title("Fig. 10.29 — Fator de potência (reprodução)")
    plt.grid(True); plt.savefig(f"{save_prefix}_10_29_PF.png", dpi=200, bbox_inches="tight")


# ---- CLI hook: adicionar flags opcionais para também gerar o trifásico ----
try:
    _prev_parse_args = parse_args  # preserva o original se existir
except NameError:
    _prev_parse_args = None

# if __name__ == "__main__":
#     try:
#         args  # se o main original definiu 'args'
#     except NameError:
#         pass
#     else:
#         if getattr(args, "tri_csv", None):
#             df_tri = compute_curves_tri(args.xmin, args.xmax, args.pts)
#             # salva CSV e figuras
#             import pandas as pd
#             df_tri.to_csv(args.tri_csv, index=False)
#             make_plots_and_save_tri(df_tri, args.tri_prefix)

if __name__ == "__main__":
    args = parse_args()

    # --- bloco monofásico ---
    df_mono = compute_curves(
        x_min=args.xmin,
        x_max=args.xmax,
        n_pts=args.pts,
        nmax=args.nmax,
        harmonics_mode=args.harmonics
    )
    df_mono.to_csv(args.csv, index=False)
    make_plots_and_save(df_mono, args.prefix)

    # --- bloco trifásico (só roda se você passar os parâmetros) ---
    df_tri = compute_curves_tri(args.xmin, args.xmax, args.pts, nmax=args.nmax, harmonics_mode=args.harmonics)
    df_tri.to_csv(args.csv, index=False)
    make_plots_and_save_tri(df_tri, args.prefix)

    plt.show()

# ======================  FIM DO BLOCO NOVO — TRIFÁSICO  ======================
