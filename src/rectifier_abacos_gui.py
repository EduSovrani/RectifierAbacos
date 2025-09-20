#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Retificador + Filtro C — Projeto via Ábacos (Cap. 10)
"""

import os
import tkinter as tk
from tkinter import ttk
import sys
import importlib # Lazy imports (adiam NumPy e Matplotlib/TkAgg)

import sys, os, ctypes
from pathlib import Path

def _resource_path(relpath: str) -> str:
    """Resolve caminho no bundle (PyInstaller) ou em dev."""
    base = getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)
    return str(Path(base) / relpath)

def _read_app_version() -> str:
    """
    Ordem:
    1) VERSION (bundle ou projeto)
    2) Recurso de versão do EXE (apenas quando frozen)
    3) 'DEV'
    """
    # — 1) procurar VERSION em vários lugares —
    here = Path(__file__).resolve().parent
    candidates = [
        Path(_resource_path("VERSION")),                # dentro do bundle
        Path(_resource_path(os.path.join("assets","VERSION"))),
        here / "VERSION",                               # mesmo dir do .py
        here.parent / "VERSION",                        # raiz do projeto (subindo 1)
        here.parent.parent / "VERSION",                 # (subindo 2, por garantia)
    ]
    for p in candidates:
        try:
            if p.is_file():
                v = p.read_text(encoding="utf-8").strip()
                if v:
                    return v
        except Exception:
            pass

    # — 2) se estiver empacotado, ler “FileVersion” do próprio EXE —
    if getattr(sys, "frozen", False) and os.name == "nt":
        try:
            GetFileVersionInfoSizeW = ctypes.windll.version.GetFileVersionInfoSizeW
            GetFileVersionInfoW     = ctypes.windll.version.GetFileVersionInfoW
            VerQueryValueW          = ctypes.windll.version.VerQueryValueW

            exe = sys.executable
            size = GetFileVersionInfoSizeW(exe, None)
            if size:
                data = ctypes.create_string_buffer(size)
                if GetFileVersionInfoW(exe, 0, size, data):
                    lp = ctypes.c_void_p()
                    ln = ctypes.c_uint()
                    if VerQueryValueW(data, u"\\VarFileInfo\\Translation", ctypes.byref(lp), ctypes.byref(ln)):
                        arr = (ctypes.c_ushort * 2).from_address(lp.value)
                        lang, codepage = arr[0], arr[1]
                        sub = u"\\StringFileInfo\\%04x%04x\\FileVersion" % (lang, codepage)
                        if VerQueryValueW(data, sub, ctypes.byref(lp), ctypes.byref(ln)):
                            return ctypes.wstring_at(lp.value, ln.value).strip()
        except Exception:
            pass

    # — 3) fallback —
    return "DEV"

class _LazyMod:
    """Carrega o módulo na 1ª vez que um atributo é acessado."""
    def __init__(self, modname):
        self._modname = modname
        self._loaded = None
    def _load(self):
        if self._loaded is None:
            self._loaded = importlib.import_module(self._modname)
        return self._loaded
    def __getattr__(self, name):
        return getattr(self._load(), name)

# Atrasamos NumPy até o 1º uso
np = _LazyMod("numpy")

# Helpers para carregar Matplotlib apenas quando necessário e já com TkAgg
def _ensure_mpl_with_tk():
    """
    Garante backend TkAgg e retorna:
      (plt, FigureCanvasTkAgg, NavigationToolbar2Tk, mpl)
    """
    # Seleciona backend ANTES de importar pyplot
    mpl = importlib.import_module("matplotlib")
    try:
        # Só troca backend se ainda não definido explicitamente
        if getattr(mpl, "get_backend", None) and mpl.get_backend().lower() != "tkagg":
            mpl.use("TkAgg")
    except Exception:
        # se falhar, segue — TkAgg pode já estar ativo
        pass

    plt = importlib.import_module("matplotlib.pyplot")
    tkagg = importlib.import_module("matplotlib.backends.backend_tkagg")
    FigureCanvasTkAgg = getattr(tkagg, "FigureCanvasTkAgg")
    NavigationToolbar2Tk = getattr(tkagg, "NavigationToolbar2Tk")

    # Aplica RCParams (equivalente ao bloco antigo)
    rc = {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "mathtext.fontset": "dejavusans",
        "axes.titlepad": 8.0,
        "axes.labelpad": 6.0,
    }
    try:
        mpl.rcParams.update(rc)
    except Exception:
        pass

    return plt, FigureCanvasTkAgg, NavigationToolbar2Tk, mpl


# Strings LaTeX-like (MathText) — usar SEMPRE os mesmos em toda atualização
LBL_t       = r"$t\,[\mathrm{s}]$"
LBL_V       = r"$v\,[\mathrm{V}]$"
LBL_I       = r"$i\,[\mathrm{A}]$"
LBL_wRC     = r"$\omega RC$"
LBL_x       = r"$V_{L,\min}/V_{pk}$"
LBL_RIcVpk  = r"$R\cdot I_{c,\mathrm{rms}}/V_{pk}$"
LBL_THD     = r"$\mathrm{THD}$"
LBL_phi1    = r"$\varphi_1\;[\mathrm{graus}]$"
LBL_PF      = r"$\mathrm{FP}$"

TTL_top     = r"Tensões: $v_{\mathrm{out,ret}}(t)$ e $v_{\mathrm{out}}(t)$"
TTL_mid     = r"Corrente de carga: $i_{\mathrm{load}}(t)$"
TTL_bot     = r"Corrente no diodo #2: $i_{D2}(t)$"
TTL_bot_na  = r"$i_{D2}(t)$ — N/A (6 pulsos)"

# ========== Parâmetros gerais ==========
FS = 1e6                      # Fs de simulação
DT = 1/FS                     # dt de simulação
CYCLES_SIM = 6                # total de ciclos simulados
WIN_START = 2.0               # início da janela (em ciclos)
WIN_END   = 6.0               # fim da janela (em ciclos)
VHYST = 0.5                   # histerese em volts

# ---------- ÁBACOS MONO (fácil de ajustar) ----------
AB_XMIN = 0.30
AB_XMAX = 0.985
AB_PTS  = 150
AB_NMAX = 50
AB_HARM = "odd"   # "odd", "all", "even"

# ---------- ÁBACOS TRI (fácil de ajustar) ----------
AB_XMIN_TRI = 0.87
AB_XMAX_TRI = 0.98
AB_PTS_TRI  = 100
AB_NMAX_TRI = 50
AB_HARM_TRI = "odd"   # "odd", "all", "even"

# ========== NOVOS LIMITES DOS ÁBACOS (por modo) ==========
# Convenção: cada entrada tem "x":(xmin,xmax) e "y":(ymin,ymax)
AB_LIMITS = {
    "MONO": {
        "wRC_vs_VLmin":            {"x": (0.30, 1.00), "y": (0.0, 200.0)},
        "R_Icrms_over_Vpk_vs_wRC": {"x": (0.0, 200.0), "y": (0.0, 5.0)},
        "THD_vs_RIc":              {"x": (0.3, 1.0),   "y": (0.0, 4.0)},
        "phi1_vs_RIc":             {"x": (0.3, 1.0),   "y": (0.0, 35.0)},
        "FP_vs_RIc":               {"x": (0.3, 1.0),   "y": (0.20, 0.80)},
    },
    "TRI": {
        "wRC_vs_VLmin":            {"x": (0.86, 1.00), "y": (0.0, 45.0)},
        "R_Icrms_over_Vpk_vs_wRC": {"x": (0.0, 45.0),  "y": (0.5, 2.5)},
        "THD_vs_RIc":              {"x": (0.86, 1.0),  "y": (0.5, 2.25)},
        "phi1_vs_RIc":             {"x": (0.86, 1.0),  "y": (7.0, 13.5)},
        "FP_vs_RIc":               {"x": (0.86, 1.0),  "y": (0.40, 0.85)},
    },
}

# ========== Utilidades ==========
def trapezoid(y, x):
    try:
        return np.trapezoid(y, x)
    except AttributeError:
        return np.trapz(y, x)

def slice_window(t, arr, f):
    """Recorta a janela comum a todos os sinais."""
    T = 1.0/f
    t0 = WIN_START*T
    t1 = WIN_END*T
    m = (t >= t0) & (t < t1)
    return t[m], arr[m]

# ========== Import do módulo de ábacos (rectifier_abacos_core.py) ==========
import importlib.util
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CAPS_PATH = os.path.join(THIS_DIR, "rectifier_abacos_core.py")
spec = importlib.util.spec_from_file_location("caps", CAPS_PATH)
caps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(caps)

# Trifásico (novo)
CAPS_TRI_PATH = os.path.join(THIS_DIR, "rectifier_abacos_core.py")
spec_tri = importlib.util.spec_from_file_location("caps_tri", CAPS_TRI_PATH)
caps_tri = importlib.util.module_from_spec(spec_tri)
spec_tri.loader.exec_module(caps_tri)

# Utilitário: carrega o conjunto de ábacos conforme a topologia e padroniza colunas
def _load_abacos(topology: str):
    if topology.strip().lower().startswith("trif"):
        try:
            df = caps.compute_curves_tri(x_min=AB_XMIN_TRI, x_max=AB_XMAX_TRI, n_pts=AB_PTS_TRI, nmax=AB_NMAX_TRI, harmonics_mode=AB_HARM_TRI)
        except TypeError:
            df = caps.compute_curves_tri()
        return df
    else:
        try:
            df = caps.compute_curves(x_min=AB_XMIN, x_max=AB_XMAX, n_pts=AB_PTS, nmax=AB_NMAX, harmonics_mode=AB_HARM)
        except TypeError:
            df = caps.compute_curves()
        return df

# ========== Simuladores ==========
def simulate_single_phase(f, Vph_rms, C, I_target):
    fs = FS
    w = 2*np.pi*f
    T = 1.0/f
    dt = 1.0/fs
    t = np.arange(0, CYCLES_SIM*T, dt)

    Vm1 = np.sqrt(2)*Vph_rms
    s = np.sin(w*t)
    v_rect = np.abs(Vm1*s)  # retificação de onda completa
    R = Vm1 / max(I_target, 1e-9)

    v_c   = np.zeros_like(t)
    cond  = np.zeros_like(t, dtype=bool)
    is_on = False
    for k in range(1, len(t)):
        vcp   = v_c[k-1]
        vrk   = v_rect[k]
        dvrec = v_rect[k] - v_rect[k-1]   # ~ dv_rect/dt (sem dividir por dt)
        # epsilon baseado no decaimento numérico esperado em 1 passo
        eps_v = max(DT, 1.0 * (vcp/R)/C * dt)

        if not is_on:
            # liga se passar acima da banda OU estiver praticamente empatado e subindo
            if (vrk >= vcp + eps_v) or (abs(vrk - vcp) <= eps_v and dvrec >= 0):
                is_on = True
        else:
            # desliga só se cair abaixo da banda E estiver descendo
            if (vrk < vcp - eps_v) and (dvrec <= 0):
                is_on = False

        cond[k] = is_on
        # física inalterada:
        v_c[k] = vrk if is_on else (vcp - (vcp/R)/C * dt)

    i_load = v_c / R
    # derivada central suave do v_rect (melhor que diff)
    dv = np.empty_like(v_c)
    dv[1:-1] = (v_c[2:] - v_c[:-2]) / (2*dt)
    dv[0]    = (v_c[1] - v_c[0]) / dt
    dv[-1]   = (v_c[-1] - v_c[-2]) / dt
    i_c = C * dv
    i_rect = np.where(cond, i_c + i_load, 0.0)
    i_d2 = np.where((s >= 0.0) & cond, i_rect, 0.0)  # diodo #2 na semiciclo positivo

    # Janela consistente
    t_win, vsrc_win = slice_window(t, v_rect, f)
    _, vout_win = slice_window(t, v_c, f)
    _, iout_win = slice_window(t, i_load, f)
    _, id2_win  = slice_window(t, i_d2, f)

    duration = max(t_win[-1] - t_win[0], 1e-9)
    v_avg = trapezoid(vout_win, t_win)/duration
    v_rms = np.sqrt(trapezoid(vout_win**2, t_win)/duration)
    ripple_pp = float(vout_win.max() - vout_win.min())
    N = int(T*fs)
    duty_last = float(np.mean(cond[-N:])) if N > 0 else float(np.mean(cond))

    id2_avg = trapezoid(id2_win, t_win)/duration
    id2_rms = np.sqrt(trapezoid(id2_win**2, t_win)/duration)

    Vpk = Vm1
    x_sim = (vout_win.min()/Vpk) if Vpk > 0 else 0.0
    wrc_sim = (w * (Vpk/max(I_target,1e-12)) * C)

    metrics = dict(
        V_out_avg=v_avg, V_out_rms=v_rms, Ripple_pp=ripple_pp, Duty=duty_last,
        Id2_avg=id2_avg, Id2_rms=id2_rms, Vpk=Vpk, x_sim=x_sim, wrc_sim=wrc_sim
    )
    return (t_win, vsrc_win, vout_win, iout_win, id2_win, R, metrics)

def simulate_three_phase(f, Vph_rms, C, I_target):
    fs = FS
    w = 2*np.pi*f
    T = 1.0/f
    dt = 1.0/fs
    t = np.arange(0, CYCLES_SIM*T, dt)

    Vm = np.sqrt(2)*Vph_rms
    va = Vm*np.sin(w*t)
    vb = Vm*np.sin(w*t - 2*np.pi/3)
    vc = Vm*np.sin(w*t + 2*np.pi/3)
    # Aproxima o 6-pulsos: diferença entre maior e menor fase
    v_rect = np.maximum.reduce([va, vb, vc]) - np.minimum.reduce([va, vb, vc])

    Vll_peak = np.sqrt(2)*np.sqrt(3)*Vph_rms
    R = Vll_peak / max(I_target, 1e-9)

    v_c   = np.zeros_like(t)
    cond  = np.zeros_like(t, dtype=bool)
    is_on = False
    for k in range(1, len(t)):
        vcp   = v_c[k-1]
        vrk   = v_rect[k]
        dvrec = v_rect[k] - v_rect[k-1]   # ~ dv_rect/dt (sem dividir por dt)
        # epsilon baseado no decaimento numérico esperado em 1 passo
        eps_v = max(DT, 1.0 * (vcp/R)/C * dt)

        if not is_on:
            # liga se passar acima da banda OU estiver praticamente empatado e subindo
            if (vrk >= vcp + eps_v) or (abs(vrk - vcp) <= eps_v and dvrec >= 0):
                is_on = True
        else:
            # desliga só se cair abaixo da banda E estiver descendo
            if (vrk < vcp - eps_v) and (dvrec <= 0):
                is_on = False

        cond[k] = is_on
        # física inalterada:
        v_c[k] = vrk if is_on else (vcp - (vcp/R)/C * dt)

    i_load = v_c / R
    # derivada central suave do v_rect (melhor que diff)
    dv = np.empty_like(v_c)
    dv[1:-1] = (v_c[2:] - v_c[:-2]) / (2*dt)
    dv[0]    = (v_c[1] - v_c[0]) / dt 
    dv[-1]   = (v_c[-1] - v_c[-2]) / dt
    i_c = C * dv
    # Corrente no retificador só quando conduz (carregando C + alimentando carga)
    i_rect = np.where(cond, i_c + i_load, 0.0)

    # Definir Id2 como corrente no diodo positivo da fase A
    # Conduz quando a fase A é a maior entre (va,vb,vc) e o retificador está em condução
    maxA = (va >= vb) & (va >= vc)
    d2_on = cond & maxA
    i_d2 = np.where(d2_on, i_rect, 0.0)

    # Janela consistente
    t_win, vsrc_win = slice_window(t, v_rect, f)
    _, vout_win = slice_window(t, v_c, f)
    _, iout_win = slice_window(t, i_load, f)
    _, id2_win  = slice_window(t, i_d2, f)

    duration = max(t_win[-1] - t_win[0], 1e-9)
    v_avg = trapezoid(vout_win, t_win)/duration
    v_rms = np.sqrt(trapezoid(vout_win**2, t_win)/duration)
    ripple_pp = float(vout_win.max() - vout_win.min())
    N = int(T*fs)
    duty_last = float(np.mean(cond[-N:])) if N > 0 else float(np.mean(cond))

    id2_avg = trapezoid(id2_win, t_win)/duration
    id2_rms = np.sqrt(trapezoid(id2_win**2, t_win)/duration)

    metrics = dict(V_out_avg=v_avg, V_out_rms=v_rms, Ripple_pp=ripple_pp, Duty=duty_last,
                   Id2_avg=id2_avg, Id2_rms=id2_rms)
    return (t_win, vsrc_win, vout_win, iout_win, id2_win, R, metrics)

# ========== Filtrinho para curvas dos ábacos ==========
def _clean_xy(x, y, x_min=None, x_max=None, yabs=1e4):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if x_min is not None: m &= (x >= x_min)
    if x_max is not None: m &= (x <= x_max)
    x = x[m]; y = y[m]
    if x.size > 1:
        idx = np.argsort(x, kind="mergesort")
        x = x[idx]; y = y[idx]
    if y.size > 3:
        dy = np.diff(y)
        q = np.nanpercentile(np.abs(dy), 97.5)
        bad = np.where(np.abs(dy) > max(q*4, 1e-9))[0]
        if bad.size:
            cut = bad[0]+1
            x = x[:cut]; y = y[:cut]
    y = np.clip(y, -yabs, yabs)
    return x, y

# ========== App Tk ==========
class App(tk.Tk):
    def __init__(self):
        super().__init__()

        # ---- título com versão ----
        self.version = _read_app_version()
        self.title(f"Retificador + Filtro C — Ábacos (v{self.version})")

        # ---- ícone da janela (procura em dev e no bundle one-file) ----
        import sys
        here = os.path.dirname(os.path.abspath(__file__))   # ...\projeto\src
        root = os.path.dirname(here)                        # ...\projeto

        candidates = []
        # quando congelado (PyInstaller one-file): arquivos extraídos em _MEIPASS
        if getattr(sys, "frozen", False):
            base = getattr(sys, "_MEIPASS", here)
            candidates += [
                os.path.join(base, "assets", "icon.ico"),
                os.path.join(base, "icon.ico"),
            ]
        # em desenvolvimento (rodando do source) e fallbacks gerais
        candidates += [
            os.path.join(root, "assets", "icon.ico"),  # onde seu ícone realmente está
            os.path.join(here, "assets", "icon.ico"),
            os.path.join(root, "icon.ico"),
            os.path.join(here, "icon.ico"),
        ]

        for p in candidates:
            if os.path.exists(p):
                try:
                    # Windows: .ico direto
                    self.iconbitmap(default=p)
                    break
                except Exception:
                    # último recurso: tentar como PhotoImage (se o Tk aceitar .ico)
                    try:
                        _img = tk.PhotoImage(file=p)
                        self.wm_iconphoto(True, _img)
                        break
                    except Exception:
                        pass

        # ---- dimensões da janela ----
        self.geometry("1200x800")

        # ---- carrega Matplotlib/Tk na 1ª necessidade ----
        (self._plt,
        self._FigureCanvasTkAgg,
        self._NavigationToolbar2Tk,
        self._mpl) = _ensure_mpl_with_tk()

        # ===== Banner de Avisos (UX sem pop-up) =====
        self._warn_var = tk.StringVar(value="")
        self._warnings = []          # acumula mensagens desta execução
        self._banner = ttk.Label(self, textvariable=self._warn_var,
                                background="#FFF3CD", foreground="#8A6D3B",
                                anchor="w", padding=(8,4))
        # escondido inicialmente
        self._banner_visible = False

        # Ponto persistente nos ábacos
        self._ab_point = (None, None)

        # Notebook
        nb = ttk.Notebook(self)
        # o banner vai no topo; notebook abaixo
        self._build_banner()
        nb.pack(fill="both", expand=True)
        self.nb = nb

        # Tabs
        self.tab_proj = ttk.Frame(nb, padding=8)
        self.tab_abacos = ttk.Frame(nb, padding=8)
        nb.add(self.tab_proj, text="Projeto + Simulação")
        nb.add(self.tab_abacos, text="Ábacos")

        self._build_tab_projeto(self.tab_proj)
        self._build_tab_abacos(self.tab_abacos)

    # ===== Banner helpers =====
    def _build_banner(self):
        # não mostra nada até ter aviso
        pass

    def _show_banner(self, text: str):
        self._warn_var.set(text)
        if not self._banner_visible:
            self._banner.pack(side="top", fill="x")
            self._banner_visible = True

    def _hide_banner(self):
        if self._banner_visible:
            self._banner.pack_forget()
            self._banner_visible = False
        self._warn_var.set("")

    def _add_warning(self, msg: str):
        if msg not in self._warnings:
            self._warnings.append(msg)

    def _flush_warnings(self):
        if self._warnings:
            self._show_banner(" | ".join(self._warnings))
        else:
            self._hide_banner()

    # ====== Helpers de limites e clamp ======
    def _ab_mode(self) -> str:
        topo = (self.topology.get() or "").lower()
        return "TRI" if "tri" in topo else "MONO"

    def _set_lims(self, ax, key: str):
        mode = self._ab_mode()
        lims = AB_LIMITS[mode][key]
        ax.set_xlim(*lims["x"])
        ax.set_ylim(*lims["y"])
        ax.grid(True, which="both", linestyle=":", linewidth=0.7)

    def _r_at_w(self, wrc_value, df):
        # interpola R*Ic_rms/Vpk em função de ωRC
        xs = np.asarray(df["omegaRC"], float)
        ys = np.asarray(df["R_Ic_eff_over_Vpk"], float)
        if wrc_value <= xs.min(): return float(ys[xs.argmin()])
        if wrc_value >= xs.max(): return float(ys[xs.argmax()])
        return float(np.interp(wrc_value, xs, ys))

    def _clamp_to_range(self, val, rng):
        lo, hi = float(rng[0]), float(rng[1])
        return max(lo, min(hi, float(val)))

    def _clamp_point_to_axes(self, x, w):
        """
        Ajusta o ponto (x, ωRC) aos limites do 1º gráfico dos ábacos.
        Se houver clamp, registra um aviso na faixa amarela.
        """
        mode = self._ab_mode()
        x_rng = AB_LIMITS[mode]["wRC_vs_VLmin"]["x"]
        y_rng = AB_LIMITS[mode]["wRC_vs_VLmin"]["y"]
        cx = self._clamp_to_range(x, x_rng)
        cw = self._clamp_to_range(w, y_rng)
        if (cx != x) or (cw != w):
            self._add_warning("Ponto (x, ωRC) ajustado aos limites do gráfico.")
        return cx, cw

    # ---------- GUI: Projeto ----------
    def _build_tab_projeto(self, parent):
        top = ttk.Frame(parent); top.pack(fill="both", expand=True)
        left = ttk.Frame(top); left.pack(side="left", fill="y", padx=(0,8))
        right = ttk.Frame(top); right.pack(side="left", fill="both", expand=True)

        # Entradas
        self.topology = tk.StringVar(value="Monofásico (ponte)")
        self.f = tk.DoubleVar(value=60.0)
        self.Vph = tk.DoubleVar(value=127.0)
        self.I = tk.DoubleVar(value=10.0)
        self.ripple = tk.DoubleVar(value=31.0)

        row = 0
        ttk.Label(left, text="Parâmetros", font=("TkDefaultFont", 10, "bold")).grid(row=row, column=0, columnspan=2, sticky="w", pady=(0,6)); row+=1
        ttk.Label(left, text="Topologia:").grid(row=row, column=0, sticky="e")
        ttk.Combobox(left, textvariable=self.topology, values=["Monofásico (ponte)", "Trifásico 6 pulsos"], state="readonly", width=20).grid(row=row, column=1, sticky="w"); row+=1
        ttk.Label(left, text="f [Hz]:").grid(row=row, column=0, sticky="e")
        ttk.Entry(left, textvariable=self.f, width=12).grid(row=row, column=1, sticky="w"); row+=1
        ttk.Label(left, text="V_fase_rms [V]:").grid(row=row, column=0, sticky="e")
        ttk.Entry(left, textvariable=self.Vph, width=12).grid(row=row, column=1, sticky="w"); row+=1
        ttk.Label(left, text="I alvo [A]:").grid(row=row, column=0, sticky="e")
        ttk.Entry(left, textvariable=self.I, width=12).grid(row=row, column=1, sticky="w"); row+=1
        ttk.Label(left, text="ΔV alvo [Vpp]:").grid(row=row, column=0, sticky="e")
        ttk.Entry(left, textvariable=self.ripple, width=12).grid(row=row, column=1, sticky="w"); row+=1

        ttk.Button(left, text="Calcular via Ábaco + Simular", command=self.run_project).grid(row=row, column=0, columnspan=2, pady=8); row+=1

        lf = ttk.LabelFrame(left, text="Resultados (janela de 4 ciclos)")
        lf.grid(row=row, column=0, columnspan=2, sticky="nsew", pady=(6,0))
        row+=1

        
        # Labels idênticos aos finais (com placeholders "—")
        self.res_vars = {
            "Topologia": tk.StringVar(value="Topologia: —"),
            "Entrada": tk.StringVar(value="Entrada: —"),
            "x = $V_{L,\\min}/V_{pk}$": tk.StringVar(value="x = V_(L_min)/V_pk = -"),
            "ωRC (ábaco)": tk.StringVar(value="ωRC = —"),
            "C (dimensionado)": tk.StringVar(value="C = —"),
            "R (adotado)": tk.StringVar(value="R ≈ — Ω"),
            "Ic_rms (ábaco)": tk.StringVar(value="Ic_rms (ábaco) ≈ — A"),
            "THD (ábaco)": tk.StringVar(value="THD (ábaco) ≈ — pu"),
            "φ1 (ábaco)": tk.StringVar(value="φ1 (ábaco) ≈ — °"),
            "FP (ábaco)": tk.StringVar(value="FP (ábaco) ≈ —"),
            "Vout_médio": tk.StringVar(value="V̄_out = — V"),
            "Vout_RMS": tk.StringVar(value="V_out,rms = — V"),
            "Ripple_sim": tk.StringVar(value="ΔV_pp = — V"),
            "I_diodo2_médio": tk.StringVar(value="I_D2, méd = — A"),
            "I_diodo2_RMS": tk.StringVar(value="I_D2, rms = — A"),
        }
        rlf = 0
        for k in self.res_vars:
            ttk.Label(lf, textvariable=self.res_vars[k]).grid(row=rlf, column=0, sticky="w")
            rlf += 1

        # Plots (lazy com _ensure_mpl_with_tk já chamado no __init__)
        plot_frame = ttk.Frame(right); plot_frame.pack(fill="both", expand=True, padx=6, pady=6)

        # Top
        self.fig_top = self._plt.Figure(figsize=(7.8, 2.6), dpi=100, constrained_layout=True)
        self.ax_top = self.fig_top.add_subplot(111); self.ax_top.grid(True)
        self.ax_top.set_title(TTL_top)
        self.ax_top.set_xlabel(LBL_t); self.ax_top.set_ylabel(LBL_V)
        self.tb_top_frame = ttk.Frame(plot_frame); self.tb_top_frame.pack(side="top", fill="x")
        self.canvas_top = self._FigureCanvasTkAgg(self.fig_top, master=plot_frame)
        self.toolbar_top = self._NavigationToolbar2Tk(self.canvas_top, self.tb_top_frame)
        self.toolbar_top.update()
        self.canvas_top.get_tk_widget().pack(fill="both", expand=True)

        # Linha de baixo (meio + inferior)
        plot_row = ttk.Frame(right); plot_row.pack(fill="both", expand=True, padx=6, pady=(0,6))

        # Meio
        self.fig_mid = self._plt.Figure(figsize=(3.8, 2.5), dpi=100, constrained_layout=True)
        self.ax_mid = self.fig_mid.add_subplot(111); self.ax_mid.grid(True)
        self.ax_mid.set_title(TTL_mid); self.ax_mid.set_xlabel(LBL_t); self.ax_mid.set_ylabel(LBL_I)
        self.tb_mid_col = ttk.Frame(plot_row); self.tb_mid_col.pack(side="left", fill="both", expand=True, padx=(0,6))
        self.tb_mid_frame = ttk.Frame(self.tb_mid_col); self.tb_mid_frame.pack(side="top", fill="x")
        self.canvas_mid = self._FigureCanvasTkAgg(self.fig_mid, master=self.tb_mid_col)
        self.toolbar_mid = self._NavigationToolbar2Tk(self.canvas_mid, self.tb_mid_frame)
        self.toolbar_mid.update()
        self.canvas_mid.get_tk_widget().pack(side="top", fill="both", expand=True)

        # Inferior
        self.fig_bot = self._plt.Figure(figsize=(3.8, 2.5), dpi=100, constrained_layout=True)
        self.ax_bot = self.fig_bot.add_subplot(111); self.ax_bot.grid(True)
        self.ax_bot.set_title(TTL_bot); self.ax_bot.set_xlabel(LBL_t); self.ax_bot.set_ylabel(LBL_I)
        self.tb_bot_col = ttk.Frame(plot_row); self.tb_bot_col.pack(side="left", fill="both", expand=True)
        self.tb_bot_frame = ttk.Frame(self.tb_bot_col); self.tb_bot_frame.pack(side="top", fill="x")
        self.canvas_bot = self._FigureCanvasTkAgg(self.fig_bot, master=self.tb_bot_col)
        self.toolbar_bot = self._NavigationToolbar2Tk(self.canvas_bot, self.tb_bot_frame)
        self.toolbar_bot.update()
        self.canvas_bot.get_tk_widget().pack(side="top", fill="both", expand=True)


    # ---------- GUI: Ábacos ----------
    def _build_tab_abacos(self, parent):
        frame = ttk.Frame(parent); frame.pack(fill="both", expand=True, padx=6, pady=6)

        # Figura principal dos ábacos (lazy via _ensure_mpl_with_tk)
        self.fig_ab = self._plt.Figure(figsize=(11.0, 6.8), dpi=100, constrained_layout=True)
        axs = self.fig_ab.subplots(2, 3)
        self.ax_ab = [axs[0,0], axs[0,1], axs[0,2], axs[1,0], axs[1,1]]
        axs[1,2].set_visible(False)

        titles = [
            r"Fig. 10.9 — $\omega RC$ × $V_{L,\min}/V_{pk}$",
            r"Fig. 10.10 — $(R\cdot I_{c,\mathrm{rms}}/V_{pk})$ × $\omega RC$",
            r"Fig. 10.27 — THD vs $V_{L,\min}/V_{pk}$",
            r"Fig. 10.28 — $\varphi_1$ (graus) vs $V_{L,\min}/V_{pk}$",
            r"Fig. 10.29 — FP vs $V_{L,\min}/V_{pk}$"
        ]
        for ax, title in zip(self.ax_ab, titles):
            ax.grid(True)
            ax.set_title(title, fontsize=10)

        self.canvas_ab = self._FigureCanvasTkAgg(self.fig_ab, master=frame)
        self.tb_ab_frame = ttk.Frame(frame); self.tb_ab_frame.pack(side="top", fill="x")
        self.toolbar_ab = self._NavigationToolbar2Tk(self.canvas_ab, self.tb_ab_frame)
        self.toolbar_ab.update()
        self.canvas_ab.get_tk_widget().pack(fill="both", expand=True)


    # ---------- Desenho dos ábacos ----------
    def _interp_safe(self, x, xs, ys):
        x = float(x)
        xs = np.asarray(xs); ys = np.asarray(ys)
        if x <= xs.min(): return float(ys[xs.argmin()])
        if x >= xs.max(): return float(ys[xs.argmax()])
        return float(np.interp(x, xs, ys))

    def _interp_clip(self, val, xs, ys, *, label="valor"):
        """
        Interpola com saturação nos limites do ábaco e registra aviso.
        """
        xs = np.asarray(xs, float); ys = np.asarray(ys, float)
        v  = float(val)
        xmin, xmax = float(xs.min()), float(xs.max())
        if v <= xmin:
            self._add_warning(f"{label} abaixo do limite do ábaco; saturando em {xmin:.4g}.")
            return float(ys[xs.argmin()])
        if v >= xmax:
            self._add_warning(f"{label} acima do limite do ábaco; saturando em {xmax:.4g}.")
            return float(ys[xs.argmax()])
        return float(np.interp(v, xs, ys))

    def _clear(self, ax):
        ax.cla(); ax.grid(True)

    def redraw_abacos(self, x_point=None, w_point=None):
        df = getattr(self, "df_ab", None)
        if df is None or not len(df):
            return

        # ωRC vs x
        ax = self.ax_ab[0]; self._clear(ax)
        ax.set_xlabel(LBL_x); ax.set_ylabel(LBL_wRC)
        X, Y = _clean_xy(df["Vcmin_over_Vpk"], df["omegaRC"], None, None, yabs=5e3)
        ax.plot(X, Y, lw=2)
        self._set_lims(ax, "wRC_vs_VLmin")
        if (x_point is not None) and (w_point is not None):
            cx, cw = self._clamp_point_to_axes(x_point, w_point)
            ax.axvline(cx, ls="--"); ax.axhline(cw, ls="--"); ax.plot([cx],[cw],"o")

        # (R·Ic_rms/Vpk) vs ωRC
        ax = self.ax_ab[1]; self._clear(ax)
        ax.set_xlabel(LBL_wRC); ax.set_ylabel(LBL_RIcVpk)
        Xw, Yr = _clean_xy(df["omegaRC"], df["R_Ic_eff_over_Vpk"], None, None, yabs=1e2)
        ax.plot(Xw, Yr, lw=2)
        self._set_lims(ax, "R_Icrms_over_Vpk_vs_wRC")
        if w_point is not None:
            # Aqui só guia horizontal/vertical em w_point previsto; a curva já está limitada por set_lims
            y10 = self._interp_safe(w_point, df["omegaRC"], df["R_Ic_eff_over_Vpk"])
            # Clampa a posição visual do ponto (w, y10)
            mode = self._ab_mode()
            w_rng = AB_LIMITS[mode]["R_Icrms_over_Vpk_vs_wRC"]["x"]
            y_rng = AB_LIMITS[mode]["R_Icrms_over_Vpk_vs_wRC"]["y"]
            cw = self._clamp_to_range(w_point, w_rng)
            cy = self._clamp_to_range(y10, y_rng)
            if (cw != w_point) or (cy != y10):
                self._add_warning("Ponto (ωRC, R·Ic/Vpk) ajustado aos limites do gráfico.")
            ax.axvline(cw, ls="--"); ax.axhline(cy, ls="--"); ax.plot([cw],[cy],"o")

        # THD vs x
        ax = self.ax_ab[2]; self._clear(ax)
        ax.set_xlabel(LBL_x); ax.set_ylabel(LBL_THD)
        Xt, Yt = _clean_xy(df["Vcmin_over_Vpk"], df["THD_pu"], None, None, yabs=10)
        ax.plot(Xt, Yt, lw=2)
        self._set_lims(ax, "THD_vs_RIc")
        if x_point is not None:
            y27 = self._interp_safe(x_point, df["Vcmin_over_Vpk"], df["THD_pu"])
            mode = self._ab_mode()
            x_rng = AB_LIMITS[mode]["THD_vs_RIc"]["x"]
            y_rng = AB_LIMITS[mode]["THD_vs_RIc"]["y"]
            cx = self._clamp_to_range(x_point, x_rng)
            cy = self._clamp_to_range(y27, y_rng)
            if (cx != x_point) or (cy != y27):
                self._add_warning("Ponto (x, THD) ajustado aos limites do gráfico.")
            ax.axvline(cx, ls="--"); ax.axhline(cy, ls="--"); ax.plot([cx],[cy],"o")

        # phi1 vs x
        ax = self.ax_ab[3]; self._clear(ax)
        ax.set_xlabel(LBL_x); ax.set_ylabel(LBL_phi1)
        Xp, Yp = _clean_xy(df["Vcmin_over_Vpk"], df["phi1_deg"], None, None, yabs=360)
        ax.plot(Xp, Yp, lw=2)
        self._set_lims(ax, "phi1_vs_RIc")
        if x_point is not None:
            y28 = self._interp_safe(x_point, df["Vcmin_over_Vpk"], df["phi1_deg"])
            mode = self._ab_mode()
            x_rng = AB_LIMITS[mode]["phi1_vs_RIc"]["x"]
            y_rng = AB_LIMITS[mode]["phi1_vs_RIc"]["y"]
            cx = self._clamp_to_range(x_point, x_rng)
            cy = self._clamp_to_range(y28, y_rng)
            if (cx != x_point) or (cy != y28):
                self._add_warning("Ponto (x, φ1) ajustado aos limites do gráfico.")
            ax.axvline(cx, ls="--"); ax.axhline(cy, ls="--"); ax.plot([cx],[cy],"o")

        # FP vs x
        ax = self.ax_ab[4]; self._clear(ax)
        ax.set_xlabel(LBL_x); ax.set_ylabel(LBL_PF)
        Xf, Yf = _clean_xy(df["Vcmin_over_Vpk"], df["PF"], None, None, yabs=1.0)
        ax.plot(Xf, Yf, lw=2)
        self._set_lims(ax, "FP_vs_RIc")
        if x_point is not None:
            y29 = self._interp_safe(x_point, df["Vcmin_over_Vpk"], df["PF"])
            mode = self._ab_mode()
            x_rng = AB_LIMITS[mode]["FP_vs_RIc"]["x"]
            y_rng = AB_LIMITS[mode]["FP_vs_RIc"]["y"]
            cx = self._clamp_to_range(x_point, x_rng)
            cy = self._clamp_to_range(y29, y_rng)
            if (cx != x_point) or (cy != y29):
                self._add_warning("Ponto (x, FP) ajustado aos limites do gráfico.")
            ax.axvline(cx, ls="--"); ax.axhline(cy, ls="--"); ax.plot([cx],[cy],"o")

        try:
            self.canvas_ab.draw_idle()
        except Exception:
            pass

    # ---------- Lógica de projeto + simulação ----------
    def run_project(self):
        # Limpa avisos desta execução
        self._warnings = []
        self._hide_banner()

        topology = self.topology.get()
        f = max(self.f.get(), 1e-3)
        Vph = max(self.Vph.get(), 0.0)
        I   = max(self.I.get(), 0.0)
        ripple = max(self.ripple.get(), 1e-9)
        w = 2*np.pi*f

        if topology == "Monofásico (ponte)":
            self.df_ab = _load_abacos(self.topology.get())

            Vpk = np.sqrt(2)*Vph
            R   = Vpk / max(I, 1e-9)
            x   = max(0.0, min(1.0, 1.0 - ripple/max(Vpk,1e-12)))

            df  = self.df_ab
            omegaRC = self._interp_clip(x, df["Vcmin_over_Vpk"], df["omegaRC"], label="x = V_{L,min}/V_{pk}")
            C = omegaRC / (w * R)

            # Guarda ponto no ábaco e redesenha (clamp visual feito no redraw)
            self._ab_point = (x, omegaRC)
            self.redraw_abacos(*self._ab_point)

            # Simula com o C do ábaco
            t, vsrc, vout, iout, id2, R_used, mets = simulate_single_phase(f, Vph, C, I)

            # Ábacos (valores)
            y10 = self._interp_clip(omegaRC, df["omegaRC"], df["R_Ic_eff_over_Vpk"], label="ωRC")
            Ic_rms_abaco = y10 * Vpk / max(R,1e-12)
            THD_ab  = self._interp_clip(x, df["Vcmin_over_Vpk"], df["THD_pu"],   label="x = V_{L,min}/V_{pk}")
            phi1_ab = self._interp_clip(x, df["Vcmin_over_Vpk"], df["phi1_deg"], label="x = V_{L,min}/V_{pk}")
            PF_ab   = self._interp_clip(x, df["Vcmin_over_Vpk"], df["PF"],       label="x = V_{L,min}/V_{pk}")

            # Resultados
            self.res_vars["Topologia"].set("Topologia: Monofásico (ponte)")
            self.res_vars["Entrada"].set(f"Entrada: {Vph:.3f} V_fase,rms @ {f:.3f} Hz")
            self.res_vars["x = $V_{L,\\min}/V_{pk}$"].set(f"x = V_(L_min)/V_pk = {x:.4f}")
            self.res_vars["ωRC (ábaco)"].set(f"ωRC = {omegaRC:.3f}")
            self.res_vars["C (dimensionado)"].set(f"C = {C*1e3:.3f} mF ({C:.6f} F)")
            self.res_vars["R (adotado)"].set(f"R ≈ {R:.3f} Ω")
            self.res_vars["Ic_rms (ábaco)"].set(f"Ic_rms (ábaco) ≈ {Ic_rms_abaco:.3f} A")
            self.res_vars["THD (ábaco)"].set(f"THD (ábaco) ≈ {THD_ab:.3f} pu")
            self.res_vars["φ1 (ábaco)"].set(f"φ1 (ábaco) ≈ {phi1_ab:.2f} °")
            self.res_vars["FP (ábaco)"].set(f"FP (ábaco) ≈ {PF_ab:.3f}")
            self.res_vars["Vout_médio"].set(f"V̄_out = {mets['V_out_avg']:.3f} V")
            self.res_vars["Vout_RMS"].set(f"V_out,rms = {mets['V_out_rms']:.3f} V")
            self.res_vars["Ripple_sim"].set(f"ΔV_pp = {mets['Ripple_pp']:.3f} V")
            self.res_vars["I_diodo2_médio"].set(f"I_D2, méd = {mets['Id2_avg']:.3f} A")
            self.res_vars["I_diodo2_RMS"].set(f"I_D2, rms = {mets['Id2_rms']:.3f} A")

            # Plots
            self._plot_all(t, vsrc, vout, iout, id2)

        else:
            self.df_ab = _load_abacos(self.topology.get())

            Vll_pk = np.sqrt(2)*np.sqrt(3)*Vph
            R = Vll_pk / max(I, 1e-9)
            x   = max(0.0, min(1.0, 1.0 - ripple/max(Vll_pk,1e-12)))

            df  = self.df_ab
            omegaRC = self._interp_clip(x, df["Vcmin_over_Vpk"], df["omegaRC"], label="x = V_{L,min}/V_{pk}")
            C = omegaRC / (w * R)

            # Guarda ponto no ábaco e redesenha (clamp visual feito no redraw)
            self._ab_point = (x, omegaRC)
            self.redraw_abacos(*self._ab_point)

            # Simula com o C do ábaco
            t, vsrc, vout, iout, id2, R_used, mets = simulate_three_phase(f, Vph, C, I)

            # Ábacos (valores)
            y10 = self._interp_clip(omegaRC, df["omegaRC"], df["R_Ic_eff_over_Vpk"], label="ωRC")
            Ic_rms_abaco = y10 * Vll_pk / max(R,1e-12)
            THD_ab  = self._interp_clip(x, df["Vcmin_over_Vpk"], df["THD_pu"],   label="x = V_{L,min}/V_{pk}")
            phi1_ab = self._interp_clip(x, df["Vcmin_over_Vpk"], df["phi1_deg"], label="x = V_{L,min}/V_{pk}")
            PF_ab   = self._interp_clip(x, df["Vcmin_over_Vpk"], df["PF"],       label="x = V_{L,min}/V_{pk}")

            # Resultados
            self.res_vars["Topologia"].set("Topologia: Trifásico 6 pulsos")
            self.res_vars["Entrada"].set(f"Entrada: {Vph:.3f} V_fase,rms @ {f:.3f} Hz")
            self.res_vars["x = $V_{L,\\min}/V_{pk}$"].set(f"x = V_(L_min)/V_pk = {x:.4f}")
            self.res_vars["ωRC (ábaco)"].set(f"ωRC = {omegaRC:.3f}")
            self.res_vars["C (dimensionado)"].set(f"C = {C*1e3:.3f} mF ({C:.6f} F)")
            self.res_vars["R (adotado)"].set(f"R ≈ {R:.3f} Ω")
            self.res_vars["Ic_rms (ábaco)"].set(f"Ic_rms (ábaco) ≈ {Ic_rms_abaco:.3f} A")
            self.res_vars["THD (ábaco)"].set(f"THD (ábaco) ≈ {THD_ab:.3f} pu")
            self.res_vars["φ1 (ábaco)"].set(f"φ1 (ábaco) ≈ {phi1_ab:.2f} °")
            self.res_vars["FP (ábaco)"].set(f"FP (ábaco) ≈ {PF_ab:.3f}")
            self.res_vars["Vout_médio"].set(f"V̄_out = {mets['V_out_avg']:.3f} V")
            self.res_vars["Vout_RMS"].set(f"V_out,rms = {mets['V_out_rms']:.3f} V")
            self.res_vars["Ripple_sim"].set(f"ΔV_pp = {mets['Ripple_pp']:.3f} V")
            self.res_vars["I_diodo2_médio"].set(f"I_D2, méd = {mets['Id2_avg']:.3f} A")
            self.res_vars["I_diodo2_RMS"].set(f"I_D2, rms = {mets['Id2_rms']:.3f} A")

            self._plot_all(t, vsrc, vout, iout, id2)

        # Mostra (ou oculta) a faixa de aviso
        self._flush_warnings()

    def _plot_all(self, t, vsrc, vout, iout, id2=None):
        # Tensões
        self.ax_top.clear(); self.ax_top.grid(True)
        self.ax_top.set_title(TTL_top)
        self.ax_top.set_xlabel(LBL_t); self.ax_top.set_ylabel(LBL_V)
        self.ax_top.plot(t, vsrc, label="$v_{out,ret}$")
        self.ax_top.plot(t, vout, label="$v_{out}$")
        self.ax_top.legend(loc="upper right", frameon=True, framealpha=0.9)
        self.canvas_top.draw_idle()

        # i_load
        self.ax_mid.clear(); self.ax_mid.grid(True)
        self.ax_mid.set_title(TTL_mid); self.ax_mid.set_xlabel(LBL_t); self.ax_mid.set_ylabel(LBL_I)
        self.ax_mid.plot(t, iout, label="i_load")
        self.canvas_mid.draw_idle()

        # i_D2, se houver
        self.ax_bot.clear(); self.ax_bot.grid(True)
        if id2 is None:
            self.ax_bot.set_title(TTL_bot_na)
        else:
            self.ax_bot.set_title(TTL_bot)
            self.ax_bot.plot(t, id2, label="i_D2")
        self.ax_bot.set_xlabel(LBL_t); self.ax_bot.set_ylabel(LBL_I)
        self.canvas_bot.draw_idle()


if __name__ == "__main__":
    App().mainloop()
