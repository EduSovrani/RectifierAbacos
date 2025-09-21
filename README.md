# RectifierAbacos

Tool for **design & simulation** of rectifiers with capacitive filters (single-phase and three-phase), plus **abacus** generation/lookup.

## Folder structure

```text
RectifierAbacos/
├─ assets/
│  └─ icon.ico                    # icon used by the EXE and the window
├─ dist/
│  └─ RectifierAbacos_vX.Y.exe    # built executable (after running the build)
├─ src/
│  ├─ rectifier_abacos_gui.py     # Tkinter application (main GUI)
│  └─ rectifier_abacos_core.py    # computation/abacus core
├─ tools/
│  ├─ build_exe_onefile.ps1       # build script (PyInstaller onefile)
│  └─ create_env.ps1              # creates/updates the virtual env (.venv)
├─ requirements.txt               # Python dependencies (pip)
└─ VERSION                        # app version (e.g., 1.0, 1.1, 1.2.3)
```

### What lives where

| Path                           | What         | Notes                                                                                     |
| ------------------------------ | ------------ | ----------------------------------------------------------------------------------------- |
| `assets/icon.ico`              | App icon     | Embedded into the EXE and shown in the GUI window.                                        |
| `dist/`                        | Build output | After building you’ll see `RectifierAbacos_vX.Y.exe`. You may commit it if you want.      |
| `src/rectifier_abacos_gui.py`  | GUI          | Shows the version in the window title (read from `VERSION`).                              |
| `src/rectifier_abacos_core.py` | Core         | Abacus math and helpers.                                                                  |
| `tools/build_exe_onefile.ps1`  | Build        | Produces a **single-file** EXE and injects Windows file version info.                     |
| `tools/create_env.ps1`         | Env          | Creates `.venv` and installs dependencies.                                                |
| `requirements.txt`             | Pip          | Used by the env scripts to install deps.                                                  |
| `VERSION`                      | Version      | Controls the EXE name (`RectifierAbacos_vX.Y.exe`) and the Windows file version metadata. |

> **Note:** `dist/` appears **only after** you build.

## Environment & build

> Prereqs: **Windows + PowerShell**, Python 3.x.

### Option A — Build the one-file EXE (auto-creates the env if missing)

You **don’t** need to create the virtual environment beforehand. If `.venv` is missing, the build script will call `tools/create_env.ps1` automatically.

```powershell
pwsh -File tools/build_exe_onefile.ps1
```

* Output: `dist/RectifierAbacos_vX.Y.exe`
* The EXE’s Windows *Properties → Details → File version* comes from `VERSION`.
* The GUI window title shows `Retificador + Filtro C — Ábacos (vX.Y)`.

### Option B — Create the environment **only** (for local runs)

If you want to run the GUI from source without building:

```powershell
# Create/refresh the local virtual environment
pwsh -File tools/create_env.ps1

# Activate it (PowerShell)
.\.venv\Scripts\Activate.ps1

# Run the app
python src\rectifier_abacos_gui.py
```

> To update dependencies later, just re-run `tools/create_env.ps1`.

## Versioning the app

* Put the desired version in the root **`VERSION`** file, e.g.:

  ```
  1.0
  ```
* The build will:

  * name the binary `RectifierAbacos_v1.0.exe`, and
  * embed `1.0` into the EXE’s Windows version metadata.

---

# Single-Phase Full-Wave Rectifier with Capacitor Filter — Complete Theory

> **Source used (paraphrased from your photos):** Ivo Barbi, *Eletrônica de Potência*, Chap. 10 (figs. 10.1–10.12; 10.25–10.29; 10.10).
> Notation $x$ and $Y$ matches the app/abaci.

## 0) Circuit, symbols, and assumptions

**Bridge D1–D4 → C → R**.
Input: $v_{ac}(t)=V_{pk}\sin(\omega t)$ with $V_{pk}=\sqrt{2}\,V_{\mathrm{rms}}$.
Capacitor voltage: $v_C(t)$; its minimum per cycle: $V_{C\min}$.

**Key normalized parameters**

* Ripple ratio (also the x-axis of the abaci):

$$
x=\frac{V_{C\min}}{V_{pk}}\in(0,1),\qquad
\Delta V_{pp}\approx V_{pk}-V_{C\min}=V_{pk}(1-x).
$$

* Adimensional group:

$$
Y=\omega RC.
$$

* Angles (one electrical cycle, $\theta=\omega t$):
  $\theta_1$ = **start of next recharge** (diodes turn on),
  $\theta_2$ = **end of recharge** (input current goes to zero),
  $\theta_3$ = start of the previous ramp. Common definitions:

$$
\alpha=\tfrac{3\pi}{2}-\theta_1,\quad
\beta=\theta_2-\tfrac{\pi}{2},\quad
\gamma=\theta_1-\theta_2,\quad
\alpha+\beta+\gamma=\pi.
$$

**Assumptions:** ideal diodes, **no** line inductance, resistive load, steady state.

---

## 1) Qualitative operation (three stages per half-cycle)

Let $\theta=\omega t$.

1. **S1 – Recharge** $(\theta_3\le \theta\le \pi/2)$
   Diodes conduct; the source **clamps** $v_C$: $v_C=V_{pk}\sin\theta$.

$$
i_C=C\frac{dv_C}{dt}=\omega C V_{pk}\cos\theta,\qquad
i_R=\frac{v_C}{R}=\frac{V_{pk}}{R}\sin\theta,
$$

$$
i_{\text{line}}=i_1=i_C+i_R.
$$

2. **S2 – Conduction tail** $(\pi/2<\theta<\theta_2)$
   $i_C$ becomes negative and exactly **cancels** $i_R$ at $\theta=\theta_2$ ⇒ diodes **turn off**.

3. **S3 – Free discharge** $(\theta_2\le \theta\le \theta_1)$
   All diodes **block**; $C$ alone feeds $R$. $v_C$ decays down to $V_{C\min}$.
   Next recharge begins when the sinusoid **catches up** with the decayed capacitor voltage: $\theta=\theta_1$.

---

## 2) Conduction window — quick vs exact boundaries

### 2.1 **Quick** boundary (handy approximation)

During S1/S2, take $i_R\approx V_{pk}/R$ **constant**. The zero-current condition at $\theta_2$ gives

$$
0=\omega C V_{pk}\cos\theta_2+\frac{V_{pk}}{R}
\quad\Rightarrow\quad
\boxed{\cos\theta_2\approx-\frac{1}{Y}}\qquad(Y>1).
$$

### 2.2 **Exact** boundary (detailed model)

Using the true $i_R=\tfrac{V_{pk}}{R}\sin\theta$ at $i_1(\theta_2)=0$:

$$
\omega C V_{pk}\cos\theta_2+\frac{V_{pk}}{R}\sin\theta_2=0
\ \Rightarrow\
\boxed{\tan\theta_2=-Y,\ \ \theta_2\in(\tfrac{\pi}{2},\pi)}.
$$

Other geometric relations:

$$
\boxed{\theta_1=\pi+\arcsin(x)},\qquad
\alpha=\tfrac{\pi}{2}-\arcsin x,\quad
\beta=\tfrac{\pi}{2}-\arctan Y,\quad
\gamma=\pi-\alpha-\beta.
$$

---

## 3) Energy balance (two equivalent routes)

### 3.1 **Route A — One-line energy sizing**

Energy **delivered by $C$** to the load every **half cycle**:

$$
\boxed{W_C=\tfrac{1}{2}C\left(V_{pk}^2-V_{C\min}^2\right)}.
$$

Steady state ⇒ match load energy for the same interval:

$$
\boxed{W_C=\frac{P_{out}}{2f}}
\quad\Rightarrow\quad
\boxed{C=\frac{P_{out}}{f\left(V_{pk}^2-V_{C\min}^2\right)}}.
$$

> Excellent for **fast capacitor sizing** from $P_{out}$, $f$, and the chosen ripple $x=V_{C\min}/V_{pk}$.

### 3.2 **Route B — Area/angle model (source of the abacus $Y=f(x)$)**

Split the transferred charge into the **three slices** S1–S3 (see the classic sketch):

* **S1** (recharge slice):

$$
\boxed{S_1=\omega C V_{pk}\left(1-\cos\alpha\right)}.
$$

* **S2** (load discharge while diodes are off):

$$
\boxed{S_2=\frac{\omega RC V_{pk}\cos\beta}{R}\left(1-e^{\gamma/Y}\right)}.
$$

* **S3** (small triangular correction around $\theta_2$):

$$
\boxed{S_3=\frac{\beta V_{pk}\cos\beta}{2R}}.
$$

Charge balance over a half cycle is **zero**: $S_1-|S_2|-|S_3|=0$. Gathering terms:

$$
\boxed{Y (1-\cos\alpha)-\frac{\beta\cos\beta}{2}-Y\cos\beta\left(1-e^{\gamma/Y}\right)=0}
$$

with $\alpha,\beta,\gamma$ functions of $x$ and $Y$ as above.
For a **given $x$**, solve **numerically** to obtain the **abacus**:

$$
\boxed{Y=\omega RC=f(x)}\quad\text{(Fig. 10.9)}.
$$

Then compute **$C=\dfrac{Y}{\omega R}$**.

> **Which route to prefer?**
>
> * **Energy one-liner**: quick $C$ sizing from $P_{out}$.
> * **Areas/angles**: ensures **coherence** with conduction window, peaks, $I_{C,\mathrm{rms}}$, THD, PF — and underlies the **abaci**.

---

## 4) Peak recharge current and effective window

A practical “equivalent charge” estimate for the peak is

$$
\boxed{I_p\approx\frac{C\Delta V}{t_c}},\qquad \Delta V\approx V_{pk}-V_{C\min}.
$$

> The book notes the **actual peak** (cosine/triangular-shaped) is typically **≈ 2× higher** than the simple **rectangular** estimate for the same $\Delta Q$. Use this as a **sizing margin** for diodes/transformer.

With the **exact** boundary, the conduction time follows directly from $\theta_2,\theta_1$ (or from $x$ via $\theta_1$).

---

## 5) Capacitor RMS current $I_{C,\mathrm{rms}}$

### 5.1 **Quick** rectangular approximation

Model the recharge as a rectangle of height $I_p$ and width $t_c$:

$$
\boxed{I_{C\mathrm{rms}}\approx I_p\sqrt{t_c f}}.
$$

> Good for **order-of-magnitude**. Remember the **≈2×** factor for the **true peak**.

### 5.2 **Exact/normalized** form (leads to **Fig. 10.10 abacus**)

Split the RMS into **conduction** and **block** contributions and normalize by $V_{pk}/R$:

$$
\boxed{\left(\frac{R I_{C\mathrm{rms}}}{V_{pk}}\right)^2=
\underbrace{\frac{Y^2}{\pi}\int_{\theta_3}^{\theta_2}\cos^2\theta d\theta}_{\text{conduction}}
+
\underbrace{\frac{\cos^2\beta}{\pi} \int_{0}^{\gamma}e^{2\theta/Y} d\theta}_{\text{block}}}.
$$

With

$$
\boxed{\theta_1=\pi+\arcsin x,\quad \theta_2=\pi+\arctan(-Y),\quad \theta_3=\arcsin x},
$$

and the closed form

$$
\int_{0}^{\gamma} e^{2\theta/Y} d\theta=\frac{Y}{2}\left(e^{2\gamma/Y}-1\right),
$$

you finally get the **normalized** value

$$
\boxed{\frac{R\,I_{C,\mathrm{rms}}}{V_{pk}}=F(Y)}
\quad\Rightarrow\quad
\boxed{I_{C,\mathrm{rms}}=\frac{V_{pk}}{R}\,F(Y)}.
$$

The book plots $F(Y)$ in **Fig. 10.10**. This is essential for **ESR/thermal** checks.

---

## 6) Input current spectrum, **THD**, advance $\varphi_1$, and **PF**

Inside the slice $[\theta_3,\theta_2]$:

$$
i_1(\theta)=\omega C V_{pk}\cos\theta+\frac{V_{pk}}{R}\sin\theta,\qquad
i_1(\theta)=0\ \text{outside that slice (half-cycle view)}.
$$

Half-wave symmetry ⇒ only **odd harmonics**. Using the integrals

$$
A1=\int_{\theta_3}^{\theta_2}\cos\theta\,\cos(n\theta) d\theta,\quad
A2=\int_{\theta_3}^{\theta_2}\sin(n\theta) d\theta,
$$

$$
B1=\int_{\theta_3}^{\theta_2}\cos\theta \sin(n\theta) d\theta,\quad
B2=\int_{\theta_3}^{\theta_2}\sin\theta \sin(n\theta) d\theta,
$$

and the book’s normalization:

$$
\bar a_n=Y\cdot A1+A2,\qquad \bar b_n=Y\cdot B1+B2,\qquad
c_n=\sqrt{\bar a_n^2+\bar b_n^2}.
$$

you obtain:

* **Harmonic amplitudes** $c_n(x)$ — Figs. **10.25** and **10.26**.
* **THD** vs $x$ — Fig. **10.27** (**rises** as ripple shrinks, i.e., $x\to1$).
* **Advance angle** $\varphi_1(x)$ — Fig. **10.28** (fundamental current **leads** the voltage).
* **Power factor (PF)** — Fig. **10.29** (**decreases** as $x$ increases).

**Qualitative takeaway**: larger $C$ ⇒ **better DC** mas **corrente de linha mais aguda** ⇒ **THD↑**, **PF↓**, $\varphi_1$ capacitiva.

---

## 7) Capacitor design — two side-by-side recipes

### (A) **Energy one-liner** (you know $P_{out}$)

1. Specify $V_{rms}$, $f$, $P_{out}$ e alvo $x=V_{C\min}/V_{pk}$ (ou $\Delta V_{pp}$).
2. $V_{pk}=\sqrt{2}V_{rms}$, $x=1-\Delta V_{pp}/V_{pk}$.

$$
\boxed{C=\dfrac{P_{out}}{f\left(V_{pk}^2-V_{C\min}^2\right)}}.
$$

4. Sanity-check com a janela de condução (rápida ou exata) e com $I_p$.

### (B) **Full abacus flow** (you know $R$ or $I$)

1. Especifique $V_{rms}$, $f$, $R$ (ou $I$) e $x$.
2. Leia $Y=\omega RC$ de **$Y=f(x)$** (Fig. 10.9) ⇒

$$
\boxed{C=\dfrac{Y}{\omega R}}.
$$

3. Com o **mesmo $Y$**, use **Fig. 10.10** para obter

$$
\boxed{I_{C,\mathrm{rms}}=\dfrac{V_{pk}}{R} F(Y)}.
$$

4. Se preciso, use Figs. **10.27–10.29** no seu $x$ para **THD / $\varphi_1$ / PF**.
5. Verifique picos em diodos/trafo (lembre do **≈2×** vs retangular).

---

## 8) Practical notes & limitations

* Real **line inductance** widens the conduction interval and **reduces** peak/THD.
* For a **non-purely resistive** load (e.g., downstream converter), use the **equivalent impedance** at the ripple frequency.
* Outside the abacus range, solve the implicit $Y$ equation **numerically** (or **simulate**).
* Here we treated **bridge + single capacitor**. The book also covers the **voltage doubler**; the energy ideas carry over, but the angle relations change.

---

## 9) Quick formula sheet

* $x=\dfrac{V_{C\min}}{V_{pk}}$, $Y=\omega RC$.
* $\theta_1=\pi+\arcsin x$, **exato** $\tan\theta_2=-Y$; **rápido** $\cos\theta_2\approx-1/Y$.
* $\alpha=\tfrac{\pi}{2}-\arcsin x,\ \beta=\tfrac{\pi}{2}-\arctan Y,\ \gamma=\pi-\alpha-\beta$.
* **Capacitor from energy**: $\displaystyle C=\dfrac{P_{out}}{f\left(V_{pk}^2-V_{C\min}^2\right)}$.
* **Capacitor from abacus**: leia $Y(x)$ ⇒ $C=\dfrac{Y}{\omega R}$.
* **Capacitor RMS**: $\displaystyle \frac{R\,I_{C,\mathrm{rms}}}{V_{pk}}=F(Y)$ (Fig. 10.10).
* **Peak estimate**: $I_p\approx\dfrac{C\,\Delta V}{t_c}$ (**true peak ≈ 2×** retangular).
* **Slice current**: $i_1(\theta)=\omega C V_{pk}\cos\theta+\dfrac{V_{pk}}{R}\sin\theta$ dentro de $[\theta_3,\theta_2]$.

---

# Three-Phase Capacitor-Input Rectifier (6-pulse, no line inductance)

*Based on Ivo Barbi, **Eletrônica de Potência**, cap. 10 (figs. 10.16–10.34).*

> It follows the same style as the single-phase section, but now for the three-phase Graetz bridge with a purely capacitive DC filter and **no** series line inductance.

---

## 1) Circuit, assumptions and symbols

* Topology: three-phase Graetz bridge feeding a large capacitor **C** and a resistive load **R** (Fig. 10.16).
* Line voltages (peak): `V_Lp` (e.g., for 380 V\_LL,rms, `V_Lp = √2 · 380 ≃ 535 V`).
* Mains electrical pulsation: `ω = 2πf`.
* Output capacitor voltage: `v_C(θ)` with a **sawtooth ripple** around its minimum `V_Cmin`.
* With no input inductance, each 60° interval (π/3 rad) transfers charge from the AC side to **C** through **two** diodes (Fig. 10.17).
* We parameterize ripple by the ratio `ρ ≜ V_Cmin / V_Lp` (typically 0.87…0.98).

Angles used in the detailed model (see Figs. 10.18 and 10.21):

* `θ1, θ2, θ3` delimit the conduction window in one 60° sector.
* Auxiliary angles:

  * `α` = portion (within the sector) where `i_C` grows linearly (pre-peak),
  * `β` = portion where the source still conducts after the current peak,
  * `γ` = gap with **no** conduction (source blocked).
  * In one 60° sector: `α + β + γ = π/3`.

Closed relations used throughout (from the geometry of Fig. 10.21):

$$
\begin{aligned}
&\frac{V_{C\min}}{V_{Lp}} = \sin\bigl(\theta_1 - \tfrac{\pi}{3}\bigr), &&\Rightarrow\quad 
\theta_1 = \tfrac{\pi}{3} + \sin^{-1}\Bigl(\tfrac{V_{C\min}}{V_{Lp}}\Bigr), \\
&\theta_2 = \pi + \tan^{-1}(-\omega RC), \\
&\theta_3 = \sin^{-1}\Bigl(\tfrac{V_{C\min}}{V_{Lp}}\Bigr), \\
&\alpha = \tfrac{5\pi}{6} - \theta_1 = \tfrac{\pi}{2} - \sin^{-1}\Bigl(\tfrac{V_{C\min}}{V_{Lp}}\Bigr), \\
&\beta = \tfrac{\pi}{2} + \tan^{-1}(-\omega RC), \qquad
\gamma = \tfrac{\pi}{3} - \alpha - \beta .
\end{aligned}
$$

> Rule-of-thumb boundary: the special case **γ = 0** (i.e., no blocking gap inside a sector) happens at
> `ωRC ≃ 1.73` and `V_Cmin/V_Lp ≃ 0.87`.
> Equivalently, for design checks:
> $R \le \dfrac{1.73}{\omega C}.$

---

## 2) Qualitative operation (per 60° sector)

1. **Start of the sector (`θ = θ3`)** – the line-to-line driving that sector exceeds `v_C`; the conducting diode pair turns on; `i_C` rises.
2. **During `α`** – the source charges **C**; with `v_C = V_Lp·\sin θ`, capacitor current is approximately
   $i_C(θ) = \omega C V_{Lp} \cos θ.$
   The load simultaneously draws `i_R = v_C/R`.
3. **Peak** – occurs near `θ ≈ θ1` (inside the sector).
4. **During `β`** – current decays; when the driving line voltage falls below `v_C`, diodes turn off.
5. **During `γ`** – **C** alone feeds **R**; `v_C` decreases almost linearly until the next sector begins.

This repeats six times per mains period (6-pulse rectification). Compared with single-phase, ripple is much smaller for the **same** C because recharge happens every 60° instead of every π radians.

---

## 3) Energy balance and first sizing of C

### 3.1 Simple energetic sizing (recommended first pass)

Within each **60°** sector the source must replenish the energy that the load removed while the diodes were off. From Fig. 10.18 (and Sec. 10.6.B), the **energy delivered to the capacitor each 60° interval** is

$$
E_\text{sec} = \tfrac{1}{2} C \bigl( V_{Lp}^2 - V_{C\min}^2 \bigr).
$$

Over one mains period there are **six** such intervals. Equating **energy provided by the source** to **energy dissipated by the load**:

$$
6 \cdot \tfrac{1}{2} C \bigl( V_{Lp}^2 - V_{C\min}^2 \bigr) = \frac{P_o}{f}.
$$

Therefore,

$$
\boxed{C = \frac{P_o}{6 f (V_{Lp}^2 - V_{C\min}^2)}}
\qquad\text{(energy method, 6-pulse)}
$$

This is the three-phase counterpart of the single-phase “difference-of-squares” sizing and is excellent for design targets stated as a **minimum allowed DC voltage** (`V_Cmin`).
If you prefer to parametrize by `ρ = V_Cmin/V_Lp`:

$$
C = \frac{P_o}{6\ f V_{Lp}^2 (1 - \rho^2)} .
$$

> Tip – If you already chose `R` and the ripple ratio `ρ`, you can also use
> `Y ≜ ωRC` from the abacus of Sec. 4 and then **C = Y/(ωR)**. Both routes are consistent.

---

## 4) Accurate boundaries and the implicit relation for `Y = ωRC`

A more precise model (Sec. 10.6.D) tracks **three** sub-intervals inside each 60° sector (`α, β, γ`). Charge conservation over the sector leads to the implicit **nonlinear** equation (tri-phase analogue of the single-phase one):

$$
\boxed{\omega RC (1-\cos\alpha)-\frac{\beta\cos\beta}{2}-\omega RC \cos\beta \Bigl(1 - e^{(\theta_1-\theta_2)/(\omega RC)}\Bigr)=0}
$$

with the angle ties listed in Sec. 1 and `γ = θ_1 - θ_2 = π/3 - α - β`.
Solving (★) for `Y = ωRC` as a function of `ρ = V_Cmin/V_Lp` produces the **abacus** in Fig. 10.22:

* As `ρ → 0.98`, required `Y` (hence **C**) grows steeply.
* The point `ρ ≃ 0.87` corresponds to `Y ≃ 1.73` (the **γ = 0** boundary).

In code or a spreadsheet, is easily solved with a 1-D root finder (`brentq/newton`).

---

## 5) Conduction window, recharge peak and input peak current

From the geometry:

$$
\alpha = \cos^{-1}\Bigl(\tfrac{V_{C\min}}{V_{Lp}}\Bigr), 
\qquad \Delta t_\text{cond} = \frac{\alpha}{\omega}.
$$

* **Capacitor peak recharge current** (at `θ≈θ1`):

$$
\boxed{i_{CP} = \omega C V_{Lp} \sin \alpha = \omega C \sqrt{V_{Lp}^2 - V_{C\min}^2}}
$$

* **Resistive current peak** at the top of the sector:

$$
i_R^\text{pk} = \frac{V_{Lp}}{R}.
$$

* **Input line peak current** (useful for diode/bridge surge checking):

$$
\boxed{i_p \approx i_{CP} + i_R^\text{pk}}
$$

---

## 6) RMS current in the filter capacitor

Split the RMS calculation in the two conducting sub-intervals within the sector (Fig. 10.23):

$$
\begin{aligned}
I_{C1,\text{ef}}^2 &= \frac{(\omega RC)^2}{\pi}\int_{\theta_3}^{\theta_2}\cos^2\theta d\theta, \\
I_{C2,\text{ef}}^2 &= \frac{\cos^2\beta}{\pi}\int_{0}^{\theta_1-\theta_2} e^{\tfrac{2\theta}{\omega RC}} d\theta, \\
I_{C,\text{ef}}    &= \frac{I_{C\text{ef}}^{\text{(line)}} R}{V_{Lp}} \quad\text{(normalization used in Fig. 10.24)}.
\end{aligned}
$$

Using those expressions yields the parametric curve of **Fig. 10.24**, which gives the normalized RMS:

$$
\boxed{\frac{R I_{C,\text{ef}}}{V_{Lp}} = f(\omega RC)}
\quad\Rightarrow\quad I_{C\text{ef}} = \frac{V_{Lp}}{R} f(\omega RC).
$$

### Special closed form (very useful): case **γ = 0**

When `γ = 0` (i.e., `ωRC ≃ 1.73`, `ρ ≃ 0.87`), integration simplifies to

$$
\boxed{I_{C,\text{ef}} \approx 0.30\,\omega\,C\,V_{Lp}}
$$

— a remarkably accurate shortcut for quick ripple-current rating of **C**.

---

## 7) Harmonic content, THD, phase advance and power factor

Define the line current in phase-**a** during its sector (Fig. 10.30):

$$
i_a(\theta) = \omega C V_p \cos \theta + \frac{V_p}{R}\sin \theta,\qquad \theta\in(\theta_3 \theta_2),
$$

and in phase-**b** (shifted by −π/3):

$$
i_b(\theta) = \omega C V_p \cos\bigl(\theta-\tfrac{\pi}{3}\bigr) + \frac{V_p}{R}\sin\bigl(\theta-\tfrac{\pi}{3}\bigr),
$$

(similar for phase-**c**). Over one period the waveform has **half-wave symmetry**, hence only **odd** harmonics. The Fourier coefficients are evaluated over the two conducting sub-intervals of each phase. Using the book’s normalization,

$$
\bar a_n = \frac{\pi R}{2V} a_n, 
\qquad
\bar b_n = \frac{\pi R}{2V} b_n,
\qquad 
\bar c_n = \sqrt{\bar a_n^2 + \bar b_n^2}.
$$

From these, the following **abaci** are read (Figs. 10.31–10.34) as functions of `ρ = V_Cmin/V_p`:

* **Harmonic amplitudes** $\bar c_n$ for $n = 1,5,7,11,13,17$ (Fig. 10.31).
* **THD** of the line current (up to 41st harmonic) increases as ripple is reduced (Fig. 10.32).
* **Fundamental phase advance** $\phi$ (current leads voltage) has a broad maximum around $\rho\approx0.89$ and then decreases (Fig. 10.33).
* **Power factor** $FP$ degrades as $\rho \to 1$ (Fig. 10.34). Typical values: with $\rho\approx0.90$, $FP\approx0.7$.

**Takeaways.** Compared with single-phase, the 6-pulse rectifier has **much lower output ripple** for the same C, and exhibits **higher PF** for the same ripple target, but the input current still has significant harmonics (n = 5, 7, 11, …).

---

## 8) Design recipe (step-by-step)

1. **Specs.** Choose mains `f` and `V_Lef`, load spec (`P_o` or `R`), and ripple target `ρ = V_Cmin/V_Lp`. Compute `V_Lp = √2·V_Lef`.
2. **First C estimate (energy method).**
   $C = \dfrac{P_o}{6 f (V_{Lp}^2 - V_{C\min}^2)} = \dfrac{P_o}{6 f V_{Lp}^2 (1-\rho^2)}.$
3. **Refine with abacus or solver.**
   Pick `Y = ωRC` from **Fig. 10.22** (or solve (★)). Then **reconcile** with your `R` by setting `C = Y/(ωR)`. Iterate with step 2 if needed to hit the exact `ρ`.
4. **Angles.**
   `α = arccos(ρ)`, `β = π/2 + atan(−Y)`, `γ = π/3 − α − β`. Check `γ ≥ 0`.
5. **Currents.**

   * Recharge peak: `i_CP = ω C V_Lp sin α`.
   * Resistive peak: `i_R^pk = V_Lp/R`.
   * Input peak: `i_p ≈ i_CP + i_R^pk`.
6. **Capacitor RMS (thermal rating).**

   * Quick: if you are near `γ ≈ 0`, use `I_{C,ef} ≈ 0.30 ω C V_Lp`.
   * Precise: read from **Fig. 10.24** or integrate numerically. Choose capacitor series/ESR accordingly.
7. **Harmonics & PF.**
   From **Figs. 10.31–10.34**, read `\bar c_n`, **THD**, `φ` and **PF** for your `ρ`. Verify grid-code/EMI constraints if relevant.
8. **Diodes/bridge selection.**

   * Average and RMS follow from the sector waveform (triangles every 60°).
   * Check surge with `i_p` and repetitive peaks with `i_CP`. Provide snubbers if needed.
9. **Sanity checks.**

   * If `R ≤ 1.73/(ωC)`, the design sits at (or beyond) the **γ = 0** boundary; waveforms match Fig. 10.23.a and the RMS shortcut is very good.
   * If `R` is larger (or **C** smaller), there will be a gap `γ>0` (Fig. 10.23.b), conduction is shorter and peaks grow—**recheck diode/bridge ratings**.

---

## 9) Worked examples (matching the book)

### 9.1 “General” example via Fig. 10.24

Given `V_Lef = 380 V` → `V_Lp ≃ 535 V`, `ω = 377 rad/s`, `C = 100 µF`, `R = 160 Ω`.
Then `ωRC ≃ 6`. Reading **Fig. 10.24** gives `(R·I_Cef)/V_Lp ≃ 1.1`, hence

```
I_Cef ≃ (V_Lp/R)·1.1 ≃ (535/160)·1.1 ≃ 3.7 A_rms
```

(consistency with the book’s result).

### 9.2 Special case **γ = 0**

For `ω = 377 rad/s`, `C = 100 µF`, `V_Lp = 535 V`, using the shortcut:

```
I_Cef ≃ 0.30 · ω · C · V_Lp
     ≃ 0.30 · 377 · 100e-6 · 535 ≃ 6 A_rms
```

(also confirmed in the text when the operating point is near `ρ ≃ 0.87`, `ωRC ≃ 1.73`).

### 9.3 Example with energy sizing + peaks

Suppose `P_o = 1800 W`, `V_Lef = 380 V`, `f = 60 Hz`, target `V_Cmin = 485 V` (`ρ ≃ 0.907`).
Using the **energy method**:

$$
C = \frac{P_o}{6 f (V_{Lp}^2 - V_{C\min}^2)}
  = \frac{1800}{6\cdot60\cdot(535^2 - 485^2)} \ \text{F}.
$$

Angles: `α = arccos(485/535) ≃ 25°`, so the **conduction time** is `Δt = α/ω ≃ 1.17 ms`.
Peaks:
`i_CP = ω C V_Lp sin α`, `i_R^pk = V_Lp/R`, and `i_p ≃ i_CP + i_R^pk` (≈ 12 A in the example from the book).

---

## 10) Practical notes & limitations

* **Start-up/inrush.** With a large **C**, the bridge can see very high turn-on surge. Use NTC/series resistor (bypass later) or an active pre-charge.
* **Capacitor selection.** Check **ripple current rating** with `I_Cef`, ESR heating, and lifetime at the expected hotspot temperature.
* **EMI and harmonics.** Although 6-pulse is smoother than single-phase, THD can still be high at small ripple (see Fig. 10.32). Consider an input reactor or an active front-end if PF/THD limits apply.
* **Semiconductor drops & parasitics.** Diode forward drop and transformer leakage slightly reduce the actual `V_Cmin` compared with the ideal formulas—give margin.
* **Model scope.** All equations above assume **no line inductance** and **resistive load**. If the load is another converter (dynamic), or if there’s line inductance, conduction angles and spectra change.

---

### What you can “read from the abaci”

* **Fig. 10.22** → `ωRC = f(V_Cmin/V_Lp)` (use to refine **C** once `R` is known).
* **Fig. 10.24** → `R·I_Cef/V_Lp` vs `ωRC` (capacitor thermal sizing).
* **Figs. 10.31–10.34** → line-current **harmonics**, **THD**, **phase advance** `φ`, and **PF** vs `V_Cmin/V_p`.
