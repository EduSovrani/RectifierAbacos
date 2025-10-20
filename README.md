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

## Practical notes & limitations

* **Start-up/inrush.** With a large **C**, the bridge can see very high turn-on surge. Use NTC/series resistor (bypass later) or an active pre-charge.
* **Capacitor selection.** Check **ripple current rating** with `I_Cef`, ESR heating, and lifetime at the expected hotspot temperature.
* **EMI and harmonics.** Although 6-pulse is smoother than single-phase, THD can still be high at small ripple (see Fig. 10.32). Consider an input reactor or an active front-end if PF/THD limits apply.
* **Semiconductor drops & parasitics.** Diode forward drop and transformer leakage slightly reduce the actual `V_Cmin` compared with the ideal formulas—give margin.
* **Model scope.** All equations above assume **no line inductance** and **resistive load**. If the load is another converter (dynamic), or if there’s line inductance, conduction angles and spectra change.