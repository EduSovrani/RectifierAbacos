param(
    [string]$Name = "RectifierAbacos",
    [string]$Icon = "assets\icon.ico",
    [switch]$UseUPX,
    [string]$UpxDir = ""
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$toolsDir  = $PSScriptRoot
$rootDir   = Split-Path $toolsDir -Parent
$venvAct   = Join-Path $rootDir ".venv\Scripts\Activate.ps1"
$createEnv = Join-Path $toolsDir "create_env.ps1"

# nomes dos arquivos (ATUAIS) em src/
$mainScript = "src\rectifier_abacos_gui.py"
$coreModule = "src\rectifier_abacos_core.py"

$mainAbs = Join-Path $rootDir $mainScript
$coreAbs = Join-Path $rootDir $coreModule
$iconAbs = Join-Path $rootDir "assets\icon.ico"

Write-Host "[1/8] Verificando ambiente virtual (.venv)..."
if (-not (Test-Path $venvAct)) {
    if (Test-Path $createEnv) {
        Write-Host "(.venv não encontrado) Criando com tools\create_env.ps1..."
        powershell -ExecutionPolicy Bypass -File $createEnv
    } else { Write-Error "(.venv) não encontrado e tools\create_env.ps1 ausente."; exit 1 }
    if (-not (Test-Path $venvAct)) { Write-Error "Falha ao criar .venv."; exit 1 }
}
. $venvAct
Write-Host "[OK] Ambiente virtual ativado."

Write-Host "[2/8] Instalando/atualizando deps..."
python -m pip install --upgrade pip
python -m pip install pyinstaller numpy matplotlib pandas

Write-Host "[3/8] Limpando build/, dist/ e .spec..."
Set-Location $rootDir
Remove-Item -Recurse -Force build, dist -ErrorAction SilentlyContinue
Remove-Item -Force *.spec -ErrorAction SilentlyContinue
Remove-Item -Force (Join-Path $toolsDir "*.spec") -ErrorAction SilentlyContinue

Write-Host "[4/8] Preparando argumentos do PyInstaller..."
$patterns = @("*.py","*.png","*.ico","*.csv","*.json","*.txt")

$args = @(
    '--clean',
    '--onedir',
    '--noconsole',
    '--noconfirm',
    '--name', $Name,
    '--hidden-import', 'matplotlib.backends.backend_tkagg',
    '--collect-submodules', 'pandas',
    '--specpath', 'tools'
)

if (Test-Path $iconAbs) { $args += @('--icon', $iconAbs) }
else { Write-Host "(!) Ícone não encontrado: $iconAbs — sem --icon." }

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PONTO CRÍTICO: colocar o core exatamente AO LADO da GUI no bundle.
# No onefile isso vira: ...\_MEIxxxxx\rectifier_abacos_core.py
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (Test-Path $coreAbs) {
    $args += @('--add-data', "$coreAbs;.")
}

# Excludes
$excludes = @(
    'matplotlib.tests',
    'matplotlib.backends._macosx',
    'matplotlib.backends.backend_qt5',
    'matplotlib.backends.backend_qt4',
    'matplotlib.backends.backend_qtagg',
    'matplotlib.backends.qt_compat',
    'PyQt5','PySide2','PySide6',
    'tkinter.test',
    'numpy.tests',
    'pandas.tests'
)
foreach ($m in $excludes) { $args += @('--exclude-module', $m) }

# Incluir outros arquivos úteis, mas:
# - NÃO incluir o main (já vai como script)
# - NÃO incluir o core de novo (evita cair em 'src' por engano)
$included = @()
foreach ($pat in $patterns) {
    Get-ChildItem -Path $rootDir -Include $pat -Recurse -File |
        Where-Object {
            $_.FullName -ne $mainAbs -and
            $_.FullName -ne $coreAbs -and
            ($_.FullName -notlike (Join-Path $rootDir "tools\*")) -and
            ($_.FullName -notlike (Join-Path $rootDir "build\*")) -and
            ($_.FullName -notlike (Join-Path $rootDir "dist\*")) -and
            ($_.FullName -notlike (Join-Path $rootDir ".venv\*"))
        } |
        ForEach-Object {
            $originAbs = $_.FullName
            $relPath   = $originAbs.Substring($rootDir.Length + 1)
            $destDir   = Split-Path $relPath -Parent
            if (-not $destDir) { $destDir = "." }
            $args += @('--add-data', "$originAbs;$destDir")
            $included += $relPath
        }
}
Write-Host "Arquivos incluídos (add-data):"
$included | Sort-Object -Unique | ForEach-Object { Write-Host "  $_" }

# script principal por último (ABSOLUTO)
$args += $mainAbs

Write-Host "[5/8] Gerando executável (onefile) com PyInstaller..."
& pyinstaller @args

Write-Host "[6/8] Verificando saída..."
$exePath = Join-Path $rootDir ("dist\" + $Name + ".exe")
if (Test-Path $exePath) {
    $sizeMB = [Math]::Round((Get-Item $exePath).Length / 1MB, 2)
    Write-Host ("   OK: {0}  ({1} MB)" -f $exePath, $sizeMB)
} else {
    Write-Host "   Aviso: dist\$Name.exe não encontrado (veja os logs)."
}

Write-Host "[7/8] Caminho de saída:"
Write-Host ("   {0}" -f $exePath)

Write-Host "[8/8] Finalizado!"
