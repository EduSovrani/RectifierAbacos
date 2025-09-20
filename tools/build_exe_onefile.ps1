param(
    [string]$Name = "RectifierAbacos",
    [string]$Icon = "assets\icon.ico",
    [switch]$UseUPX,
    [string]$UpxDir = ""
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

# --- caminhos base (script está em tools/) ---
$toolsDir  = $PSScriptRoot
$rootDir   = Split-Path $toolsDir -Parent
$venvAct   = Join-Path $rootDir ".venv\Scripts\Activate.ps1"
$createEnv = Join-Path $toolsDir "create_env.ps1"

# --- nomes dos arquivos (atuais) em src/ ---
$mainScript    = "src\rectifier_abacos_gui.py"
$coreModule    = "src\rectifier_abacos_core.py"
$coreTriModule = "src\rectifier_abacos_core_tri.py"   # opcional: só entra se existir

$mainAbs    = Join-Path $rootDir $mainScript
$coreAbs    = Join-Path $rootDir $coreModule
$coreTriAbs = Join-Path $rootDir $coreTriModule
$iconAbs    = Join-Path $rootDir $Icon

# ---------------- [1/8] Verificar/criar venv ----------------
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

# ---------------- (novo) Ler versão e preparar nomes ----------
# Lê a versão do arquivo VERSION (na raiz). Usa "1.0" se não existir.
$versionFile = Join-Path $rootDir "VERSION"
$Version = if (Test-Path $versionFile) { (Get-Content $versionFile -Raw).Trim() } else { "1.0" }
if ($Version -notmatch '^\d+(\.\d+){0,3}$') {
    Write-Host "(!) Versão inválida em VERSION: '$Version' — usando 1.0"
    $Version = "1.0"
}
$SanitizedVersion = ($Version -replace '[^0-9\.]', '-')
$OutputName = "{0}_v{1}" -f $Name, $SanitizedVersion  # nome final do .exe (sem extensão)

# ---------------- [2/8] Dependências de build ---------------
Write-Host "[2/8] Instalando/atualizando deps (pyinstaller + numpy + matplotlib + pandas)..."
python -m pip install --upgrade pip
python -m pip install pyinstaller numpy matplotlib pandas

# ---------------- [3/8] Limpeza anterior -------------------
Write-Host "[3/8] Limpando build/, dist/ e *.spec antigos..."
Set-Location $rootDir
Remove-Item -Recurse -Force build, dist -ErrorAction SilentlyContinue
Remove-Item -Force *.spec -ErrorAction SilentlyContinue
Remove-Item -Force (Join-Path $toolsDir "*.spec") -ErrorAction SilentlyContinue

# ---------------- (novo) Gerar version_info.txt --------------
# Gera em build\version_info.txt (build/ já foi limpa acima)
$buildDir = Join-Path $rootDir "build"
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null

# Converte "1.2" -> (1,2,0,0)
$verParts = $Version.Split('.')
while ($verParts.Count -lt 4) { $verParts += '0' }
$verTuple = $verParts[0..3] | ForEach-Object { [int]$_ }
$FileVersionStr = ($verTuple -join '.')

$versionInfoPath = Join-Path $buildDir "version_info.txt"
@"
# -*- coding: utf-8 -*-
VSVersionInfo(
  ffi=FixedFileInfo(
    filevers=($($verTuple[0]), $($verTuple[1]), $($verTuple[2]), $($verTuple[3])),
    prodvers=($($verTuple[0]), $($verTuple[1]), $($verTuple[2]), $($verTuple[3])),
    mask=0x3f,
    flags=0x0,
    OS=0x40004,
    fileType=0x1,
    subtype=0x0,
    date=(0, 0)
    ),
  kids=[
    StringFileInfo([
      StringTable(
        '041604B0',
        [
          StringStruct('CompanyName',        'Your Company'),
          StringStruct('FileDescription',    'Rectifier + Filtro C — Ábacos e GUI'),
          StringStruct('FileVersion',        '$FileVersionStr'),
          StringStruct('InternalName',       '$OutputName'),
          StringStruct('LegalCopyright',     '© $(Get-Date -Format yyyy)'),
          StringStruct('OriginalFilename',   '$OutputName.exe'),
          StringStruct('ProductName',        '$Name'),
          StringStruct('ProductVersion',     '$FileVersionStr')
        ]
      )]),
    VarFileInfo([VarStruct('Translation', [0x0416, 0x04B0])])
  ]
)
"@ | Set-Content -Encoding UTF8 $versionInfoPath

# ---------------- [4/8] Preparar argumentos ----------------
Write-Host "[4/8] Preparando argumentos do PyInstaller..."
$patterns = @("VERSION","*.py","*.png","*.ico","*.csv","*.json","*.txt")

$args = @(
    '--clean',
    '--onefile',
    '--noconsole',
    '--noconfirm',
    '--name', $OutputName,                 # usa nome com versão
    '--hidden-import', 'matplotlib.backends.backend_tkagg',
    '--collect-data', 'pandas',
    '--collect-submodules', 'pandas._libs',
    '--collect-submodules', 'pandas.core',
    '--version-file', $versionInfoPath,    # caminho ABSOLUTO em build\
    '--specpath', 'tools'
)
if (Test-Path $iconAbs) { $args += "--icon=$iconAbs" } else { Write-Host "(!) Ícone não encontrado: $iconAbs — sem --icon." }

# Ícone (usar ABSOLUTO)
if (Test-Path $iconAbs) {
    $args += @('--icon', $iconAbs)
} else {
    Write-Host "(!) Ícone não encontrado: $iconAbs — sem --icon."
}

# UPX (opcional)
if ($UseUPX) {
    $upxExe = ""
    if ($UpxDir -ne "" -and (Test-Path (Join-Path $UpxDir "upx.exe"))) {
        $upxExe = (Join-Path $UpxDir "upx.exe")
    } else {
        $upxExe = (Get-Command upx -ErrorAction SilentlyContinue)?.Source
    }
    if ($upxExe) {
        $args += @('--upx-dir', (Split-Path $upxExe -Parent))
        Write-Host "[UPX] Usando: $upxExe"
    } else {
        Write-Host "(!) UPX não encontrado — build seguirá sem compressão."
    }
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

# e empacote o arquivo para a GUI achar em runtime:
$args += "--add-data=$(Join-Path $rootDir 'assets\icon.ico');assets"

# PONTO CRÍTICO: colocar o(s) core(s) AO LADO da GUI no bundle (onefile -> _MEI\)
if (Test-Path $coreAbs)    { $args += "--add-data=$coreAbs;." }
if (Test-Path $coreTriAbs) { $args += "--add-data=$coreTriAbs;." }

# empacota o arquivo VERSION na raiz do bundle (._MEI\)
$verAbs = Join-Path $rootDir 'VERSION'
if (Test-Path $verAbs) { $args += "--add-data=$verAbs;." }

# Incluir outros arquivos úteis, mas:
# - NÃO incluir o main
# - NÃO incluir os core(s) novamente
$included = @()
Get-ChildItem -Path $rootDir -Recurse -File -Include $patterns |
    Where-Object {
        $_.FullName -ne $mainAbs -and
        $_.FullName -ne $coreAbs -and
        $_.FullName -ne $coreTriAbs -and
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
        $args += "--add-data=$originAbs;$destDir"
        $included += $relPath
    }

Write-Host "Arquivos incluídos (add-data):"
$included | Sort-Object -Unique | ForEach-Object { Write-Host "  $_" }

# script principal por último (ABSOLUTO)
$args += $mainAbs

# ---------------- [5/8] Build ------------------------------
Write-Host "[5/8] Gerando executável (onefile) com PyInstaller..."
& pyinstaller @args

# ---------------- [6/8] Verificar saída --------------------
Write-Host "[6/8] Verificando saída..."
$exePath = Join-Path $rootDir ("dist\" + $OutputName + ".exe")
if (Test-Path $exePath) {
    $sizeMB = [Math]::Round((Get-Item $exePath).Length / 1MB, 2)
    Write-Host ("   OK: {0}  ({1} MB)" -f $exePath, $sizeMB)
} else {
    Write-Host ("   Aviso: dist\{0}.exe não encontrado (veja os logs)." -f $OutputName)
}

# ---------------- [7/8] Caminho de saída -------------------
Write-Host "[7/8] Caminho de saída:"
Write-Host ("   {0}" -f $exePath)

# ---------------- [8/8] Fim --------------------------------
Write-Host "[8/8] Finalizado!"

# (opcional) limpar o version_info depois do build:
# Remove-Item -Force $versionInfoPath -ErrorAction SilentlyContinue
