# create_env.ps1
# Uso:
#   . .\create_env.ps1                 # dot-source (mantém a venv ativada na sessão atual)
#   # ou
#   powershell -ExecutionPolicy Bypass -File .\create_env.ps1

[CmdletBinding()]
param(
  [string]$VenvPath = '.venv',
  [string]$PythonExe = ''  # Se vazio, detecta automaticamente
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Write-Step([string]$msg) { Write-Host ("[SETUP] {0}" -f $msg) -ForegroundColor Cyan }
function Write-Ok  ([string]$msg) { Write-Host ("[ OK ] {0}" -f $msg)   -ForegroundColor Green }
function Write-Warn([string]$msg) { Write-Host ("[WARN] {0}" -f $msg)   -ForegroundColor Yellow }
function Write-Err ([string]$msg) { Write-Host ("[FAIL] {0}" -f $msg)   -ForegroundColor Red }

# --- Detecta Python ---
if (-not $PythonExe) {
  $candidates = @('py -3.13','py -3.12','py -3.11','python','python3')
  foreach ($c in $candidates) {
    try {
      $v = & $c -c 'import sys;print(sys.version)' 2>$null
      if ($LASTEXITCODE -eq 0 -and $v) { $PythonExe = $c; break }
    } catch {}
  }
  if (-not $PythonExe) {
    Write-Err 'Python não encontrado. Instale Python 3.x (recomendo 3.13) e tente novamente.'
    exit 1
  }
}
Write-Step 'Usando Python: $PythonExe'

# --- Cria a venv (se não existir) ---
if (-not (Test-Path $VenvPath)) {
  Write-Step 'Criando ambiente virtual em '$VenvPath'...'
  & $PythonExe -m venv $VenvPath
  if ($LASTEXITCODE -ne 0) { Write-Err 'Falha ao criar a venv.'; exit 1 }
  Write-Ok 'Venv criada.'
} else {
  Write-Step 'Venv já existe em '$VenvPath'.'
}

# --- Ativa a venv na sessão atual (se o script for dot-sourced) ---
$activate = Join-Path $VenvPath 'Scripts\Activate.ps1'
if (Test-Path $activate) {
  try {
    . $activate
    Write-Ok 'Venv ativada.'
  } catch {
    Write-Warn 'Não consegui ativar automaticamente. Ative manualmente: `. $activate`'
  }
} else {
  Write-Err 'Arquivo de ativação não encontrado: $activate'
  exit 1
}

# --- Garante pip e atualiza ferramentas base ---
Write-Step 'Atualizando pip / setuptools / wheel...'
python -m ensurepip --upgrade | Out-Null
python -m pip install --upgrade pip setuptools wheel
if ($LASTEXITCODE -ne 0) { Write-Err 'Falha ao atualizar pip/setuptools/wheel.'; exit 1 }
Write-Ok 'Ferramentas atualizadas.'

# --- Escolhe arquivo de requisitos ---
$reqTxt = 'requirements.txt'
$reqIn  = 'requirements.in'
$reqToUse = $null

if (Test-Path $reqTxt) {
  $reqToUse = $reqTxt
} elseif (Test-Path $reqIn) {
  $reqToUse = $reqIn
} else {
  Write-Warn 'Nenhum requirements.txt ou requirements.in encontrado. Nada para instalar.'
}

# --- Instala dependências ---
if ($reqToUse) {
  Write-Step 'Instalando dependências de '$reqToUse'...'
  python -m pip install -r $reqToUse
  if ($LASTEXITCODE -ne 0) { Write-Err 'Falha ao instalar dependências de '$reqToUse'.'; exit 1 }
  Write-Ok 'Dependências instaladas.'
}

# --- Smoke test opcional dos imports mais pesados (sem heredoc) ---
try {
  Write-Step 'Rodando smoke test rápido de imports (matplotlib, numpy, pandas, scipy)...'
  $smoke = @'
import importlib
for m in ('numpy','pandas','matplotlib','scipy'):
    importlib.import_module(m)
print('Smoke test OK.')
'@
  $tmp = New-TemporaryFile
  Set-Content -Path $tmp -Value $smoke -Encoding UTF8
  python $tmp
  Remove-Item $tmp -Force
  if ($LASTEXITCODE -eq 0) { Write-Ok 'Imports básicos OK.' }
} catch {
  Write-Warn 'Smoke test falhou (pode ainda assim estar tudo certo para seu projeto).'
}

Write-Host ''
Write-Ok 'Ambiente pronto!'
Write-Host 'Para ativar em novas sessões:'
Write-Host '.\.venv\Scripts\Activate.ps1'
