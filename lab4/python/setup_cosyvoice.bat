@echo off
setlocal enableextensions enabledelayedexpansion

set "PYDIR=%~dp0"
if "%PYDIR:~-1%"=="\" set "PYDIR=%PYDIR:~0,-1%"
set "VENV=%PYDIR%\.cosyvenv"

echo ===== Lab 4 setup: CosyVoice + kNN-VC =====
echo Workdir: %PYDIR%
echo Venv:    %VENV%
echo.

py -3.11 --version >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python 3.11 not found via py launcher. Install it from python.org and retry.
    exit /b 1
)

if not exist "%VENV%" (
    echo [1/4] Creating venv...
    py -3.11 -m venv "%VENV%"
    if errorlevel 1 (
        echo ERROR: venv creation failed.
        exit /b 1
    )
) else (
    echo [1/4] Venv already exists, reusing.
)

call "%VENV%\Scripts\activate.bat"
if errorlevel 1 (
    echo ERROR: failed to activate venv at %VENV%
    exit /b 1
)

echo [2/4] Updating pip toolchain...
python -m pip install --upgrade pip wheel setuptools

echo [3/4] Running setup_cosyvoice.py (torch CUDA + clone + deps + models)...
python "%PYDIR%\setup_cosyvoice.py" %*
set "EC=%errorlevel%"
if not "%EC%"=="0" (
    echo ERROR: setup_cosyvoice.py exited with code %EC%
    exit /b %EC%
)

echo.
echo [4/4] Done. The venv is at:
echo     %VENV%\Scripts\python.exe
echo Lab 4 panel will use it automatically.
echo.
endlocal
