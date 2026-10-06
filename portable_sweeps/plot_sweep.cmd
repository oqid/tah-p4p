@echo off
setlocal
if not exist "%~dp0.venv\Scripts\python.exe" (
    echo Run setup_windows.cmd first.
    pause
    exit /b 1
)
"%~dp0.venv\Scripts\python.exe" "%~dp0plot_sweep.py" %*
set "sweep_exit=%errorlevel%"
pause
exit /b %sweep_exit%
