@echo off
rem One-off Windows Task Scheduler wrapper for the Stage 4 permutation run.
rem Runs (or resumes) the run, logs status, then deletes its own scheduled task.
cd /d "%~dp0.."
echo %date% %time% task start >> results\logs\perm_task.log
".venv\Scripts\python.exe" -I -u src\permutation.py run --b-full 100 --b-sub 50 --n-cpus 8 >> results\logs\perm_run.out 2>&1
echo %date% %time% runner exited with code %errorlevel% >> results\logs\perm_task.log
".venv\Scripts\python.exe" -I src\permutation.py status >> results\logs\perm_task.log 2>&1
schtasks /Delete /TN "TCGA_BRCA_permutation" /F >> results\logs\perm_task.log 2>&1
