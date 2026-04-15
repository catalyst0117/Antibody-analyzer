@echo off
setlocal

set ROOT=%~dp0

echo Starting backend...
start "Antibody Backend" cmd /k "cd /d %ROOT%backend && call .venv\Scripts\activate.bat && uvicorn app.main:app --reload"

echo Starting frontend...
start "Antibody Frontend" cmd /k "cd /d %ROOT%frontend && npm run dev"

echo.
echo Backend:  http://localhost:8000/api
echo Frontend: http://localhost:5173
echo.
echo Two terminal windows were opened.
pause

@REM #!/usr/bin/env bash
@REM set -e

@REM read -p "Enter path to input txt file: " INPUT_TXT
@REM read -p "Enter disease name (press Enter for broad disease search): " DISEASE

@REM if [ -z "$DISEASE" ]; then
@REM     python gene_report.py "$INPUT_TXT"
@REM else
@REM     python gene_report.py "$INPUT_TXT" --disease "$DISEASE"
@REM fi