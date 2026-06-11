@echo off
setlocal

set ROOT=%~dp0

echo Checking backend virtual environment...
if not exist "%ROOT%backend\.venv\Scripts\activate.bat" (
    echo Virtual environment not found. Creating .venv...
    cd /d "%ROOT%backend"
    python -m venv .venv
)

echo Installing backend requirements...
call "%ROOT%backend\.venv\Scripts\activate.bat"
python -m pip install --upgrade pip
pip install -r "%ROOT%backend\requirements.txt"

echo Starting backend...
start "Antibody Backend" cmd /k "cd /d %ROOT%backend && call .venv\Scripts\activate.bat && python -muvicorn app.main:app --reload"

echo Installing frontend dependencies...
cd /d "%ROOT%frontend"
call npm install

echo Starting frontend...
start "Antibody Frontend" cmd /k "cd /d %ROOT%frontend && npm run dev"

echo.
echo Backend:  http://localhost:8000/api
echo Frontend: http://localhost:5173
echo.
echo Two terminal windows were opened.
pause
