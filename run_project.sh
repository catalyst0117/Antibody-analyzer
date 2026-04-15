#!/bin/bash

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME="cs180_env"

osascript <<EOF
tell application "Terminal"
    activate

    do script "cd \"$ROOT_DIR/backend\"; source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate $ENV_NAME; uvicorn app.main:app --reload"

    do script "cd \"$ROOT_DIR/frontend\"; npm run dev"
end tell
EOF

echo "Backend:  http://localhost:8000/api"
echo "Frontend: http://localhost:5173"