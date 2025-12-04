#!/usr/bin/env bash
# =============================================================================
# LncRTPred — Entorno de dependencias oficial (Python 3.8 + libs exactas)
# Crea o repara el entorno 'lncrtpred_env' conforme a la documentación oficial.
# Idempotente: si el entorno ya existe, sólo verifica e instala faltantes.
# =============================================================================
set -euo pipefail

ENV_NAME="lncrtpred_env"
PY_VER="3.8"
REQS=(
  "pandas==1.5.3"
  "numpy==1.20.3"
  "biopython==1.78"
  "scikit-learn==0.24.2"
  "pickleshare==0.7.5"
  "lightgbm==2.3.1"
)

log(){ echo -e "[`date +'%F %T'`] $*"; }

# -------------------- Detecta conda --------------------
if ! command -v conda >/dev/null 2>&1; then
  echo "[ERR] No se encontró conda. Instala Miniconda o Anaconda primero." >&2
  exit 1
fi

# -------------------- Crea entorno si no existe --------------------
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  log "Creando entorno ${ENV_NAME} (Python ${PY_VER})…"
  conda create -y -n "$ENV_NAME" python="$PY_VER"
else
  log "Entorno ${ENV_NAME} ya existe. Verificando paquetes…"
fi

# -------------------- Activar entorno --------------------
# shellcheck disable=SC1091
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# -------------------- Instalar dependencias --------------------
log "Instalando dependencias exactas según documentación…"
for pkg in "${REQS[@]}"; do
  if python -c "import pkg_resources; pkg_resources.require('$pkg')" >/dev/null 2>&1; then
    log "[OK] $pkg ya instalado."
  else
    log "[+] Instalando $pkg"
    conda install -y "$pkg" || pip install "$pkg"
  fi
done

# -------------------- Verificación final --------------------
log "Verificando versiones instaladas:"
python - <<'EOF'
import pandas, numpy, Bio, sklearn, pickleshare, lightgbm, sys
print("Python:", sys.version.split()[0])
print("pandas:", pandas.__version__)
print("numpy:", numpy.__version__)
print("BioPython:", Bio.__version__)
print("scikit-learn:", sklearn.__version__)
print("pickleshare:", pickleshare.__version__)
print("lightgbm:", lightgbm.__version__)
EOF

log "✅ Entorno ${ENV_NAME} listo. Usa:"
echo ""
echo "    conda activate ${ENV_NAME}"
echo "    cd ~/LncRTPred/Human"
echo "    python predict_lncrna_mrna_human.py"
echo ""
