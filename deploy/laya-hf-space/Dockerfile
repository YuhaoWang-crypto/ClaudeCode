# Hugging Face Docker Space: serves convaiinnovations/laya over HTTP (laya-serve).
FROM python:3.11-slim-bookworm

ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN pip install "torch==2.14.0" --index-url https://download.pytorch.org/whl/cpu \
    && pip install "laya[serve]==0.3.11"

# HF Spaces run the container as UID 1000 and route traffic to port 7860.
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    HF_HOME=/home/user/.cache/huggingface \
    PYTHONUNBUFFERED=1 \
    TOKENIZERS_PARALLELISM=false \
    OMP_NUM_THREADS=2 \
    LAYA_HOST=0.0.0.0 \
    LAYA_PORT=7860 \
    LAYA_DEVICE=cpu \
    LAYA_PRELOAD=1 \
    LAYA_THREADS=2
WORKDIR /home/user

EXPOSE 7860
# LAYA_API_KEY comes from the Space's secrets; when set, callers must send
# "Authorization: Bearer <key>".
CMD ["laya-serve"]
