# Latch base image (provides the SDK runtime). Pin to a current tag when registering.
FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

# --- system deps: Java (FlashFry) ---
RUN apt-get update -y && \
    apt-get install -y default-jre-headless wget && \
    rm -rf /var/lib/apt/lists/*

# --- FlashFry (single JAR; builds its own genome index, computes Doench2016 + CFD) ---
RUN wget -q https://github.com/mckenna-lab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar \
        -O /opt/FlashFry.jar

# --- python deps for the pipeline (stdlib only, but rs3 optional for Rule Set 3) ---
RUN python3 -m pip install --no-cache-dir rs3 || true

# --- pipeline code ---
COPY scripts /root/scripts

# Latch SDK metadata expectations
COPY wf /root/wf
COPY version /root/version

WORKDIR /root
