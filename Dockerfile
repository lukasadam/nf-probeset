FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libgomp1 \
    libglib2.0-0 \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/nf-probeset

COPY requirements.txt ./requirements.txt
RUN python -m pip install --upgrade pip \
    && python -m pip install -r requirements.txt

COPY bin ./bin

ENTRYPOINT ["python", "/opt/nf-probeset/bin/spapros.py"]
