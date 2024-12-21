FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bedtools \
        python3.12 \
        python3.12-distutils \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python3.12 -m pip install --no-cache-dir --upgrade pip

WORKDIR /app

COPY . /app

RUN pip install --no-cache-dir .


CMD ["bash"]