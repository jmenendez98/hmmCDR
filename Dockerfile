FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y bedtools && \
    apt-get install -y python==3.12 pip

RUN python -m pip install --upgrade pip && \
    pip install .

CMD ["bash"]