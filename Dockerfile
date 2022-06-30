FROM ubuntu:focal

WORKDIR /app

RUN apt update && apt install -y \
    python3.8 \
    python3-pip

RUN python3 -m pip install numpy scipy matplotlib pyamg networkx==2.5

ENTRYPOINT ["python3"]
