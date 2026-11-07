FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ cmake make ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cloth
COPY CMakeLists.txt CMakePresets.json ./
COPY cmake/ cmake/
COPY third_party/ third_party/
COPY include/ include/
COPY src/ src/
COPY apps/ apps/
COPY tests/ tests/
RUN cmake --preset cpu-release \
    && cmake --build --preset cpu-release --parallel 2

CMD ["ctest", "--preset", "cpu-release"]
