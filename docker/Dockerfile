FROM ubuntu:20.04 AS builder

RUN set -e; \
    \
    apt-get update; \
    apt-get install -y \
      gcc \
      git \
      make \
      meson \
      ninja-build \
      python3 \
      python3-setuptools \
      python3-wheel

COPY . /sider

WORKDIR /sider

RUN set -e; \
    \
    meson setup \
      --buildtype=release \
      --default-library=static \
      --wrap-mode=forcefallback \
      build; \
    ninja -C build

FROM ubuntu:20.04

LABEL maintainer="tmiller@mochsl.org.br"

COPY --from=builder /sider/build/src/sider /usr/local/bin/

ENTRYPOINT ["sider"]
