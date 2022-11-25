FROM ubuntu:20.04

LABEL maintainer="tmiller@mochsl.org.br"

ENV SIDER_ROOT="/sider"

COPY . ${SIDER_ROOT}

RUN apt-get update \
  && apt-get install -y \
    gcc \
    git \
    make \
    meson \
    ninja-build \
    python3 \
    python3-setuptools \
    python3-wheel \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN cd ${SIDER_ROOT} \
  && meson build \
  && ninja -C build \
  && ninja -C build install

ENTRYPOINT ["sider"]
