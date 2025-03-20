ARG VERSION


FROM julia:${VERSION} AS system_stage

ARG SPE11_DIR=/opt/spe11csp

RUN apt update
RUN apt install -y git git-lfs
RUN apt install -y openssh-server ca-certificates vim

FROM system_stage AS spe11csp_stage

RUN git clone https://github.com/jafranc/CSP11_JutulDarcy.jl -b feat/docker ${SPE11_DIR}

WORKDIR ${SPE11_DIR}
RUN git-lfs install
RUN git-lfs pull 
RUN cp -v docker/Project.toml .

