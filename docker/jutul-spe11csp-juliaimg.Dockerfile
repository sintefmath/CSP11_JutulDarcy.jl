ARG VERSION


FROM julia:${VERSION} AS system_stage

ARG SPE11_DIR=/opt/spe11csp

RUN apt update
RUN apt install -y git git-lfs
RUN apt install -y openssh-server ca-certificates vim

FROM system_stage AS spe11csp_stage

#ENV SPE_CASE=a
#ENV INPUT_PATH=/opt/spe11csp/examples/
#ENV CASEFILE=csp11${SPE_CASE}

#RUN apt install -y python3 python3-venv python3-pip git
RUN git clone https://github.com/jafranc/CSP11_JutulDarcy.jl -b feat/docker ${SPE11_DIR}

WORKDIR ${SPE11_DIR}
RUN git-lfs install
RUN git-lfs pull 

#RUN julia --project=${SPE11_DIR} \
#	-e "import Pkg; Pkg.develop(path=\"${SPE11_DIR}/CSP11\");Pkg.instantiate();"

#CMD ["julia","--project=\".\"" ,"-e", "\"import Pkg;Pkg.instantiate();Pkg.add([\"Jutul\",\"JutulDarcy\",\"HYPRE\"]);Pkg.develop(path=./CSP11);include(\"docker_run_mrst_grid_spe11.jl\");\""]
#CMD ["sh", "-c", "/usr/local/bin/pyopmcsp11 -i /opt/spe11csp/examples/csp11${SPE_CASE}.txt -o output_csp11${SPE_CASE}"]
