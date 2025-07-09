# bkyanjo
# docker build --build-arg MATLAB_RELEASE=R2025a
#              --build-arg MATLAB_PRODUCT_LIST="MATLAB Deep_Learning_Toolbox Symbolic_Math_Toolbox"
#              --build-arg MATLAB_INSTALL_LOCATION="/opt/matlab/R2025a"
#              --build-arg LICENSE_SERVER=12345@hostname.com
#              -t my_matlab_image_name .


# Use uppercase to specify the release, for example: ARG MATLAB_RELEASE=R2021b
ARG MATLAB_RELEASE=R2025a

ARG MATLAB_RELEASE=R2025a

# Specify the list of products to install into MATLAB.
ARG MATLAB_PRODUCT_LIST="MATLAB"

# Specify MATLAB Install Location.
ARG MATLAB_INSTALL_LOCATION="/opt/matlab/${MATLAB_RELEASE}"

# Specify license server information using the format: port@hostname
ARG LICENSE_SERVER=1711@matlablic.ecs.gatech.edu

FROM mathworks/matlab-deps:${MATLAB_RELEASE}

# Declare build arguments to use at the current build stage.
ARG MATLAB_RELEASE
ARG MATLAB_PRODUCT_LIST
ARG MATLAB_INSTALL_LOCATION
ARG LICENSE_SERVER

# Install mpm dependencies
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install --no-install-recommends --yes \
        wget git bash gcc gfortran g++ make python3-dev libtool cmake autoconf automake autotools-dev \
        curl python3 python3-numpy python3-pip python3-tk python3-venv \
        python3-requests zlib1g-dev libboost-dev sudo gmsh \
        bison flex ninja-build \
        libocct-ocaf-dev libocct-data-exchange-dev \
        swig graphviz \
        libcurl4-openssl-dev libxml2-dev libssl-dev \
        vim docker.io openssh-client build-essential unzip patch ca-certificates libopenblas-serial-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab \
&& echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab \
&& chmod 0440 /etc/sudoers.d/matlab

# Set user and work directory.
USER matlab
WORKDIR /home/matlab

# Pass in $HOME variable to install support packages into the user's HOME folder.
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && sudo HOME=${HOME} ./mpm install \
    --release=${MATLAB_RELEASE} \
    --destination=${MATLAB_INSTALL_LOCATION} \
    --products ${MATLAB_PRODUCT_LIST} \
    || (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) \
    && sudo rm -rf mpm /tmp/mathworks_root.log \
    && sudo ln -s ${MATLAB_INSTALL_LOCATION}/bin/matlab /usr/local/bin/matlab


# Add MATLAB license
ENV MLM_LICENSE_FILE=$LICENSE_SERVER

# Add MATLAB to path
ENV PATH="${MATLAB_INSTALL_LOCATION}/bin:${PATH}"

# Use a more sane locale
ENV LC_ALL=C.UTF-8

# Change the `ubuntu` user to `ISSM`
# and ensure that we do not run as root on self-hosted systems
USER root 
RUN usermod -d /opt/ISSM -m ubuntu && \
    usermod -l ISSM ubuntu && \
    groupmod -n ISSM ubuntu && \
    usermod -aG sudo ISSM && \
    echo "ISSM:docker" | chpasswd && \
    echo "ISSM ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    ldconfig

# Remove /opt/ISSM if it exists and clone ISSM repository as root
RUN rm -rf /opt/ISSM \
    && git clone https://github.com/ISSMteam/ISSM.git /opt/ISSM \
    && chown -R ISSM:ISSM /opt/ISSM

# Set ISSM user’s home directory to /opt/ISSM
RUN usermod -d /opt/ISSM -m ISSM

# Switch to ISSM user
USER ISSM
WORKDIR /opt/ISSM

# set ISSM variable
ENV ISSM_DIR=/opt/ISSM
RUN ["/bin/bash", "-c", "source $ISSM_DIR/etc/environment.sh"]

# Install external packages
WORKDIR $ISSM_DIR/externalpackages/autotools
RUN ./install-linux.sh

WORKDIR $ISSM_DIR/externalpackages/cmake
RUN ./install.sh 

# WORKDIR $ISSM_DIR/externalpackages/petsc
# RUN bash -c 'VER="3.22.3"; \
#     PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src"; \
#     PREFIX="${ISSM_DIR}/externalpackages/petsc/install"; \
#     LDFLAGS="${LDFLAGS:-}"; \
#     "${ISSM_DIR}/scripts/DownloadExternalPackage.sh" "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz" "petsc-${VER}.tar.gz"; \
#     tar -zxvf "petsc-${VER}.tar.gz"; \
#     rm -rf "${PREFIX}" "${PETSC_DIR}"; \
#     mkdir -p "${PETSC_DIR}"; \
#     mv "petsc-${VER}"/* "${PETSC_DIR}/"; \
#     rm -rf "petsc-${VER}"; \
#     cd "${PETSC_DIR}"; \
#     ./configure \
#         --prefix="${PREFIX}" \
#         --PETSC_DIR="${PETSC_DIR}" \
#         --LDFLAGS="${LDFLAGS}" \
#         --COPTFLAGS="-O3 -march=native -mtune=native" \
#         --CXXOPTFLAGS="-O3 -march=native -mtune=native" \
#         --FOPTFLAGS="-O3 -march=native -mtune=native" \
#         --with-c2html=0 \
#         --with-debugging=0 \
#         --with-valgrind=0 \
#         --with-x=0 \
#         --with-ssl=0 \
#         --with-pic=1 \
#         --with-zlib \
#         --with-fortran-bindings=0 \
#         --with-make-np=4 \
#         --with-shared-libraries=1 \
#         --download-fblaslapack=1 \
#         --download-metis=1 \
#         --download-mpich=1 \
#         --download-mpich-device=ch3:sock \
#         --download-mumps=1 \
#         --download-parmetis=1 \
#         --download-scalapack=1 \
#         --download-hdf5=1 \
#         --download-hwloc=1 \
#         --download-zlib=1 \
#         PETSC_ARCH=install; \
#     mv install/include/petscconf.h install/include/old_petscconf.nope;'

# # build default ISSM PETSC
# RUN bash -c 'export PACKAGES="/opt/ISSM/externalpackages/petsc/install"; \
#     PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src"; \
#     PREFIX="${ISSM_DIR}/externalpackages/petsc/install"; \
#     cd ${PETSC_DIR}; \
#     ./configure \
#         --prefix="${PREFIX}" \
#         --PETSC_DIR="${PETSC_DIR}" \
#         --LDFLAGS="${LDFLAGS:-}" \
#         --with-debugging=0 \
#         --with-valgrind=0 \
#         --with-x=0 \
#         --with-ssl=0 \
#         --with-pic=1 \
#         --with-zlib-dir="${PACKAGES}" \
#         --with-fblaslapack-dir="${PACKAGES}" \
#         --with-metis-dir="${PACKAGES}" \
#         --with-mpi-dir="${PACKAGES}" \
#         --with-mumps-dir="${PACKAGES}" \
#         --with-parmetis-dir="${PACKAGES}" \
#         --with-scalapack-dir="${PACKAGES}"'

# RUN bash -c '\
#     set -e; \
#     PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src"; \
#     cd "${PETSC_DIR}"; \
#     make -j4 all && \
#     make install'

WORKDIR $ISSM_DIR/externalpackages/petsc
ENV PETSC_DIR=${ISSM_DIR}/externalpackages/petsc/src
ENV PREFIX=${ISSM_DIR}/externalpackages/petsc/install
ENV VER=3.22.0

# Download and unpack PETSc
RUN wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz \
    && tar -zxvf petsc-${VER}.tar.gz \
    && rm -rf ${PREFIX} ${PETSC_DIR} \
    && mkdir -p ${PETSC_DIR} \
    && mv petsc-${VER}/* ${PETSC_DIR} \
    && rm -rf petsc-${VER} petsc-${VER}.tar.gz

# Configure, compile, and install PETSc
WORKDIR ${PETSC_DIR}
RUN ./configure \
    --prefix="${PREFIX}" \
    --PETSC_DIR="${PETSC_DIR}" \
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-pic=1 \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-mpich=1 \
    --download-mpich-device=ch3:sock \
    --download-mumps=1 \
    --download-parmetis=1 \
    --download-scalapack=1 \
    --download-zlib=1 \
    && make \
    && make install

WORKDIR $ISSM_DIR/externalpackages/triangle
RUN ./install-linux.sh

WORKDIR $ISSM_DIR/externalpackages/m1qn3
RUN ./install-linux.sh

# souce ISSM environment
RUN ["/bin/bash", "-c", "source $ISSM_DIR/etc/environment.sh"]

WORKDIR ${ISSM_DIR}
RUN autoreconf -ivf

# configure and build ISSM
RUN bash -c './configure \
   --prefix=${ISSM_DIR} \
   --with-numthreads=2 \
   --with-matlab-dir="/opt/matlab/R2025a" \
   --with-fortran-lib="-L/usr/lib/x86_64-linux-gnu -lgfortran" \
   --with-mpi-include="${ISSM_DIR}/externalpackages/petsc/install/include" \
   --with-mpi-libflags="-L${ISSM_DIR}/externalpackages/petsc/install/lib -lmpi -lmpicxx -lmpifort" \
   --with-metis-dir="${ISSM_DIR}/externalpackages/petsc/install" \
   --with-parmetis-dir="${ISSM_DIR}/externalpackages/petsc/install" \
   --with-blas-lapack-dir="${ISSM_DIR}/externalpackages/petsc/install" \
   --with-scalapack-dir="${ISSM_DIR}/externalpackages/petsc/install" \
   --with-mumps-dir="${ISSM_DIR}/externalpackages/petsc/install" \
   --with-petsc-dir="${ISSM_DIR}/externalpackages/petsc/install" \
   --with-triangle-dir="${ISSM_DIR}/externalpackages/triangle/install" \
   --with-m1qn3-dir="${ISSM_DIR}/externalpackages/m1qn3/install"'

RUN bash -c 'make -j4 && make install;'

# Clean up build artifacts
RUN rm -rf /opt/ISSM/externalpackages/petsc/**/externalpackages \
    && rm -rf /opt/ISSM/externalpackages/petsc/src/docs \
    && rm -rf /opt/ISSM/externalpackages/petsc/src/**/tutorials/output/* \
    && rm -rf /opt/ISSM/externalpackages/petsc/src/**/tests/output/* 

# SET ENVIRONMENT VARIABLES
ENV ISSM_DIR=/opt/ISSM \
    MPICH_DIR=/opt/ISSM/externalpackages/petsc/install/bin \
    HDF5_DIR=/opt/ISSM/externalpackages/petsc/install \
    HDF5_MPI=ON \
    PATH=$PATH:/opt/ISSM