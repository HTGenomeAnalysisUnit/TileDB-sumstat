```dockerfile
# Based on your requested base
FROM cgr.dev/chainguard/wolfi-base:latest

# Where micromamba will install environments
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV CONDA_LOCK=/tmp/conda-linux-64.lock

# Fetch micromamba (static tarball) and install it to /usr/local/bin
# We use ADD to fetch remote URL at build time (Docker supports it)
ADD https://micromamba.snakepit.net/api/micromamba/linux-64/latest /tmp/micromamba.tar.bz2

RUN mkdir -p $MAMBA_ROOT_PREFIX \
    && tar -xjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 \
    && chmod +x /usr/local/bin/micromamba \
    && micromamba --version

# copy the lock file into the image (make sure this file is present next to Dockerfile)
COPY conda-linux-64.lock $CONDA_LOCK

# Create the environment using the lockfile.
# We create the env under $MAMBA_ROOT_PREFIX/envs/scqtl
RUN micromamba create -y -p $MAMBA_ROOT_PREFIX/envs/scqtl --file $CONDA_LOCK \
    && micromamba clean --all -y

# Put the env's bin first in PATH for runtime
ENV PATH=$MAMBA_ROOT_PREFIX/envs/scqtl/bin:$PATH

# Copy repository source into /src and run make install inside the environment
COPY . /src
WORKDIR /src

# Run your make install inside the created environment.
# Use micromamba run so we don't need to "activate" a shell.
RUN micromamba run -p $MAMBA_ROOT_PREFIX/envs/scqtl -- make install

# Default command — starts a shell with the environment available via micromamba run
CMD ["micromamba", "run", "-p", "/opt/conda/envs/scqtl", "--", "sh", "-c", "echo 'Container ready. Use micromamba run -p /opt/conda/envs/scqtl -- <command>' && sleep infinity"]
```