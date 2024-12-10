# Start with a base image that includes Mambaforge
FROM quay.io/condaforge/mambaforge

# Set environment variables to configure Conda/Mamba behavior
ENV PATH=/opt/conda/bin:$PATH \
    CONDA_AUTO_UPDATE_CONDA=false

COPY env.yaml /tmp/env.yaml
COPY env2.yaml /tmp/env2.yaml

# Create the environments using mamba for speed
RUN mamba env create -f /tmp/env.yaml -y && \
    mamba env create -f /tmp/env2.yaml -y && \
    mamba clean --all --yes

# Add a script to make it easy to activate either environment
RUN echo 'source /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc && \
    echo 'mamba activate base' >> ~/.bashrc && \
    echo 'function activate_env() { source /opt/conda/etc/profile.d/conda.sh && conda activate $1; }' >> ~/.bashrc

# Set the working directory
WORKDIR /workspace

# Define the default command for the container
CMD ["bash"]

