FROM python:3.9-slim

# Install system dependencies including RDKit prerequisites and Graphviz
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    gfortran \
    make \
    gawk \
    libnetcdf-dev \
    libnetcdff-dev \
    graphviz \
    libxrender1 \
    libxext6 \
    libcairo2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements and install
COPY gecko_web/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY gecko_web ./gecko_web
COPY samples ./samples

# Copy GECKO-A source and compile
COPY docker/gecko_source ./gecko_source
RUN cd gecko_source/OBJ && make clean && make
RUN mkdir -p gecko_source/RUN/OUT

# Configure GECKO-A paths
RUN sed -i "s|dirgecko=.*|dirgecko='/app/gecko_source/'|" gecko_source/INPUT/gecko.nml && \
    sed -i "s|dirout=.*|dirout='/app/gecko_source/RUN/OUT/'|" gecko_source/INPUT/gecko.nml

# Copy BOXMODEL source and compile
COPY docker/boxmodel_source ./boxmodel_source
COPY docker/docker_arch.gnu ./boxmodel_source/SCRIPTS/arch/docker_arch.gnu
RUN cd boxmodel_source && ./build.sh --arch docker_arch.gnu

# Create data directories
RUN mkdir -p /data/output /data/mechanisms

# Expose port
EXPOSE 8000

# Command to run the application
CMD ["uvicorn", "gecko_web.main:app", "--host", "0.0.0.0", "--port", "8000"]
