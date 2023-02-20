FROM continuumio/miniconda3:latest

# Install necessary system packages, including aspera
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    samtools \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get -y update && apt-get -y install nginx

RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

# Add Bioconda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda

RUN conda create -n myenv fastapi uvicorn minimap2 sra-tools fastq-dl seqtk

# Activate the environment
RUN echo "conda activate myenv" > ~/.bashrc
ENV PATH /opt/conda/envs/myenv/bin:$PATH

# Copy the FastAPI app into the container
COPY app.py .
COPY ref.fa .
COPY count_lines.py .
COPY servers.sh .
# copy from nginx_config to /etc/nginx/sites-enabled/default
COPY nginx_conf /etc/nginx/sites-enabled/default

# COPY js/build directory into the container
COPY js/build /var/www/html

# Expose the necessary port
EXPOSE 80
EXPOSE 81


# Start FastAPI and nginx
CMD ["bash", "servers.sh"]