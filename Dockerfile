# Build image
FROM node:14-alpine AS build

WORKDIR /app

# Copy the frontend source code and install dependencies
COPY js/package.json js/yarn.lock ./
RUN yarn install --frozen-lockfile

COPY js ./
RUN yarn build

# Final image
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

RUN conda create -n myenv fastapi uvicorn minimap2 sra-tools fastq-dl==2.0.4 seqtk

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

# Copy the built frontend from the build image
COPY --from=build /app/build /var/www/html

# Expose the necessary port
EXPOSE 80
EXPOSE 81

# Start FastAPI and nginx
CMD ["bash", "servers.sh"]
