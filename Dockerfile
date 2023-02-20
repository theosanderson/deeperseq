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

# Add Bioconda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda

# Create a new Conda environment and install necessary packages
RUN conda create -n myenv fastapi uvicorn minimap2 sra-tools fastq-dl

# Activate the environment
RUN echo "conda activate myenv" > ~/.bashrc
ENV PATH /opt/conda/envs/myenv/bin:$PATH

# Copy the FastAPI app into the container
COPY app.py .
COPY ref.fa .
COPY count_lines.py .

# Expose the necessary port
EXPOSE 80

# Start the FastAPI app using uvicorn
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "80"]