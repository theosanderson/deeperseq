import subprocess
import logging
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler("app.log"),
        logging.StreamHandler()
    ]
)


def raw_to_bam(accession: str):
    logging.info(f"Processing SRA accession {accession}")
    yield f"Processing SRA accession {accession}"

    # Download FASTQ files
    cmd_fastq_dump = f'fastq-dump --split-files {accession}'
    subprocess.run(cmd_fastq_dump, shell=True, check=True)

    logging.info(f"Downloaded FASTQ files for {accession}")
    yield f"Downloaded FASTQ files for {accession}"

    # Align to reference genome using minimap2 and pipe to samtools view
    cmd_minimap2 = f'minimap2 -a ./ref.fa {accession}_1.fastq {accession}_2.fastq | samtools view -bS - > {accession}.bam'
    subprocess.run(cmd_minimap2, shell=True, check=True)

    logging.info(f"Aligned reads for {accession}")
    yield f"Aligned reads for {accession}"


    cmd_sort_bam = f'samtools sort -o {accession}.sorted.bam {accession}.bam'
    subprocess.run(cmd_sort_bam, shell=True, check=True)

    logging.info(f"Sorted BAM file for {accession}")
    yield f"Sorted BAM file for {accession}"

    cmd_index_bam = f'samtools index {accession}.sorted.bam'
    subprocess.run(cmd_index_bam, shell=True, check=True)

    logging.info(f"Indexed BAM file for {accession}")
    yield f"Indexed BAM file for {accession}"

    return f"Finished processing {accession}"

@app.get("/align/{accession}")
async def align(accession: str):
    logging.info(f"Aligning {accession} to reference genome")
    return StreamingResponse(raw_to_bam(accession))
    

@app.get("/dl_align/{accession}")
async def index(accession: str):
    logging.info(f"Returning BAM index file for {accession}")

    # Return contents of BAM index file
    def generate():
        with open(f"{accession}.sorted.bam", "rb") as file:
            yield from file

    return StreamingResponse(generate(), media_type="application/octet-stream")

@app.get("/index/{accession}")
async def index(accession: str):
    logging.info(f"Returning BAM index file for {accession}")

    # Return contents of BAM index file
    def generate():
        with open(f"{accession}.sorted.bam.bai", "rb") as file:
            yield from file

    return StreamingResponse(generate(), media_type="application/octet-stream")