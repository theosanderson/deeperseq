import asyncio
import logging
import subprocess
import uuid

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

# Store running tasks and logs
tasks = {}
logs = {}

async def raw_to_bam(accession: str, task_id: str):
    logging.info(f"Processing SRA accession {accession}")
    logs[task_id].append(f"Processing SRA accession {accession}")

    # Download FASTQ files
    cmd_fastq_dump = f'fasterq-dump --split-files {accession}'
    proc = await asyncio.create_subprocess_shell(cmd_fastq_dump)
    await proc.wait()
    logging.info(f"Downloaded FASTQ files for {accession}")
    logs[task_id].append(f"Downloaded FASTQ files for {accession}")

    # Align to reference genome using minimap2 and pipe to samtools view
    cmd_minimap2 = f'minimap2 -a ./ref.fa {accession}_1.fastq {accession}_2.fastq | samtools view -bS - > {accession}.bam'
    proc = await asyncio.create_subprocess_shell(cmd_minimap2)
    await proc.wait()
    logging.info(f"Aligned reads for {accession}")
    logs[task_id].append(f"Aligned reads for {accession}")

    cmd_sort_bam = f'samtools sort -o {accession}.sorted.bam {accession}.bam'
    proc = await asyncio.create_subprocess_shell(cmd_minimap2)
    await proc.wait()
    logging.info(f"Sorted BAM file for {accession}")
    logs[task_id].append(f"Sorted BAM file for {accession}")

    cmd_index_bam = f'samtools index {accession}.sorted.bam'
    proc = await asyncio.create_subprocess_shell(cmd_minimap2)
    await proc.wait()
    logging.info(f"Indexed BAM file for {accession}")
    logs[task_id].append(f"Indexed BAM file for {accession}")

    logs[task_id].append(f"Finished processing {accession}")

async def start_task(accession: str):
    task_id = str(uuid.uuid4())
    logs[task_id] = []
    task = asyncio.create_task(raw_to_bam(accession, task_id))
    tasks[task_id] = task
    return task_id

@app.post("/align/{accession}")
async def align(accession: str):
    logging.info(f"Aligning {accession} to reference genome")
    task_id = await start_task(accession)
    return {"task_id": task_id}

@app.get("/poll/{task_id}")
async def poll(task_id: str):
    if task_id not in tasks:
        return {"error": "Invalid task ID"}

    if not tasks[task_id].done():
        return {"status": "processing", "log": logs[task_id]}

    if task_id in logs:
        return {"status": "complete", "log": logs[task_id]}
    else:
        return {"status": "complete"}
    


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