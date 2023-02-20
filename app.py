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

import glob, os



async def raw_to_bam(accession: str, task_id: str):

    def do_log(s):
        logs[task_id].append(s)
        logging.info(s)

    do_log(f"Processing SRA accession {accession}")

    # Download FASTQ files
    #cmd_fastq_dump = f'fasterq-dump --split-files {accession}'
    # using Aspera from ENA
    cmd_fastq_dump = f'fastq-dl -a {accession}'
    proc = await asyncio.create_subprocess_shell(cmd_fastq_dump)
    await proc.wait()
    
    # check what files were downloaded
    files = glob.glob(f"{accession}*.fastq.gz")
    do_log(f"Downloaded {len(files)} files: for {accession}")
    do_log(f"Files: {files}")
    type = None
    for file in files:
        if "1.fastq.gz" in file:
            type = "paired"
            break
        else:
            type = "single"

    do_log(f"Type: {type}")

    # find how many cores are available
    cmd_nproc = f'nproc'
    proc = await asyncio.create_subprocess_shell(cmd_nproc, stdout=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate()
    nproc = int(stdout.decode('utf-8').strip())
    do_log(f"Found {nproc} cores")

    


    if type == "paired":
        cmd_minimap2 = f'minimap2 -a ./ref.fa {accession}_1.fastq.gz {accession}_2.fastq.gz -t {nproc} | python count_lines.py {task_id}.lines |  samtools view -bS - > {accession}.bam'
    else:   
        cmd_minimap2 = f'minimap2 -a ./ref.fa {accession}.fastq.gz | python count_lines.py {task_id}.lines | samtools view -bS -t {nproc} - > {accession}.bam'
    
    proc = await asyncio.create_subprocess_shell(cmd_minimap2)
    await proc.wait()
    do_log(f"Aligned {accession} to reference genome")

    cmd_sort_bam = f'samtools sort -o {accession}.sorted.bam {accession}.bam'
    proc = await asyncio.create_subprocess_shell(cmd_minimap2)
    await proc.wait()
    do_log(f"Sorted BAM file for {accession}")

    cmd_index_bam = f'samtools index {accession}.sorted.bam'
    proc = await asyncio.create_subprocess_shell(cmd_minimap2)
    await proc.wait()
    do_log(f"Indexed BAM file for {accession}")

    do_log(f"Finished processing {accession}")

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
        # read task_id.lines to get current line count
        try:
            with open(f"{task_id}.lines", "r") as f:
                lines = int(f.read())
        except:
            lines = 0
        
        return {"status": "processing", "log": logs[task_id], "lines": lines}

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