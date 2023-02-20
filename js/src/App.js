import React, { useState,useRef , useEffect} from 'react';
import axios from 'axios';
import ClipLoader from 'react-spinners/ClipLoader';
import igv from "igv"
import './App.css';
const backend = ""



function Alignment() {
  const [accession, setAccession] = useState('ERR8254282');

  
  const [taskID, setTaskID] = useState('');
  const [logs, setLogs] = useState([]);
  const [lines, setLines] = useState([]);
  const [status, setStatus] = useState('');
  const [error, setError] = useState('');
  const [bamURL, setBamURL] = useState('');
  const [baiURL, setBaiURL] = useState('');
  const igvDiv = useRef(null);


  useEffect(() => {
    if (status=="complete") {
      const bamURL = `${backend}/${accession}.sorted.bam`;
      const baiURL = `${backend}/${accession}.sorted.bam.bai`;
      console.log("creating IGV");
      igv.createBrowser(igvDiv.current, {
        //ref is /ref.fa
        reference: {
          genome: 'ASM985889v3',
        
        },
        locus: 'NC_045512.2:1-29903',
        tracks: [
          {
            type: 'alignment',
            format: 'bam',
            url: bamURL,
            indexURL: baiURL,
            name: 'Aligned Reads',
            displayMode: 'SQUISHED',

          },
        ],
      });
    }
  }, [status]);







  const handleAccessionChange = (event) => {
    setAccession(event.target.value);
  };

  const doAlign = async (accession) => {
    try {
      // check accession starts with SRR or ERR
      accession = accession.trim();
      if (!accession.startsWith('SRR') && !accession.startsWith('ERR')) {
        setError('Accession must start with SRR or ERR');
        return;
      }
      const response = await axios.post(`${backend}/align/${accession}`);
      const taskID = response.data.task_id;
      setTaskID(taskID);
      setStatus('processing');
      setError('');
      setLogs([]);
      pollTask(taskID);
    } catch (error) {
      console.error(error);
      setError('Error starting alignment task');
    }
  };

  const handleAlignmentSubmit = async (event) => {
    event.preventDefault();
    doAlign(accession);
   
  };

  // if the query string has an accession, start the alignment
  useEffect(() => {
    const urlParams = new URLSearchParams(window.location.search);
    const acc = urlParams.get('acc');
    if (acc) {
      setAccession(acc);
      doAlign(acc);
    }
  }, []);


  const setAcc = (acc) => {
    setBaiURL(`${backend}/index/${acc}`);
    setBamURL(`${backend}/dl_align/${acc}`);
    setStatus('complete');
  };

  const pollTask = async (taskID) => {
    try {
      const response = await axios.get(`${backend}/poll/${taskID}`);
      const taskStatus = response.data.status;
      if (taskStatus === 'processing') {
        setLogs(response.data.log);
        setLines(response.data.lines);
        setTimeout(() => pollTask(taskID), 5000);
      } else if (taskStatus === 'complete') {
        setLogs(response.data.log);
        setStatus('complete');
        setBamURL(`${backend}/dl_align/${accession}`);
        setBaiURL(`${backend}/index/${accession}`);
      }
    } catch (error) {
      console.error(error);
      setError('Error polling alignment task');
      setTimeout(() => pollTask(taskID), 5000);
    }
  };

  return (
    // overflow-y-auto
    <div className="bg-gray-100 min-h-screen overflow-y-auto pb-3">
      <div className="pb-3">
    <div className="flex flex-col items-center space-y-4">
      
      <h1 className="text-2xl font-bold mt-2">DeepSeqer</h1>
      <p className="text-gray-500">A tool for viewing SARS-CoV-2 deep sequencing data</p>
      {status !== 'complete' && status!== 'processing' && (
        <>
         <p
        className='text-center'>Enter a SRR/ERR accession. This will be mapped to the Hu-1 reference and then displayed<br />(It will be downsampled to 100,000 reads).</p>
      <form className="flex flex-col space-y-2" onSubmit={handleAlignmentSubmit}>
       

        <label className="font-semibold">
          SRA/ENA Accession ID:
          <input
            type="text"
            value={accession}
            onChange={handleAccessionChange}
            className="border border-gray-300 rounded-md px-2 py-1 focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </label>
        <button type="submit" className="bg-blue-500 text-white px-4 py-2 rounded-md">
          Start Alignment
        </button>
      </form>
      </>
      )}
      {error && <p className="text-red-500">{error}</p>}
      {status === 'processing' && (
        <div className="flex flex-col items-center space-y-4">
          <ClipLoader color="#000" loading={true} size={150} />
          <h2 className="text-xl font-bold">Alignment Processing</h2>
          <p>Lines: {lines}</p>
          <h3 className="text-lg font-bold">Logs</h3>
          <ul className="divide-y divide-gray-300 w-full max-w-md">
            {logs.map((log, index) => log && <li key={index}><pre>{log}</pre></li>)}
          </ul>
        </div>
      )}

      {status === 'complete' && (
        <div className="flex flex-col items-center space-y-4 w-full">
          
          <h3 className="text-lg font-bold">{accession} reads</h3>
          <div
            id="igv-div"
            className="bg-white border border-gray-300 rounded-md  mb-3"
            style={{ width: 'calc(100vw -5 em)'}}
            ref={igvDiv}
          />
          <p className="text-center text-gray-500 text-sm">Reads may be downsampled.</p>
         
        </div>
      )}
    </div>
    </div>
    </div>
  );
};

export default Alignment;