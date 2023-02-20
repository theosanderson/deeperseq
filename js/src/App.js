import React, { useState,useRef , useEffect} from 'react';
import axios from 'axios';
import igv from "igv"
const backend = "http://localhost"

function Alignment() {
  const [accession, setAccession] = useState('');
  const [taskID, setTaskID] = useState('');
  const [logs, setLogs] = useState([]);
  const [lines, setLines] = useState([]);
  const [status, setStatus] = useState('');
  const [error, setError] = useState('');
  const [bamURL, setBamURL] = useState('');
  const [baiURL, setBaiURL] = useState('');
  const igvDiv = useRef(null);
  useEffect(() => {
    if (bamURL && baiURL) {
      igv.createBrowser(igvDiv.current, {
        //ref is /ref.fa
        reference: {
          id: 'NC_045512.2',
          fastaURL: '/ref.fa',
        },
        locus: 'NC_045512.2:1-29903',
        tracks: [
          {
            type: 'alignment',
            format: 'bam',
            url: bamURL,
            indexURL: baiURL,
            name: 'Aligned Reads',
          },
        ],
      });
    }
  }, [bamURL, baiURL]);







  const handleAccessionChange = (event) => {
    setAccession(event.target.value);
  };

  const handleAlignmentSubmit = async (event) => {
    event.preventDefault();
    try {
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
    <div>
      <h1>Alignment</h1>
      <form onSubmit={handleAlignmentSubmit}>
        <label>
          SRA Accession ID:
          <input type="text" value={accession} onChange={handleAccessionChange} />
        </label>
        <button type="submit">Start Alignment</button>
      </form>
      {error && <p>{error}</p>}
      {status === 'processing' && 

        <div>
          <h2>Alignment Processing</h2>
          <p>Lines: {lines}</p>
          <h3>Logs</h3>
          <ul>
            {logs.map(
              (log, index) =>
                log && (
                  <li key={index}>
                    <pre>{log}</pre>
                  </li>
                )
            )}
          </ul>
        </div>
      }

      {status === 'complete' && (
        <div>
          <h2>Alignment Complete</h2>
          <h3>Aligned Reads</h3>
          <div id="igv-div" style={{ width: '100%', height: '1000px' }} ref={igvDiv} />
          
          
          
          <h3>Logs</h3>
          <ul>
            {logs.map(
              (log, index) =>
                log && (
                  <li key={index}>
                    <pre>{log}</pre>
                  </li>
                )
            )}

          </ul>
        </div>
      )}
    </div>
  );
}


export default Alignment;