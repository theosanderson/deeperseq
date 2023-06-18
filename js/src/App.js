import React, { useState,useRef , useEffect} from 'react';
import axios from 'axios';
import ClipLoader from 'react-spinners/ClipLoader';
import igv from "igv"
import './App.css';


const referenceOptions = [
   {id: "custom", label: "Custom", downsampleTo: 100000},
   {id: "-----", label: "-----", downsampleTo: 100000},

  {id: "ASM985889v3", label: "SARS-CoV-2", faUrl:"https://s3.amazonaws.com/igv.org.genomes/ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna", downsampleTo: 100000},

  // Credit to Viral MSA for sequences below
  {id: "NC_039345", label: "Bombali Virus (Bombali ebolavirus)", downsampleTo: 100000},
  {id: "NC_014373", label: "Bundibugyo Virus (Bundibugyo ebolavirus)", downsampleTo: 100000},
  {id: "NC_001477", label: "Dengue Virus 1", downsampleTo: 100000},
  {id: "NC_001474", label: "Dengue Virus 2", downsampleTo: 100000},
  {id: "NC_001475", label: "Dengue Virus 3", downsampleTo: 100000},
  {id: "NC_002640", label: "Dengue Virus 4", downsampleTo: 100000},
  {id: "NC_002549", label: "Ebola Virus (Zaire ebolavirus)", downsampleTo: 100000},
  {id: "NC_004102", label: "HCV genotype 1", downsampleTo: 100000},
  {id: "NC_038882", label: "HCV genotype 1 (isolate H77)", downsampleTo: 100000},
  {id: "NC_009823", label: "HCV genotype 2", downsampleTo: 100000},
  {id: "NC_009824", label: "HCV genotype 3", downsampleTo: 100000},
  {id: "NC_009825", label: "HCV genotype 4", downsampleTo: 100000},
  {id: "NC_009826", label: "HCV genotype 5", downsampleTo: 100000},
  {id: "NC_009827", label: "HCV genotype 6", downsampleTo: 100000},
  {id: "NC_030791", label: "HCV genotype 7", downsampleTo: 100000},
  {id: "NC_001802", label: "HIV-1", downsampleTo: 100000},
  {id: "NC_001722", label: "HIV-2", downsampleTo: 100000},
  {id: "NC_063383", label: "Monkeypox Virus", downsampleTo: 100000},
  {id: "NC_004161", label: "Reston Virus (Reston ebolavirus)", downsampleTo: 100000},
  {id: "NC_006432", label: "Sudan Virus (Sudan ebolavirus)", downsampleTo: 100000},
  {id: "NC_014372", label: "Tai Forest Virus (Tai Forest ebolavirus, Cote d'Ivoire ebolavirus)", downsampleTo: 100000},


]



const backend = ""



function Alignment() {
  const defaultAcc = "ERR8254282";
  const defaultRef = referenceOptions[2];
  const [accession, setAccession] = useState('');
  const [refGenome, setRefGenome] = useState(referenceOptions[2]);
  const [genbankId, setGenbankId] = useState('');
  const faUrl = refGenome.faUrl ? refGenome.faUrl : refGenome.id === "custom" ? `https://genbank-api.vercel.app/api/genbank/${genbankId}?rettype=fasta` : `https://genbank-api.vercel.app/api/genbank/${refGenome.id}?rettype=fasta`


  
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
      const reference = refGenome.id == "ASM985889v3" ? {
          genome: refGenome.id,
        
        
        } : {fastaURL: faUrl, indexed: false}
      igv.createBrowser(igvDiv.current, {
        showSVGButton: true,
        showCursorTrackingGuide: true,
       
        
        //ref is /ref.fa
        reference: reference,
        locus: refGenome.id == "ASM985889v3" ? 'NC_045512.2:1-29903' : undefined,
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
      const response = await axios.post(`${backend}/align/${accession}?ref=${faUrl}&downsampleTo=${refGenome.downsampleTo}`);
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
    <div className="bg-gray-100 min-h-screen overflow-y-auto p-3">
      <div className="pb-3">
    <div className="flex flex-col items-center space-y-4">
      
      <h1 className="text-2xl font-bold mt-2">DeeperSeq</h1>
      <p className="text-gray-500">A tool for viewing deep viral sequencing data</p>
      {status !== 'complete' && status!== 'processing' && (
        <>
         <p
        className='text-center'>Enter a SRR/ERR accession. This will be mapped to the selected reference and then displayed<br />(Reads may be downsampled).</p>
      <form className="flex flex-col space-y-2" onSubmit={handleAlignmentSubmit}>
       <div>
        <label className="font-semibold block text-sm text-gray-600">
          SRA/ENA Accession ID:
          <input
            type="text"
            value={accession}
            onChange={handleAccessionChange}
            className="form-input mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-300 focus:ring focus:ring-blue-200 focus:ring-opacity-50"
          />
        </label>
        </div>
        <div>
        <label className="font-semibold block text-sm text-gray-600">
          Reference Genome:
          <select
            value={refGenome.id}
            onChange={(e) => setRefGenome(referenceOptions.find((o) => o.id === e.target.value))}
            className="form-select mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-300 focus:ring focus:ring-blue-200 focus:ring-opacity-50"
          >
            {referenceOptions.map((o) => (
              <option key={o.id} value={o.id}>
                {o.label}
              </option>
            ))}
          </select>
        </label>
        </div>
              {refGenome.id === "custom" && (
  <div>
    <label className="font-semibold block text-sm text-gray-600">
      Reference Genbank ID:
      <input
        type="text"
        value={genbankId}
        onChange={e => setGenbankId(e.target.value.trim())}
        className="form-input mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-300 focus:ring focus:ring-blue-200 focus:ring-opacity-50"
      />
    </label>
  </div>
)}


        <button type="submit" className="bg-blue-500 text-white px-4 py-2 rounded-md mt-2">
          Start Alignment
        </button>
        <div className="text-center pt-8 text-gray-500">
         {// load default accession with underlined text
          }
          <a href="#" onClick={() => {
            setAccession('SRR10903401');
            setRefGenome(defaultRef);
          }
          } className="underline">Fill form with example accessions</a>
        </div>
        <div className="text-center pt-1 text-gray-500">

          <a href="https://github.com/theosanderson/deeperseq/" className="underline">GitHub repo</a>
        </div>



      </form>
      </>
      )}
      {error && <p className="text-red-500">{error}</p>}
      {status === 'processing' && (
        <div className="flex flex-col items-center space-y-4">
          <ClipLoader color="#000" loading={true} size={150} />
          <h2 className="text-xl font-bold">Alignment Processing</h2>
          
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
