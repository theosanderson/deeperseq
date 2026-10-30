import React, { useState,useRef , useEffect} from 'react';
import axios from 'axios';
import ClipLoader from 'react-spinners/ClipLoader';
import igv from "igv"
import './App.css';


const referenceOptions = [
   {id: "custom", label: "Custom", downsampleTo: 100000},
   {id: "-----", label: "-----", downsampleTo: 100000},

  {id: "ASM985889v3", label: "SARS-CoV-2", faUrl:"https://raw.githubusercontent.com/theosanderson/deeperseq/refs/heads/master/GCF_009858895.2_ASM985889v3_genomic.fna", downsampleTo: 100000},

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

const getFaUrl = (refGenome, genbankId) => refGenome.faUrl ? refGenome.faUrl : refGenome.id === "custom" ? `https://genbank-api.vercel.app/api/genbank/${genbankId}?rettype=fasta` : `https://genbank-api.vercel.app/api/genbank/${refGenome.id}?rettype=fasta`



function Alignment() {
  const defaultAcc = "ERR8254282";
  const defaultRef = referenceOptions[2];
  const [accession, setAccession] = useState('');
  const [inputMode, setInputMode] = useState('accession');
  const [fastqUrl, setFastqUrl] = useState('');
  const [fastqUrl2, setFastqUrl2] = useState('');
  const [displayName, setDisplayName] = useState('');
  const [needsRef, setNeedsRef] = useState(false);
  const [refGenome, setRefGenome] = useState(referenceOptions[2]);
  const [genbankId, setGenbankId] = useState('');
  const faUrl = getFaUrl(refGenome, genbankId);


  
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

  // ref and gbId are passed explicitly when starting from query params, before state has updated
  const doAlign = async (accession, ref = refGenome, gbId = genbankId) => {
    try {
      // check accession starts with SRR or ERR
      accession = accession.trim();
      if (!accession.startsWith('SRR') && !accession.startsWith('ERR')) {
        setError('Accession must start with SRR or ERR');
        return;
      }
      const response = await axios.post(`${backend}/align/${accession}?ref=${getFaUrl(ref, gbId)}&downsampleTo=${ref.downsampleTo}`);
      const taskID = response.data.task_id;
      setDisplayName(accession);
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

  const doAlignUrl = async (url1, url2, ref = refGenome, gbId = genbankId) => {
    try {
      const urls = [url1.trim(), url2.trim()].filter(u => u);
      if (urls.length === 0 || !urls.every(u => /^(https?|ftp):\/\//.test(u))) {
        setError('FASTQ URLs must start with http://, https:// or ftp://');
        return;
      }
      const params = new URLSearchParams({ url: urls[0], ref: getFaUrl(ref, gbId), downsampleTo: ref.downsampleTo });
      if (urls[1]) {
        params.append('url2', urls[1]);
      }
      const response = await axios.post(`${backend}/align_url?${params}`);
      if (response.data.error) {
        setError(response.data.error);
        return;
      }
      // the server names the sample itself; the BAM is served under that name
      setAccession(response.data.name);
      setDisplayName(urls.map(u => u.split('/').pop()).join(' + '));
      setTaskID(response.data.task_id);
      setStatus('processing');
      setError('');
      setLogs([]);
      pollTask(response.data.task_id);
    } catch (error) {
      console.error(error);
      setError('Error starting alignment task');
    }
  };

  const handleAlignmentSubmit = async (event) => {
    event.preventDefault();
    if (inputMode === 'url') {
      doAlignUrl(fastqUrl, fastqUrl2);
    } else {
      doAlign(accession);
    }
  };

  // if the query string has an accession (?acc=) or FASTQ URL(s) (?url=&url2=), start the alignment
  // ?ref= can be one of the reference ids above, or any other GenBank ID
  // FASTQ URLs without ?ref= only fill in the form, so that a reference can be chosen
  useEffect(() => {
    const urlParams = new URLSearchParams(window.location.search);
    const acc = urlParams.get('acc');
    const url = urlParams.get('url');
    const url2 = urlParams.get('url2') || '';
    const refParam = urlParams.get('ref');
    let ref = refGenome;
    let gbId = genbankId;
    if (refParam) {
      ref = referenceOptions.find((o) => o.id === refParam && o.id !== "-----") || referenceOptions[0];
      gbId = ref.id === "custom" ? refParam : '';
      setRefGenome(ref);
      setGenbankId(gbId);
    }
    if (acc) {
      setAccession(acc);
      doAlign(acc, ref, gbId);
    } else if (url) {
      setInputMode('url');
      setFastqUrl(url);
      setFastqUrl2(url2);
      if (refParam) {
        doAlignUrl(url, url2, ref, gbId);
      } else {
        setNeedsRef(true);
      }
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
        // check if "[ERROR]" in logs
        const error = response.data.log.filter((line) => line.includes('[ERROR]'));
        if (error.length > 0) {
          window.alert('There was an error processing this accession');
        }

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
        className='text-center'>Enter a SRR/ERR accession or the URL of a FASTQ file. This will be mapped to the selected reference and then displayed<br />(Reads may be downsampled).</p>
      <form className="flex flex-col space-y-2" onSubmit={handleAlignmentSubmit}>
       <div className="flex rounded-md border border-gray-300 overflow-hidden text-sm">
        {[['accession', 'SRA/ENA accession'], ['url', 'FASTQ URL']].map(([mode, label]) => (
          <button
            key={mode}
            type="button"
            onClick={() => { setInputMode(mode); setError(''); }}
            className={`flex-1 px-3 py-1 ${inputMode === mode ? 'bg-blue-500 text-white' : 'bg-white text-gray-600'}`}
          >
            {label}
          </button>
        ))}
       </div>
       {inputMode === 'accession' ? (
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
       ) : (
       <>
        <div>
        <label className="font-semibold block text-sm text-gray-600">
          FASTQ URL (optionally gzipped):
          <input
            type="url"
            value={fastqUrl}
            placeholder="https://example.com/reads.fastq.gz"
            onChange={e => setFastqUrl(e.target.value)}
            className="form-input mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-300 focus:ring focus:ring-blue-200 focus:ring-opacity-50"
          />
        </label>
        </div>
        <div>
        <label className="font-semibold block text-sm text-gray-600">
          Second FASTQ URL (optional, for paired-end reads):
          <input
            type="url"
            value={fastqUrl2}
            onChange={e => setFastqUrl2(e.target.value)}
            className="form-input mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-300 focus:ring focus:ring-blue-200 focus:ring-opacity-50"
          />
        </label>
        </div>
       </>
       )}
        {needsRef && (
          <p className="text-blue-700 text-sm">Choose the reference genome for these reads, then start the alignment.</p>
        )}
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
            setInputMode('accession');
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
          
          <h3 className="text-lg font-bold">{displayName || accession} reads</h3>
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
