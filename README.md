# pdb-profiling

![cover](docs/figs/cover.png)

Profiling Protein Structures from Protein Data Bank and integrate various resources.

## Goal

* Gather helpful/insightful indexes to evaluate a PDB structure's usefulness in:
  * Entry level
  * Model level (?)
  * Entity level
  * Chain level
  * Residue level
* Define the representative set of protein structures:
  * of a cluster with nearly identical sequences
  * of UniProt Entry
  * of UniProt Isoform
  * or any other assigned structure dataset
* Provide interface for ID/residue mapping
* Apply mature and robust API to collect well-organized data
  * Implement PDBe REST API (i.e <https://www.ebi.ac.uk/pdbe/api/doc/pdb.html>)
  * Implement PDBe Graph Database API (<https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/>)
* Release this project as a Python package

## Architecture Design

* Processer
* Scheduler
* Fethcer
* Pipelines
* Engine

![Software Architecture](docs/figs/Software_Architecture.png)

Processer -> Scheduler -> Fetcher -> Processer/Pipelines -> Scheduler

### Processer

> "..." stands for "can be extended by users"

* Accept a stream of IDs and their ID TYPE
  * UniProt 
  * PDB
  * RefSeq Nucleotide/Transcript/Protein
  * Ensembl Gene/Transcript/Protein
  * ...
* Prepared with per-defined procedures to choose
  * Towards ID Mapping
    * via `UniProt Retrieve/ID mapping`
    * via `G2S`
    * ...
  * Access specified PDB Data
    * via PDBe API
    * via Local Neo4j DataBase
    * ...
* Accept results from the `fetcher` passed by the `engine` and decide whether to pass the results to the `pipelines` directly or keep the results and generate new requests according to the procedure to the `scheduler` until the results to be integrated enough to be passed
  * Deliver requests to the `engine` (towards the `scheduler`)
    * initial
    * in the process
  * Deliver results to the `engine` (towards the `pipelines`)

### Scheduler

It collects requests to be sent to `fetcher`.

* Accept requests from the `processer`  (initial or add in the process)
  * Accept a stream of IDs and their ID TYPE
  * Accept corresponding pathway settings from the choosen procedure
* Schedule the tasks with queue (offer/poll)
* Deliver the requests to the `engine` (towards `fetcher`)

### Fethcer

* Fetch data via Internet (RESTful API) or local DataBase as defined by requests
* Deliver results to the `engine` (towards the `processer`)

### Pipelines

* Accept results from the `engine` (delivered by the `processer`)
* Perform data curation
  * Integrate data
  * Filtering
  * Ranking/Sorting
  * Scoring
* Output results
  * database
  * file system

## Copyright Notice

This project is developed by [Zefeng Zhu](https://github.com/NatureGeorge) and hold by [Minghui Group](https://lilab.jysw.suda.edu.cn/).