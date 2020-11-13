# pdb-profiling

[![Build](https://img.shields.io/travis/naturegeorge/pdb-profiling?style=for-the-badge&logo=travis)](https://github.com/naturegeorge/pdb-profiling)
[![SupportPythonVersion](https://img.shields.io/pypi/pyversions/pdb-profiling.svg?style=for-the-badge&logo=python)](https://pypi.org/project/pdb-profiling/)
[![Version](https://img.shields.io/pypi/v/pdb-profiling?style=for-the-badge&logo=PYPI)](https://github.com/naturegeorge/pdb-profiling/blob/master/pdb_profiling/__init__.py)
[![Dependencies](https://img.shields.io/librariesio/github/NatureGeorge/pdb-profiling?style=for-the-badge&logo=PYPI)](https://github.com/naturegeorge/pdb-profiling/blob/master/setup.py)
[![PyPIDownloads](https://img.shields.io/pypi/dm/pdb-profiling.svg?style=for-the-badge&logo=PYPI)](https://pypi.org/project/pdb-profiling/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge&logo=github)](https://github.com/naturegeorge/pdb-profiling/blob/master/LICENSE)
[![GitHubDownloads](https://img.shields.io/github/downloads/NatureGeorge/pdb-profiling/total?style=for-the-badge&logo=github)](https://github.com/NatureGeorge/pdb-profiling/releases/)
[![Coverage Status](https://img.shields.io/coveralls/github/NatureGeorge/pdb-profiling?style=for-the-badge&logo=coveralls)](https://coveralls.io/github/NatureGeorge/pdb-profiling?branch=master)

![cover](https://user-images.githubusercontent.com/43134199/95018149-58cfc200-0690-11eb-9e64-760faec5130f.png)

Profiling Protein Structures from Protein Data Bank and integrate various resources.

## Features

* `Collection`: Implement various API to collect the well-organized metadata of PDB in real time.
* `Compatibility`: Provide a unified call for API-interface and return-data-form as well as subsequent data processing.
* `Detection`: Reorganize metadata to evaluate a PDB structure in Entry-Assembly/Model-Entity-Chain-Residue level and integrated with UniProt-KB.
* `Interaction`: Include UniProt Isoform Interaction in Asymmetric unit plus Biological Assembly level.
* `Selection`: Define the representative set of PDB structures in Monomeric|Homomeric|Heteromeric states.
* `Mapping`: Provide interface for both entry-identifier/accession-level and residue-level bidirectional mapping.

## Install

> Notice: require Python Environment >= 3.6, Platform Independent

Install by `pip` command.

### *Before your Installation

* Make sure that your 64-bit machine is installed with 64-bit Python.
* To avoid some unexpected issues, you should upgrade your `pip` beforehand:

```bash
pip install --upgrade pip
``` 

### Official Installation

```bash
pip install pdb-profiling
```

If you have already installed an older version of `pdb-profiling`, use the following command to install the latest version:

```bash
pip install --upgrade pdb-profiling
```

## Documentation

<https://pdb-profiling.netlify.app/>

## Examples

See `examples/...`

1. [Introduction](https://nbviewer.jupyter.org/github/NatureGeorge/pdb-profiling/blob/master/examples/Introduction.ipynb)
2. [Batch](https://nbviewer.jupyter.org/github/NatureGeorge/pdb-profiling/blob/master/examples/Batch.ipynb)
3. [DisplayPDB](https://nbviewer.jupyter.org/github/NatureGeorge/pdb-profiling/blob/master/examples/DisplayPDB.ipynb)
4. ...

## Resources

* PDBe Entry-Based API
* PDBe Aggregated API (PDBe Graph API)
* PDBe ModelServer API
* SWISS-MODEL Repository API
* UniProt API
* EBI Proteins API
* Interactome3D API
* ...

> click [here](https://pdb-profiling.netlify.app/docs/5-reference/) for more details

## Related Resources

> Using similar data resources but meant to achieve different goals.

<details>

<summary>Click to view</summary>

* `RCSB`
  * [Build Customize Tabular Reports of PDB Data](https://www.rcsb.org/news?year=2020&article=5f6529e207302466657ec0e9&feature=true)
  * [RCSB PDB Search API](http://search.rcsb.org/)
    * [Documentation for New and Improved APIs](https://www.rcsb.org/news?year=2020&article=5f65165507302466657ec0e8&feature=true)
* [MolArt](https://github.com/davidhoksza/MolArt)

</details>

## Copyright Notice

This project is developed by [Zefeng Zhu](https://github.com/NatureGeorge) and hold by [Minghui Group](https://lilab.jysw.suda.edu.cn/).