# pdb-profiling

[![DOI](https://zenodo.org/badge/247475852.svg)](https://zenodo.org/badge/latestdoi/247475852)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat&logo=github&colorB=5A65B3)](https://github.com/naturegeorge/pdb-profiling/blob/master/LICENSE)
[![SupportPythonVersion](https://img.shields.io/pypi/pyversions/pdb-profiling.svg?style=flat&logo=python&colorB=5A65B3)](https://pypi.org/project/pdb-profiling/)
[![Version](https://img.shields.io/pypi/v/pdb-profiling?style=flat&logo=PYPI&colorB=5A65B3)](https://github.com/naturegeorge/pdb-profiling/blob/master/pdb_profiling/__init__.py)
[![PyPIDownloads](https://img.shields.io/pypi/dm/pdb-profiling.svg?style=flat&logo=PYPI)](https://pypi.org/project/pdb-profiling/)
[![GitHubDownloads](https://img.shields.io/github/downloads/NatureGeorge/pdb-profiling/total?style=flat&logo=github)](https://github.com/NatureGeorge/pdb-profiling/releases/)
[![Build](https://img.shields.io/travis/naturegeorge/pdb-profiling?style=flat&logo=travis)](https://github.com/naturegeorge/pdb-profiling)
[![Coverage Status](https://img.shields.io/coveralls/github/NatureGeorge/pdb-profiling?style=flat&logo=coveralls)](https://coveralls.io/github/NatureGeorge/pdb-profiling?branch=master)

![cover](https://user-images.githubusercontent.com/43134199/95018149-58cfc200-0690-11eb-9e64-760faec5130f.png)

Profiling Protein Structures from Protein Data Bank and integrate various resources.

## Features

* `Collection`: Implement various API to collect the well-organized metadata of PDB in real time.
* `Integration`: Provide a unified call for API-interface and return-data-form as well as subsequent data processing.
* `Detection`: Reorganize metadata to evaluate a PDB structure in Entry-Assembly/Model-Entity-Chain-Residue level and integrated with UniProt-KB.
* `Interaction`: Include UniProt Isoform Interaction in Asymmetric unit plus Biological Assembly level.
* `Selection`: Define the representative set of PDB structures in Monomeric|Homomeric|Heteromeric states.
* `Mapping`: Provide interface for both entry-identifier/accession-level and residue-level bidirectional mapping.

## Install

> Notice: require Python Environment >= 3.7, Platform Independent

Install by `pip` command.

### *Before your Installation

* Make sure that your 64-bit machine is installed with 64-bit Python.
* To avoid some unexpected issues, you should upgrade your `pip` beforehand:

```bash
python -m pip install --upgrade pip
``` 

### Official Installation

```bash
python -m pip install pdb-profiling
```

If you have already installed an older version of `pdb-profiling`, use the following command to install the latest version:

```bash
python -m pip install --upgrade pdb-profiling
```

### Build From Source (optional, for non-windows environment)

```bash
python -m pip install cython
git clone https://github.com/NatureGeorge/pdb-profiling.git
python setup.py build_ext --inplace  # Need GCC or Other Compiler For C
python setup.py install              # or "sudo python setup.py install" or "python setup.py install --user"
```

## Documentation

<https://pdb-profiling.netlify.app/>

## Examples

### Basic Usage

* [Command Line Example](https://github.com/NatureGeorge/pdb-profiling/discussions/2)
* [Retrieve Bound Molecule Data From PDBe](https://github.com/NatureGeorge/pdb-profiling/discussions/3)
* ...

### Large-Scale-Example

* [ExAC](https://github.com/NatureGeorge/pdb-profiling/blob/master/examples/exac_example.md)

## Resources

* PDBe Entry-Based API
* PDBe Aggregated API (PDBe Graph API)
* PDBe ModelServer API
* SWISS-MODEL Repository API
* UniProt API
* EBI Proteins API
* RCSB Data API
* RCSB Search API
* Eutils API (minimum usage)
* ...

> click [here](https://pdb-profiling.netlify.app/docs/5-reference/) for more details

## Related Resources

> Using similar data resources but meant to achieve different goals.

<details>

<summary>Click to view</summary>

* `RCSB`: [Build Customize Tabular Reports of PDB Data](https://www.rcsb.org/news?year=2020&article=5f6529e207302466657ec0e9&feature=true)
* [MolArt](https://github.com/davidhoksza/MolArt)

</details>

## Copyright Notice

This project is developed by [Zefeng Zhu](https://github.com/NatureGeorge) and hold by [Minghui Group](https://lilab.jysw.suda.edu.cn/).