<p align="center"><img src=".images/multiply-logo.png" width="500"></p>

<p align="center">Multiplex PCR design, <i>in silico</i></p>

## Overview
`multiply` facilitates the design of optimised and reporducible multiplex PCRs. The diagram below summarises the `multiply` pipeline:
<br></br>
<p align="center"><img src=".images/multiply-pipeline.png" width="600"></p>

## Install
`multiply` has several software dependencies which can be installed using conda:
```
conda update conda
conda env create
conda activate multiply2
```
Afterwards, you can install the python package locally by running...
```
pip install .
```
...in the `multiply` directory.

## Basic usage
TODO

## Resources
`multiply` uses the following external software and databases:
- `primer3`. Individual primer pair design. https://primer3.org/
- `bedtools`. Genome arithmetic. https://bedtools.readthedocs.io/en/latest/
- `blastn`. Local alignment search. https://blast.ncbi.nlm.nih.gov/Blast.cgi
- **plasmodb**. *Plasmodium* reference genome. http://plasmodb.org/plasmo/
