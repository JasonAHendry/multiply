<p align="center"><img src=".images/multiply-logo.png" width="500"></p>

<p align="center">Multiplex PCR design, <i>in silico</i></p>

## Overview
`multiply` is a command-line tool enabling the design of multiplexed PCRs for a user-specified set of target genes and/or regions. It works by first producing a set of candidate primers for each target using primer3 (`multiply generate`). It then computes the number of SNPs in each primer (`multiply snpcheck`); potential dimers betweeen pairs of primers (`multiply align`); and potential mispriming and off-target amplicons for each primer (`multiply blast`). Information from these three quality control steps is passed to a cost function, which is minimised by brute-force or with a greedy search algorithm (`multiply select`). 

The pipeline is summarised below:
<br></br>
<p align="center"><img src=".images/multiply-pipeline.png" width="700"></p>

## Install
First, clone the repository to your local machine:

```
git clone https://github.com/JasonAHendry/multiply
```
Then, install the software dependencies using conda:

```
cd multiply
conda update conda
conda env create
conda activate multiply2
```
Finally, install `multiply` itself with pip:

```
pip install .
```

Test installation by running:

```
multiply
```

## Basic usage
TODO

## Known limitations
TODO

## Resources
`multiply` uses the following external software and databases:
- `primer3`. Individual primer pair design. https://primer3.org/
- `bedtools`. Genome arithmetic. https://bedtools.readthedocs.io/en/latest/
- `blastn`. Local alignment search. https://blast.ncbi.nlm.nih.gov/Blast.cgi
- **plasmodb**. *Plasmodium* reference genome. http://plasmodb.org/plasmo/
