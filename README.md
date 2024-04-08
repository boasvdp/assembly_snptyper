# assembly_snptyper
Type assembly based on a VCF using minimap2 + samtools

This software can take a bacterial genome assembly in fasta format, align it to a reference genome and search for specific SNPs based on the alignment. It uses minimap2 and samtools to do this. You would typically only run this tool instead of a read mapping tool if you don't have access to the sequence reads.

On a laptop, the tool can type 4000 *S. pyogenes* for M1UK SNPs in 10m30s using 4 processes and ~600 Mb memory.

The tool was written with M1UK typing in mind. In this specific case, the reference genome is an M1global strain, from which M1UK can be discerned by only 27 lineage-specific SNPs. Important assumptions here are that the type of interest (here, M1UK) is closely related to the reference genome, warranting the use of the asm5 preset for minimap2.

For some pathogens, split k-mer typing (https://github.com/bacpop/ska.rust) might also be an option. However, this has the assumption that there cannot be any sequence variation around type SNPs within a kmer's length. In the case of M1UK, a MNP in positions 116162-116163 complicates this.

## Installation

Clone the repository and install locally using pip:

```
# clone repo
git clone https://github.com/boasvdp/assembly_snptyper.git
cd assembly_snptyper

# install the conda environment for minimap2, samtools and pandas
mamba env create -f env.yaml
conda activate env_assembly_snptyper

# install assembly_snptyper using pip
python -m pip install .
assembly_snptyper --version
```

## Usage

Note that: 
- the reference VCF must be based on the reference genome (i.e. identical chromosome names).
- the list of input files must be a file list of paths to assemblies. This is to accomodate running the script on large numbers of files

```
assembly_snptyper --help
usage: assembly_snptyper [-h] --vcf VCF --reference REFERENCE --list_input LIST_INPUT [-p PROCESSES] [-v] [--version]

options:
  -h, --help            show this help message and exit
  --vcf VCF             VCF file with variants that determine the type (default: None)
  --reference REFERENCE
                        Reference genome (default: None)
  --list_input LIST_INPUT
                        List of input assemblies (default: None)
  -p PROCESSES, --processes PROCESSES
                        Number of processes passed to multiprocessing (default: 1)
  -v, --verbose         Verbose output, can be used multiple times. 0 = warning, 1 = info, 2 = debug (default: 0)
  --version             show program's version number and exit
```