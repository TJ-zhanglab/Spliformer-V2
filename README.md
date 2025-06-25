# Spliformer-V2

Spliformer-V2 is a deep-learning tool based on Segment-NT architecture that predicts RNA splicing usage accross 18 tissues based on pre-mRNA sequences(more details, see [paper](). Spliformer can take a VCF file containing variants of interest with its genotype as an input and predict acceptor/donor usage. 

## Installation

Our developing environment has been tested on ```Debian 6.1.140-1 (2025-05-22) x86_64```, Python version is ```python 3.8.19```, and Pytorch GPU version is ```pytorch 2.4.1```.
We encourage users to create a new conda environment before using Spliformer, you can first download miniconda through <https://docs.conda.io/en/latest/miniconda.html>, and then you can create a new conda environment named ```spliformerv2```  with all required package through the following commands:
```
conda env create -f spliformerV2.yml
``` 

## Usage
After git clone the repository, run the example through the following command:
```sh
cd Spliformer-V2
python run.py -I ./example/input/input19.vcf -O ./output/output.vcf -R ./reference/hg19.fa -A ./reference/grch37.txt  -M 0 -T LMC -G 1
```
**Required parameters**


-   -I: Input VCF with variants of interest.
-   -O: Output VCF with prediction of Spliformer-V2
-   -R: Reference fasta file of human genome. Please download and unzip it to the reference folder first before making prediction from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) or [GRCh38/hg38](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
-   -A: Annotation file of human genome.  We created hg19/hg38 genome annotation file according to the GENCODE v39 gtf file. The files locate in the [./reference/](https://github.com/TJ-zhanglab/Spliformer-V2/tree/main/reference).
-   -T: Tissues for prediction
-   -G: Provide genotype in the VCF (1: provided; 0: not provided)
**Optional parameters**


-   -M: Mask predicted scores with annotated acceptor/donor gain and unannotated acceptor/donor loss. ```0: not masked; 1: masked``` (default: 0).

Format of Spliformer-V2 INFO field in the VCF: ALT|gene|increased acceptor usage score|decreased acceptor usage score|increased donor usage score|decreased donor usage score|Increased acceptor dis|Decreased acceptor dis|Increased donor dis|Decreased donor dis

|Name                          |Disciption                         |
|-------------------------------|-----------------------------|
|ALT            |Alternate allele            |
|gene            |Gene name            |
|increased acceptor usage score| The increased usage of the acceptor site|
|decreased acceptor usage score| The decreased usage of the acceptor site|
|increased donor usage score| The increased usage of the donor site|
|decreased donor usage score| The decreased usage of the donor site|
|Increased acceptor dis|The distance of an acceptor with maximum increased score away from the variant|
|Decreased acceptor dis|The distance of an acceptor with maximum decreased score away from the variant|
|Increased donor dis|The distance of a donor with maximum increased score away from the variant|
|Decreased donor dis|The distance of a donor with maximum decreased score away from the variant|

## Cite us
If you use Spliformer-V2 for prediction, please cite [paper]()

