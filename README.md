# 🧬 Spliformer-V2

**Spliformer-V2** is a deep learning model that is based on the Segment-NT architecture to fine-tune on paired genome/transcriptome datasets for predicting RNA splicing usage across 18 human tissues using pre-mRNA sequences (more details, see [our preprint paper]().Spliformer-v2 includes 12 CNS tissues and 6 non-CNS tissues, representing so far the most complete splicing prediciton model for human brain regions and spinal cord tissues. It supports VCF input with genotypes and outputs tissue-specific changes in splice site usage.

## ✨ Features

- Predicts splice site usage alterations from genetic variants  
- Supports genotype-aware prediction  
- Tissue-specific modeling across 18 human tissues (including 12 CNS tissues) 
- Accepts standard VCF inputs  
- Output in annotated VCF format
  
## 🛠 Installation
The development environment has been tested on:

OS: Debian 6.1.140-1 (2025-05-22) x86_64

Python: 3.8.19

PyTorch (GPU): 2.4.1

### Step-by-step:
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Clone the repository and install dependencies:
    ```bash
    git clone https://github.com/TJ-zhanglab/Spliformer-V2.git
    cd Spliformer-V2
    conda env create -f spliformerV2.yml
    conda activate spliformerv2
    ```

## 🚀 Usage

After setting up the environment, you can run an example prediction:
```bash
python run.py \
  -I ./example/input/input19.vcf \
  -O ./output/output.vcf \
  -R ./reference/hg19.fa \
  -A ./reference/grch37.txt \
  -T LMC \
  -G 1 \
  -M 0
```
**⚙️ Required parameters**


-   -I: Input VCF with variants of interest.
-   -O: Output VCF with prediction of Spliformer-V2
-   -R: Reference fasta file of human genome. Please download and unzip it to the reference folder first before making prediction from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) or [GRCh38/hg38](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
-   -A: Annotation file of human genome.  We created hg19/hg38 genome annotation file according to the GENCODE v39 gtf file. The files locate in the [./reference/](https://github.com/TJ-zhanglab/Spliformer-V2/tree/main/reference).
-   -T: Tissues for prediction
-   -G: Provide genotype in the VCF (1: yes; 0: no)
  
**⚙️ Optional parameters**


-   -M: Mask predicted scores with annotated acceptor/donor gain and unannotated acceptor/donor loss. ```0: not masked; 1: masked``` (default: 0).

## 📤 Output VCF Format
The prediction results are stored in the INFO field of the output VCF in the following format: **ALT|gene|increased acceptor usage score|decreased acceptor usage score|increased donor usage score|decreased donor usage score|Increased acceptor dis|Decreased acceptor dis|Increased donor dis|Decreased donor dis**

|Field                          |Disciption                         |
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

## 📚 Cite us
If you use Spliformer-V2 for in your work, please cite [paper]()

