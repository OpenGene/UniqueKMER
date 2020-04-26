# UniqueKMER
Generate unique KMERs for every contig in a FASTA file.  

Unique KMER is consisted of KMER keys (i.e. ATCGATCCTTAAGG) that are only presented in one contig, but not presented in any other contigs (for both forward and reverse strands).  

This tool accepts the input of a FASTA file consisting of many contigs, and extract unique KMERs for each contig.

# get this tool
## compile from source
```shell
git clone https://github.com/OpenGene/UniqueKMER.git

# step 3: build
cd UniqueKMER
make
```
# usage
Simple example:
```shell
uniquekmer -f test.fasta
```
You can get the test.fasta from: http://opengene.org/test.fasta

Options:
```shel
  -f, --fasta           FASTA input file name (string)
  -o, --outdir          Directory for output. Default is unique_kmers in the current directory. (string [=unique_kmers])
  -k, --kmer            The length k of KMER (10~32), default 25 (int [=25])
  -s, --spacing         If a key with POS is recorded, then skip [POS+1...POS+spacing] to avoid too compact result (0~100). default 0 means no skipping. (int [=0])
  -g, --genome_limit    Process up to genome_limit genomes in the FASTA input file. Default 0 means no limit. This option is for DEBUG. (int [=0])
  -?, --help            print this message
```
