# UniqueKMER
Generate unique KMERs for every contig in a FASTA file.  

Unique KMER is consisted of KMER keys (i.e. ATCGATCCTTAAGG) that are only presented in one contig, but not presented in any other contigs (for both forward and reverse strands).  

This tool accepts the input of a FASTA file consisting of many contigs, and extract unique KMERs for each contig.

The output unique KMER file and Genome file can be used for fastv: https://github.com/OpenGene/fastv, which is an ultra-fast tool to identify and visualize microbial sequences from sequencing data.

# what does UniqueKMER output?
This tool outputs a folder (folder name can be specified by `-o/--outdir`), which contains a `index.html` and a subfolder `genomes_kmers`. The subfolder `genomes_kmers` contains a KMER file and a Genome file for each contig, both in FASTA format. You can open the `index.html` with any browser, then click on the contig names to find its KMER file and Genome file.
* a small example: http://opengene.org/uniquekmer/test/index.html

# get this tool
```shell
git clone https://github.com/OpenGene/UniqueKMER.git
cd UniqueKMER
make
# install to system (may need sudo)
make install
```
# simple example:
```shell
uniquekmer -f test.fasta
```
You can get the test.fasta from: http://opengene.org/test.fasta

# more examples
### set the KMER key length
```shell
# 16-mer (i.e. ATCGATCGATCGATCG...)
uniquekmer -f test.fasta -k 16
```
### filter the KMER keys that can be mapped to a reference genome (i.e. human genome)
```shell
# KMER sequences that can be mapped to hg38 with `edit distance <=2`  will be removed
uniquekmer -f test.fasta -r hg38.fasta -e 2
```
### set the spacing to avoid many continuous KMER keys
```shell
# the spacing will be 2, which means if `key(pos)` is stored, then `key(pos+1)`  and `key(pos+2)` will be skipped
uniquekmer -f test.fasta -s 2
```

options:
```shel
  -f, --fasta            FASTA input file name (string)
  -o, --outdir           Directory for output. Default is unique_kmers in the current directory. (string [=unique_kmers])
  -k, --kmer             The length k of KMER (10~32), default 25 (int [=25])
  -s, --spacing          If a key with POS is recorded, then skip [POS+1...POS+spacing] to avoid too compact result (0~100). default 0 means no skipping. (int [=0])
  -g, --genome_limit     Process up to genome_limit genomes in the FASTA input file. Default 0 means no limit. This option is for DEBUG. (int [=0])
  -r, --ref              Reference genome FASTA file name. Specify this only when you want to filter out the unique KMER that can be mapped to reference genome. (string [=])
  -e, --edit_distance    KMER mapped to reference genome with edit distance <= edit_distance will be removed (0~16). 3 for default. (int [=3])
  -?, --help             print this message
```
