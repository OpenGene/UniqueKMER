# UniqueKMER
Generate unique KMERs for every contig in a FASTA file.

## how it works?
A `unique` KMER means a KMER key (such ATCGATCCTTAAGG) is only presented in one contig, and not presented at any other contig (in both forward and reverse strands). This tool accepts the input of a FASTA file consisting of many contigs, and extract unique KMERs for each contig.

# get this tool
## or compile from source
```shell
git clone https://github.com/OpenGene/UniqueKMER.git

# step 3: build
cd UniqueKMER
make
```
