#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

using namespace std;

class Options{
public:
    Options();
    void init();
    bool validate();

public:
    // file name of FASTA file
    string faFile;
    // file name of human genome assembly
    string hgFile;
    // output folder name
    string outdir;
    // KMER len
    int kmerKeyLen;

};

#endif