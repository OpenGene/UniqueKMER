#ifndef GENOMES_H
#define GENOMES_H

// includes
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include "common.h"
#include "fastareader.h"
#include <vector>
#include <list>
#include <set>
#include <unordered_map>
#include "options.h"

using namespace std;

class MapResult{

public:
    MapResult(){
        mapped = false;
        start = 0;
        len = 0;
        ed = 0x7FFFFFFF; // initialized with a very large ED
    }
public:
    bool mapped;
    uint32 start;
    uint32 len;
    uint32 ed; // edit distance
};

#pragma pack(2) 
// if contig is -1, means this is a dupe entry, and position will be the position in the dupList
struct ContigPos{
    ContigPos(uint16 ctg, uint32 pos) {
        contig = ctg;
        position  = pos;
    }
    uint16 contig;
    uint32 position;
};
#pragma pack() 

class Genomes
{
public:
    Genomes(string fastaFile, Options* opt);
    ~Genomes();
    void run();
    void output();

private:
    void buildKmerTable(bool reversed);
    void addKmer(uint64 key, int id, bool reversed);
    void makeUniqueKMER();
    int outputKmer(int id, string& path, string& filename);
    void outputGenome(int id, string& path, string& filename);
    string alignToDigits(int val, int digits);
    void filterReferenceGenome();

private:
    int mGenomeNum;
    FastaReader* mFastaReader;
    vector<string> mSequences;
    vector<string> mNames;
    vector<set<string>> mUniqueKmers;
    vector<vector<uint64>> mUniqueKeys;
    vector<vector<uint64>> mFilteredUniqueKeys;
    // KMER: value
    // If value is negative, means more than one contig hit this KMER
    unordered_map<uint64, int> mKmerTable;
    Options* mOptions;
};


#endif