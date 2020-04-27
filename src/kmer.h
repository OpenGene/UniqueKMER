#ifndef KMER_H
#define KMER_H

// includes
#include "common.h"
#include <vector>
#include <unordered_map>
#include "fastareader.h"
#include "options.h"
#include <fstream>

using namespace std;

class Kmer
{
public:
    Kmer(string filename, Options* opt);
    ~Kmer();
    void init(string filename);
    bool add(uint64 kmer64);
    void report();
    double getMeanHit();
    string getPlotX();
    string getPlotY();
    int getKmerCount();
    void reportJSON(ofstream& ofs);

    static uint64 seq2uint64(string& seq, uint32 pos, uint32 len, bool& valid);
    static string seqFromUint64(uint64 key, uint32 len);
    static uint64 reverseComplement(uint64 key, uint32 len);
    static bool test();

private:
    unordered_map<uint64, uint32> mKmerHits;
    FastaReader* mFastaReader;
    unordered_map<uint64, string> mNames;
    unordered_map<uint64, string> mSequences;
    Options* mOptions;
};


#endif