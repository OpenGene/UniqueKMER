#include "genomes.h"
#include "util.h"
#include "kmer.h"
#include "editdistance.h"
#include <sstream>
#include <memory.h>

// we use 512M memory
const int BLOOM_FILTER_LENGTH = (1<<29);

Genomes::Genomes(string faFile, Options* opt)
{
    mFastaReader = new FastaReader(faFile);
    mOptions = opt;
    mFastaReader->readAll();
}

Genomes::~Genomes()
{
    if(mFastaReader) {
        delete mFastaReader;
        mFastaReader = NULL;
    }
}

void Genomes::run() {
    map<string, string> genomes = mFastaReader->contigs();
    map<string, string>::iterator iter;
    mGenomeNum = 0;
    for(iter = genomes.begin(); iter != genomes.end() ; iter++) {
        if(mGenomeNum >= 100) {
            cerr << "UniqueKMER only supports up to 10 genomes, other genomes will be skipped." << endl;
            break;
        }
        mNames.push_back(iter->first);
        mSequences.push_back(iter->second);
        mUniqueKmers.push_back(vector<string>());
        mGenomeNum++;
    }

    buildKmerTable();
    makeUniqueKMER();
}

void Genomes::makeUniqueKMER() {
    unordered_map<uint64, int>::iterator iter;
    for(iter = mKmerTable.begin(); iter!=mKmerTable.end(); iter++) {
        uint64 key = iter->first;
        int id = iter->second;

        if(id >= 0) {
            string kmer = Kmer::seqFromUint64(key, mOptions->kmerKeyLen);
            mUniqueKmers[id].push_back(kmer);
        }
    }
}

void Genomes::buildKmerTable() {
    int keylen = mOptions->kmerKeyLen;
    int blankBits = 64 - 2*keylen;
    const int polyATailLen = 28;
    bool valid = true;
    for(uint32 i=0; i<mNames.size(); i++) {
        string& seq = mSequences[i];
        if(seq.length() < keylen)
            continue;
        // first calculate the first keylen-1 kmer
        // skip the polyA tail
        uint32 start = 0;
        uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        while(valid == false) {
            start++;
            key = Kmer::seq2uint64(seq, start, keylen-1, valid);
            // reach the tail
            if(start >= seq.length() - keylen - polyATailLen)
                return;
        }
        for(uint32 pos = start; pos < seq.length() - keylen - polyATailLen; pos++) {
            key = (key << 2);
            switch(seq[pos + keylen-1]) {
                case 'A':
                    key += 0;
                    break;
                case 'T':
                    key += 1;
                    break;
                case 'C':
                    key += 2;
                    break;
                case 'G':
                    key += 3;
                    break;
                case 'N':
                default:
                    // we have to skip the segments covering this N
                    pos++;
                    bool outterBreak = false;
                    key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                    while(valid == false) {
                        pos++;
                        key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                        // reach the tail
                        if(pos >= seq.length() - keylen - polyATailLen) {
                            outterBreak = true;
                            break;
                        }
                    }
                    if(outterBreak)
                        break;

                    continue;
            }
            key = (key << blankBits) >> blankBits;
            addKmer(key, i);
        }
    }
}

void Genomes::addKmer(uint64 key, int id) {
    unordered_map<uint64, int>::iterator iter = mKmerTable.find(key);
    if(iter == mKmerTable.end())
        mKmerTable[key] = id;
    else {
        if(iter->second != id)
            iter->second = -1;
    }

}

void Genomes::report() {
    for(int i=0; i<mGenomeNum; i++) {
        cerr << ">" << mNames[i] << endl;
        for(int k=0; k<mUniqueKmers[i].size(); k++) {
            cerr << k << ":" << mUniqueKmers[i][k] << endl;
        }
        cerr << endl;
    }
}

void Genomes::reportJSON(ofstream& ofs) {
}

void Genomes::reportHtml(ofstream& ofs) {
}