#include "genomes.h"
#include "util.h"
#include "kmer.h"
#include "editdistance.h"
#include <sstream>
#include <memory.h>
#include <math.h>
#include "sequence.h"

// we use 512M memory
const int BLOOM_FILTER_LENGTH = (1<<29);

Genomes::Genomes(string faFile, Options* opt)
{
    cerr << "------Load FASTA: " << faFile << endl;
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
        if(mOptions->genomeLimit>0 && mGenomeNum >= mOptions->genomeLimit) {
            cerr << "UniqueKMER only supports up to " << mOptions->genomeLimit << " genomes, other genomes will be skipped." << endl;
            break;
        }
        mNames.push_back(iter->first);
        mSequences.push_back(iter->second);
        mUniqueKmers.push_back(set<string>());
        mGenomeNum++;
    }

    cerr << "------Build KMER table" << endl;
    cerr << "------Forward >>>"<<endl;
    buildKmerTable(false);
    cerr << "------Reverse <<<"<<endl;
    buildKmerTable(true);
    cerr << "------Find unique KMERs...";
    makeUniqueKMER();
    cerr << "done." << endl;
}

void Genomes::makeUniqueKMER() {
    unordered_map<uint64, int>::iterator iter;
    for(iter = mKmerTable.begin(); iter!=mKmerTable.end(); iter++) {
        uint64 key = iter->first;
        int id = iter->second;

        if(id >= 0) {
            string kmer = Kmer::seqFromUint64(key, mOptions->kmerKeyLen);
            mUniqueKmers[id].insert(kmer);
        }
    }
}

void Genomes::buildKmerTable(bool reversed) {
    int keylen = mOptions->kmerKeyLen;
    int blankBits = 64 - 2*keylen;
    const int polyATailLen = 28;
    bool valid = true;
    for(uint32 i=0; i<mNames.size(); i++) {
        string seq = mSequences[i];
        if(reversed) {
            Sequence s(seq);
            Sequence rc = ~s;
            seq = rc.mStr;
        }
        cerr << (i+1) << "/" << mGenomeNum << ": " << mNames[i] <<  endl;
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
            addKmer(key, i, reversed);
        }
    }
}

void Genomes::addKmer(uint64 key, int id, bool reversed) {
    unordered_map<uint64, int>::iterator iter = mKmerTable.find(key);
    if(iter == mKmerTable.end()) {
        // if it's reverse complemented sequence, we just use it to mark the KMER not unique, but dot add it as a new unique
        if(!reversed)
            mKmerTable[key] = id;
    }
    else {
        if(iter->second != id)
            iter->second = -1;
    }

}

string Genomes::alignToDigits(int val, int digits) {
    string s = to_string(val);
    while(s.length() < digits)
        s = "0"+s;

    return s;
}

int Genomes::outputKmer(int id, string& path, string& kmerFilename) {
    ofstream kmerofs;
    kmerofs.open(joinpath(path, kmerFilename));

    map<int, string> posSeq;
    set<string>::iterator setiter;

    if(mUniqueKmers[id].size() < mSequences[id].length()/300) {
        for(setiter=mUniqueKmers[id].begin(); setiter!=mUniqueKmers[id].end(); setiter++) {
            string seq = *setiter;
            int pos = mSequences[id].find(seq);
            if (pos != string::npos) {
                posSeq[pos]=seq;
            } else {
                cerr << "not found: " << seq << endl;
            }
        }
    } else {
        string& genome = mSequences[id];
        for(int pos=0; pos<genome.length() - mOptions->kmerKeyLen; pos++) {
            string seq = genome.substr(pos, mOptions->kmerKeyLen);
            setiter = mUniqueKmers[id].find(seq);
            if(setiter != mUniqueKmers[id].end()) {
                posSeq[pos] =  seq;
                mUniqueKmers[id].erase(setiter);
            }
        }
    }

    int digits = ceil(log10(mSequences[id].length()+1.0));

    int last = -1;
    int count = 0;
    map<int, string>::iterator iter;
    for(iter = posSeq.begin(); iter != posSeq.end(); iter++) {
        int pos = iter->first;
        if(last >=0 && pos-last <= mOptions->kmerSpacing)
            continue;
        kmerofs << ">p" << alignToDigits(pos, digits) << endl;
        kmerofs << iter->second << endl;
        last = pos;
        count++;
    }
    kmerofs.close();

    return count;
}

void Genomes::outputGenome(int id, string& path, string& genomeFilename) {
    ofstream ofs;
    ofs.open(joinpath(path, genomeFilename));

    ofs << ">" << mNames[id];

    const int line =  80;
    size_t seqlen = mSequences[id].length();
    size_t finished = 0;
    while(finished < seqlen) {
        size_t len = min((size_t)80, seqlen - finished);
        string seq = mSequences[id].substr(finished, len);
        finished += len;
        ofs << endl << seq;
    }
    ofs.close();
}

void Genomes::output() {
    ofstream index;
    index.open( joinpath(mOptions->outdir, "index.html"), ifstream::out);

    index<<"<HTML><head><title>UniqueKMER Report</title></head><div><ul>" << endl;

    cerr << "------Output unique KMER and genome FASTA files" << endl;
    for(int i=0; i<mGenomeNum; i++) {
        int contigSize = mSequences[i].size();
        string folder = to_string(contigSize % 100);
        string path = joinpath(mOptions->outdir,"genomes_kmers");
        if(!file_exists(path)) {
            if(mkdir(path.c_str(), 0755) != 0)
                error_exit("Failed to create directory: " + path);
        } else {
            if(!is_directory(path))
                error_exit("Not a directory: " + path);
        }
        path = joinpath(path,folder);
        if(!file_exists(path)) {
            if(mkdir(path.c_str(), 0755) != 0)
                error_exit("Failed to create directory: " + path);
        } else {
            if(!is_directory(path))
                error_exit("Not a directory: " + path);
        }

        string kmerFilename = str_keep_valid_filename(mNames[i]) + ".kmer.fasta";
        string genomeFilename = str_keep_valid_filename(mNames[i]) + ".fasta";

        string color;
        int unique = mUniqueKmers[i].size();
        if(unique  > 0)
            color = "blue";
        else
            color = "#333333";

        int count = outputKmer(i, path, kmerFilename);

        index << "<li>" <<  mNames[i] <<  " (" << count << " unique)";
        index << "  &nbsp;<a style='color:" << color << "' href='genomes_kmers/" << folder << "/" << kmerFilename << "'>KMER file</a>";
        index << "&nbsp; | &nbsp;<a style='color:" << color << "' href='genomes_kmers/" << folder << "/" << genomeFilename << "'>Genome file</a>";
        index << " </li>" << endl;
        
        outputGenome(i, path, genomeFilename);
        cerr << (i+1) << "/" << mGenomeNum << ": " << mNames[i] << "." << " unique: " << count << endl;
    }

    index << "</ul></div></body></html>" << endl;

    index.close();
}
