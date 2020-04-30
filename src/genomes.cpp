#include "genomes.h"
#include "util.h"
#include "kmer.h"
#include "editdistance.h"
#include <sstream>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <atomic>
#include "sequence.h"

const int REPETITIVE_THRESHOLD = 20;

Genomes::Genomes(string faFile, Options* opt)
{
    cerr << "--- Load FASTA: " << faFile << endl;
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
        mUniqueKeys.push_back(vector<uint64>());
        mFilteredUniqueKeys.push_back(vector<uint64>());
        mGenomeNum++;
    }

    cerr << "--- Build KMER table" << endl;
    cerr << "--- Forward >>>"<<endl;
    buildKmerTable(false);
    cerr << "--- Reverse <<<"<<endl;
    buildKmerTable(true);
    cerr << "--- Find unique KMERs"<<endl;
    makeUniqueKMER();
}

void Genomes::makeUniqueKMER() {
    unordered_map<uint64, int>::iterator iter;
    for(iter = mKmerTable.begin(); iter!=mKmerTable.end(); iter++) {
        uint64 key = iter->first;
        int id = iter->second;

        if(id >= 0) {
            if(!mOptions->refFile.empty())
                mUniqueKeys[id].push_back(key);
            else {
                string kmer = Kmer::seqFromUint64(key, mOptions->kmerKeyLen);
                mUniqueKmers[id].insert(kmer);
            }
        }
    }
    // no longer needed
    mKmerTable.clear();
    
    if(!mOptions->refFile.empty())
        filterReferenceGenome();
}

void Genomes::filterReferenceGenome() {
    cerr << "--- Prepare to remove KMER keys that can be mapped to reference genome" << endl;
    cerr << "--- Load reference FASTA: " << mOptions->refFile << endl;
    FastaReader* hg = new FastaReader(mOptions->refFile);
    hg->readAll();
    map<string, string> contigs = hg->contigs();

    if(contigs.size() == 0)
        error_exit("Not a valid FASTA file: " + mOptions->refFile);

    map<string, string>::iterator ctgiter;
    vector<string> allseq;
    vector<string> allname;
    for(ctgiter = contigs.begin(); ctgiter != contigs.end() ; ctgiter++) {
        allseq.push_back(ctgiter->second);
        allname.push_back(ctgiter->first);
    }

    // init a table to indicate which key in reference need to be stored
    int keylen = min(14, mOptions->kmerKeyLen);
    cerr << "--- Calculate unique key coverage in " << keylen << " bp" << endl;
    
    size_t flagBufSize = 1L << (2*keylen);
    uint8* flagBuf = new uint8[flagBufSize];
    memset(flagBuf, 0, flagBufSize * sizeof(uint8));
    int blankBits = 64 - 2*keylen;

    uint64 mask = 0;
    for(int i=0; i<keylen*2; i++)
        mask = (mask<<1) + 0x01;

    #pragma omp parallel for
    for(int i=0; i<mGenomeNum; i++) {
        for(int j=0; j<mUniqueKeys[i].size(); j++) {
            uint64 key = mUniqueKeys[i][j];
            uint64 rckey = Kmer::reverseComplement(key, mOptions->kmerKeyLen);

            for(int offset = 0; offset <= mOptions->kmerKeyLen - keylen; offset++) {
                int shift = offset * 2;
                uint64 part = (key & (mask << shift)) >> shift;
                uint64 rcpart = (rckey & (mask << shift)) >> shift;
                flagBuf[part] = 1;
                flagBuf[rcpart] = 1;
            }
        }
    }

    cerr << "--- Index reference genome" << endl;

    unordered_map<uint32, vector<ContigPos>> keyCP;
    unordered_map<uint32, vector<ContigPos>>::iterator cpIter;

    for(int c = 0; c < allseq.size(); c++) {
        cerr << (c+1) << "/" << allseq.size() << ": " << allname[c] <<  endl;
        string seq = allseq[c];
        uint32 start = 0;
        bool valid;
        uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        while(valid == false) {
            start++;
            key = Kmer::seq2uint64(seq, start, keylen-1, valid);
            // reach the tail
            if(start >= seq.length() - keylen)
                return;
        }
        for(uint32 pos = start; pos < seq.length() - keylen; pos++) {
            key = (key << 2);
            switch(seq[pos + keylen-1]) {
                case 'A': key += 0; break;
                case 'T': key += 1; break;
                case 'C': key += 2; break;
                case 'G': key += 3; break;
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
                        if(pos >= seq.length() - keylen) {
                            outterBreak = true;
                            break;
                        }
                    }
                    if(outterBreak)
                        break;

                    continue;
            }
            key = (key << blankBits) >> blankBits;

            uint32 key32 = (uint32)key;
            if(flagBuf[key32] > 0 && flagBuf[key32]!= REPETITIVE_THRESHOLD) {
                flagBuf[key32]++;
                // this is a repetitive key, ignore and clear it
                if(flagBuf[key32]== REPETITIVE_THRESHOLD) {
                    keyCP.erase(key32);
                }  else if (flagBuf[key32]==1) {
                    // the first
                    keyCP[key32] = vector<ContigPos>();
                    keyCP[key32].push_back(ContigPos(c, pos));
                } else {
                    keyCP[key32].push_back(ContigPos(c, pos));
                }
            }
        }
    }

    cerr << "--- Filter KMER keys that can be aligned to reference genome" << endl;

    atomic_int finished;
    finished=0;
    #pragma omp parallel for
    for(int i=0; i<mGenomeNum; i++) {
        cerr << "Filtering " << (finished+1) << "/" << mGenomeNum << ": " << mNames[i] <<  endl;
        finished++;
        for(int j=0; j<mUniqueKeys[i].size(); j++) {
            uint64 key = mUniqueKeys[i][j];
            uint64 rckey = Kmer::reverseComplement(key, mOptions->kmerKeyLen);
            string keySeq = Kmer::seqFromUint64(key, mOptions->kmerKeyLen);
            string rckeySeq = Kmer::seqFromUint64(rckey, mOptions->kmerKeyLen);

            //bool debug = (keySeq ==  "GGGAAAAAAGCCGCGCGGGGGCGCC");

            bool mapped = false;
            for(int rc =0; rc<2; rc++) {
                if(mapped)
                    break;
                uint64 curkey = key;
                if(rc == 1)
                    curkey = Kmer::reverseComplement(key, mOptions->kmerKeyLen);
                
                int partHits = 0;
                bool repetitiveMode =  false;
                for(int offset = 0; offset <= mOptions->kmerKeyLen - keylen; offset++) {
                    if(mapped)
                        break;
                    int shift = offset*2;
                    uint64 part = (curkey & (mask << shift)) >> shift;
                    //if(debug)
                    // cerr << Kmer::seqFromUint64(part, keylen)<<endl;
                    uint32 part32 = (uint32) part;
                    if(flagBuf[part32] > 1) {
                        partHits++;
                        // all parts hit reference genome, mark it as mapped
                        if(partHits >= mOptions->kmerKeyLen - keylen + 1 && repetitiveMode)  {
                            mapped = true;
                            break;
                        }
                    }
                    // repetitive KMER of reference genome, skip them
                    if(flagBuf[part32] == REPETITIVE_THRESHOLD) {
                        repetitiveMode = true;
                        // all parts hit reference genome, mark it as mapped
                        if(partHits >= mOptions->kmerKeyLen - keylen + 1) {
                            mapped = true;
                            break;
                        } else
                            continue;
                    }
                    // hit the reference, do a mapping
                    if(flagBuf[part32] > 1) {
                        //if(debug)
                        // cerr << " hit"  << endl;
                        vector<ContigPos>& cplist = keyCP[part32];
                        for(int p=0; p<cplist.size(); p++) {
                            ContigPos cp = cplist[p];
                            uint32 pos = cp.position;
                            uint16 ctg = cp.contig;
                            int compStart = pos - (mOptions->kmerKeyLen - keylen - offset);
                            if(compStart <0 || compStart + mOptions->kmerKeyLen > allseq[ctg].size())
                                break;

                            /*if(debug) {
                                cerr << "compStart: " << compStart <<  endl;
                                cerr << "offset: " << offset <<  endl;
                                cerr << "p: " << p <<  endl;
                            }*/

                            int ed = 0;
                            if(rc == 1)
                                ed = edit_distance(allseq[ctg].c_str() + compStart, mOptions->kmerKeyLen, rckeySeq.c_str(), mOptions->kmerKeyLen);
                            else
                                ed =  edit_distance(allseq[ctg].c_str() + compStart, mOptions->kmerKeyLen, keySeq.c_str(), mOptions->kmerKeyLen);
                            if(ed <= mOptions->edThreshold) {
                                mapped = true;
                                break;
                            }
                        }
                    }
                }
            }

            if(!mapped) {
                mFilteredUniqueKeys[i].push_back(key);
            }
        }
        mUniqueKeys[i].clear();
    }

    delete[] flagBuf;
    delete hg;
    keyCP.clear();

    // change to string
    for(int i=0; i<mGenomeNum; i++) {
        for(int j=0; j<mFilteredUniqueKeys[i].size(); j++) {
            uint64 key = mFilteredUniqueKeys[i][j];
            string keySeq = Kmer::seqFromUint64(key, mOptions->kmerKeyLen);
            mUniqueKmers[i].insert(keySeq);
        }
    }

}

void Genomes::buildKmerTable(bool reversed) {
    int keylen = mOptions->kmerKeyLen;
    int blankBits = 64 - 2*keylen;
    bool valid = true;
    for(uint32 i=0; i<mNames.size(); i++) {
        string seq = mSequences[i];
        if(reversed) {
            Sequence s(seq);
            Sequence rc = ~s;
            seq = rc.mStr;
        }
        if(reversed)
            cerr << "<<< ";
        else
            cerr << ">>> ";
        cerr << (i+1) << "/" << mGenomeNum << ": " << mNames[i] <<  endl;
        if(seq.length() < keylen)
            continue;
        // first calculate the first keylen-1 kmer
        uint32 start = 0;
        uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        while(valid == false) {
            start++;
            key = Kmer::seq2uint64(seq, start, keylen-1, valid);
            // reach the tail
            if(start >= seq.length() - keylen)
                return;
        }
        for(uint32 pos = start; pos < seq.length() - keylen; pos++) {
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
                        if(pos >= seq.length() - keylen) {
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

string Genomes::outputKmer(int id, string& path, string& kmerFilename, int& count) {
    ofstream kmerofs;
    kmerofs.open(joinpath(path, kmerFilename));

    map<int, string> posSeq;
    set<string>::iterator setiter;

    string kmerStr;
    kmerStr += ">" + mNames[id] + "\n";

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
    count = 0;
    map<int, string>::iterator iter;
    for(iter = posSeq.begin(); iter != posSeq.end(); iter++) {
        int pos = iter->first;
        if(last >=0 && pos-last <= mOptions->kmerSpacing)
            continue;
        kmerofs << ">p" << alignToDigits(pos, digits) << endl;
        kmerofs << iter->second << endl;
        kmerStr += iter->second + "\n";
        last = pos;
        count++;
    }
    kmerofs.close();

    return kmerStr;
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
        ofs << endl <<  seq;
    }
    ofs.close();
}

extern string command;
void Genomes::output() {
    ofstream index;
    index.open( joinpath(mOptions->outdir, "index.html"), ifstream::out);
    ofstream kmercollection;
    kmercollection.open( joinpath(mOptions->outdir, "kmercollection.fasta"), ifstream::out);


    index<<"<HTML><head><title>UniqueKMER Report</title></head><body>" <<  endl;
    index<<"<h1>Unique KMER</h1>"<<endl;
    index<<"<div style='font-size:12px;color:#666666;'>Generated by <a href='https://github.com/OpenGene/UniqueKMER'>UniqueKMER v" << UNIQUEKMER_VER  << "</a>, <font color='#999999'>"<<endl;
    index<<" command: " << command  << "</font></div>"<<endl;
    index<<"<div><ul>" << endl;

    cerr << "--- Output unique KMER and genome FASTA files" << endl;
    vector<int> counts(mGenomeNum);
    vector<string> pathes(mGenomeNum);
    vector<string> kmerFilenames(mGenomeNum);
    vector<string> genomeFilenames(mGenomeNum);
    vector<string> kmerStrings(mGenomeNum);

    // generate folders
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

        pathes[i]=path;
        kmerFilenames[i]=kmerFilename;
        genomeFilenames[i]=genomeFilename;
    }

    // process and  output files
    atomic_int finished;
    finished=0;
    #pragma omp parallel for
    for(int i=0; i<mGenomeNum; i++) {
        string kmerFilename = kmerFilenames[i];
        string genomeFilename = genomeFilenames[i];
        string path = pathes[i];
        int count;
        kmerStrings[i] = outputKmer(i, path, kmerFilename, count);
        counts[i] = count;
        outputGenome(i, path, genomeFilename);
        cerr << "Output " << (finished+1) << "/" << mGenomeNum << ": " << mNames[i] << "." << " unique: " << counts[i] << endl;
        finished++;
    }

    // output index.html and kmercollection.fasta
    for(int i=0; i<mGenomeNum; i++) {
        int contigSize = mSequences[i].size();
        string folder = to_string(contigSize % 100);
        string kmerFilename = kmerFilenames[i];
        string genomeFilename = genomeFilenames[i];
        string path = pathes[i];
        int count = counts[i];

        string color;
        if(count  > 0)
            color = "blue";
        else
            color = "#333333";

        index << "<li>" <<  mNames[i] <<  " (" << count << " unique)";
        index << "  &nbsp;<a style='color:" << color << "' href='genomes_kmers/" << folder << "/" << kmerFilename << "'>KMER file</a>";
        index << "&nbsp; | &nbsp;<a style='color:" << color << "' href='genomes_kmers/" << folder << "/" << genomeFilename << "'>Genome file</a>";
        index << " </li>" << endl;

        if(count>0)
            kmercollection << kmerStrings[i];
    }

    index << "</ul></div></body></html>" << endl;

    index.close();
    kmercollection.close();
}
