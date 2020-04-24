#include <stdio.h>
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "common.h"
#include "util.h"
#include "options.h"
#include "processor.h"

// TODO: code refactoring to remove these global variables
string command;
mutex logmtx;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cerr << "UniqueKMER: Generate unique KMERs for every contig in a FASTA file." << endl << "version " << UNIQUEKMER_VER << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "uniquekmer " << UNIQUEKMER_VER << endl;
        return 0;
    }

    cmdline::parser cmd;
    cmd.add<string>("fasta", 'f', "FASTA input file name", true, "");
    cmd.add<string>("outdir", 'o', "Directory for output. Default . means current directory.", false, ".");
    cmd.add<int>("kmer", 'k', "The length k of KMER (10~32), default 25", false, 25); 

    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }

    Options opt;

    opt.faFile = cmd.get<string>("fasta");
    opt.kmerKeyLen = cmd.get<int>("kmer");
    opt.outdir = cmd.get<string>("outdir");

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    opt.validate();


    Processor p(&opt);
    p.process();
    
    time_t t2 = time(NULL);

    cerr << endl << "Output to: " << opt.outdir << endl;
    cerr << endl << command << endl;
    cerr << "uniquekmer v" << UNIQUEKMER_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
