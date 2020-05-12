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
        cerr << "UniqueKMER: Generate unique k-mers for every contig in a FASTA file." << endl << "version " << UNIQUEKMER_VER << endl;
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
    cmd.add<string>("outdir", 'o', "Directory for output. Default is unique_kmers in the current directory.", false, "unique_kmers");
    cmd.add<int>("kmer", 'k', "The length k of k-mer (10~32), default 25", false, 25); 
    cmd.add<int>("spacing", 's', "If a key with POS is recorded, then skip [POS+1...POS+spacing] to avoid too compact result (0~100). default 0 means no skipping.", false, 0); 
    cmd.add<int>("genome_limit", 'g', "Process up to genome_limit genomes in the FASTA input file. Default 0 means no limit. This option is for DEBUG.", false, 0); 
    cmd.add<string>("ref", 'r', "Reference genome FASTA file name. Specify this only when you want to filter out the unique k-mer that can be mapped to reference genome.", false, "");
    cmd.add<int>("edit_distance", 'e', "k-mer mapped to reference genome with edit distance <= edit_distance will be removed (0~16). 3 for default.", false, 3);

    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }

    Options opt;

    opt.faFile = cmd.get<string>("fasta");
    opt.kmerKeyLen = cmd.get<int>("kmer");
    opt.genomeLimit = cmd.get<int>("genome_limit");
    opt.kmerSpacing = cmd.get<int>("spacing");
    opt.outdir = cmd.get<string>("outdir");
    opt.refFile = cmd.get<string>("ref");
    opt.edThreshold = cmd.get<int>("edit_distance");

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

    cerr << endl << "Please find results (index.html, k-mer/Genome files) in folder: " << opt.outdir << endl;
    cerr << endl << command << endl;
    cerr << "uniquekmer v" << UNIQUEKMER_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
