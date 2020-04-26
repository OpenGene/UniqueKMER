#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string.h>

Options::Options(){
    faFile = "";
    hgFile = "";
    outdir = "";
    kmerKeyLen = 25;
    genomeLimit = 0;
    kmerSpacing = 0;
}

void Options::init() {
}

bool Options::validate() {
    if(faFile.empty()) {
        error_exit("FASTA file (-f) should be specified");
    } else {
        check_file_valid(faFile);
    }

    if(!faFile.empty()) {
        check_file_valid(faFile);
    }

    if(!hgFile.empty()) {
        check_file_valid(hgFile);
    }

    if(!outdir.empty()) {
        if(!file_exists(outdir)) {
            if(mkdir(outdir.c_str(), 0755)!=0)
                error_exit("Failed to create directory: " + outdir);
        }

        if(!is_directory(outdir))
            error_exit("Not a directory: " + outdir);
    }

    if(kmerKeyLen < 10 || kmerKeyLen > 32)
        error_exit("KMER length (-k) should be between 10 ~ 32, suggest 25");

    if(genomeLimit < 0)
        error_exit("Genome limit (-g) cannot be negative");

    return true;
}