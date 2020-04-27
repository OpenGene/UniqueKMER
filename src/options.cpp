#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string.h>

Options::Options(){
    faFile = "";
    refFile = "";
    outdir = "";
    kmerKeyLen = 25;
    genomeLimit = 0;
    kmerSpacing = 0;
    edThreshold = 3;
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

    if(!refFile.empty()) {
        check_file_valid(refFile);
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

    if(edThreshold < 0 || edThreshold > 16)
        error_exit("Edit distance threshold (-e) should be between 0 ~ 16, suggest 3");

    if(genomeLimit < 0)
        error_exit("Genome limit (-g) cannot be negative");

    return true;
}