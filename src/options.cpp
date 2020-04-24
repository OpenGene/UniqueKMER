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

    if(kmerKeyLen < 10 || kmerKeyLen > 32)
        error_exit("KMER length (-k) should be between 10 ~ 32, suggest 25");

    return true;
}