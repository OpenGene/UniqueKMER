#include "processor.h"
#include "kmer.h"
#include "genomes.h"

Processor::Processor(Options* opt) {
    mOptions = opt;
}

Processor::~Processor() {
}

bool Processor::process() {
	Genomes genomes(mOptions->faFile, mOptions);
	genomes.run();
	genomes.report();
    return true;
}