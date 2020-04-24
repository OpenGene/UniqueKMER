#include "processor.h"

Processor::Processor(Options* opt) {
    mOptions = opt;
}

Processor::~Processor() {
}

bool Processor::process() {
    return true;
}