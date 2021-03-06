#include "unittest.h"
#include "sequence.h"
#include "kmer.h"
#include <time.h>

UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= Kmer::test();
    printf("\n==========================\n");
    printf("%s\n\n", passed?"ALL PASSED":"FAILED");
}

bool UnitTest::report(bool result, string message) {
    printf("%s:%s\n\n", message.c_str(), result?" PASSED":" FAILED");
    return result;
}