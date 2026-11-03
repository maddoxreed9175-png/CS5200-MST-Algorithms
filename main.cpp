#include <iostream>

#include "test.h"
#include "eval.h"

using namespace std;

int main(){

    
    mst_eval::runBenchmarks(cout, "bench.csv");
    
    return mst_test::runAllTests() ? 0 : 1;

    return 0;
}
