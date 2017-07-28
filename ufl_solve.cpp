#include <iostream>
#include <chrono>
#include <map>
#include <iterator>
#include <algorithm>
#include "branching.h"
#include "branch_cut.h"

using namespace std;
using namespace std::chrono;

milliseconds start;


int main() {
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/simpleExample2.txt";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/50/50.1";
//        string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/70/70.1";
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/150/150.5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/200/200.10";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/a/gs250a-5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/b/gs250b-3";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/500/a/gs500a-2";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/750/a/gs750a-1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/kmedian/1000-10";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/kmedian/2500-10";
    readFromFile(fileName);
    start = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );
//    bendersDecomposition();

//    branchAndCut(false,true);
    branchAndCut(true, true, true, false, false);
//    branchAndCutWithLocalBranching(true, true, true, false, false);
//    branchAndCutWithInOut(false);
//    branchAndCutWithDualSimplex(false);
    milliseconds end = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    cout << "Time " << (end.count() - start.count()) << endl;
    destroy();
    return 0;
}