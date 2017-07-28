//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BRANCHING_H
#define CMIP_BRANCHING_H

#include <xprb.h>
#include <vector>

using namespace std;

int LOCAL_BRANCHING_K = 10;

struct BranchingConstraint {
    string name;
    int branchingLocationId;
    int ctrType;
    int bound;
};

struct LocalBranchingConstraint {
    string name;
//    int branchingLocationId;
    vector<int> nonZeros;
    vector<int> zeros;
    int ctrType;
    int k;
};

struct BranchingConstraintSet {
    XPRBbasis basis;
    vector<BranchingConstraint> branchingCons;
};

struct LocalBranchingConstraintSet {
    XPRBbasis basis;
    vector<LocalBranchingConstraint> branchingCons;
};

vector<BranchingConstraintSet> branchingNodes;

vector<LocalBranchingConstraintSet> localBranchingNodes;

#endif //CMIP_BRANCHING_H
