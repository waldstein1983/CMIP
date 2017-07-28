//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_LOCAL_BRANCHING_H
#define CMIP_LOCAL_BRANCHING_H

#include <string>
#include <vector>
//#include <xprb.h>
#include "ufl_model.h"

using namespace std;

int nodeNum = 0;

struct LocalBranchingConstraint {
    string name;
    vector<int> nonZeros;
    vector<int> zeros;
    int ctrType;
    int k;
};


struct LocalBranchingConstraintSet {
    XPRBbasis basis;
    vector<LocalBranchingConstraint> branchingCons;
};

vector<LocalBranchingConstraintSet> localBranchingNodes;

void buildLocalBranchingConstraintSet(Solution &best, LocalBranchingConstraintSet &targetSet, int ctrType, int k);

void localBranching(Solution &best, LocalBranchingConstraintSet &set, int k);

#endif //CMIP_LOCAL_BRANCHING_H
