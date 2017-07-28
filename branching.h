//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BRANCHING_H
#define CMIP_BRANCHING_H

//#include <xprb.h>
#include <vector>
#include <cfloat>
#include <cmath>
#include "ufl_model.h"

using namespace std;

int LOCAL_BRANCHING_K = 10;

//int nodeNum = 0;


struct BranchingConstraint {
    string name;
    int branchingLocationId;
    int ctrType;
    int bound;
};


struct BranchingConstraintSet {
    XPRBbasis basis;
    vector<BranchingConstraint> branchingCons;
};

vector<BranchingConstraintSet> branchingNodes;

void buildBranchingConstraintSet(BranchingConstraintSet &targetSet, int branchingLocationId, int ctrType);

void branching(BranchingConstraintSet &set);

#endif //CMIP_BRANCHING_H
