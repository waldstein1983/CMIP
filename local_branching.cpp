//
// Created by baohuaw on 7/28/17.
//

#include "local_branching.h"

void buildLocalBranchingConstraintSet(Solution &best, LocalBranchingConstraintSet &targetSet, int ctrType, int k) {
    LocalBranchingConstraintSet branchingSet = {nullptr, {}};

    for (auto &ctr : targetSet.branchingCons) {
        branchingSet.branchingCons.push_back(ctr);
    }

    vector<int> nonZeros;
    vector<int> zeros;
    for (int i = 1; i <= numFacility; i++) {
        if (best.selectedLocations[i] == 1) {
            nonZeros.push_back(i);
        } else {
            zeros.push_back(i);
        }
    }

    if (ctrType == 0) {
        LocalBranchingConstraint left = {"Left " + nodeNum, nonZeros, zeros, 0, k};
        nodeNum++;
        branchingSet.branchingCons.push_back(left);
    } else {
        LocalBranchingConstraint right = {"Right " + nodeNum, nonZeros, zeros, 0, k + 1};
        nodeNum++;
        branchingSet.branchingCons.push_back(right);
    }
    localBranchingNodes.insert(localBranchingNodes.begin(), branchingSet);
}


void localBranching(Solution &best, LocalBranchingConstraintSet &set, int k) {
    buildLocalBranchingConstraintSet(best, set, 0, k);
    buildLocalBranchingConstraintSet(best, set, 1, k);
}


