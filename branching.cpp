//
// Created by baohuaw on 7/28/17.
//

#include "branching.h"


void buildBranchingConstraintSet(BranchingConstraintSet &targetSet, int branchingLocationId, int ctrType) {
    BranchingConstraintSet branchingSet = {nullptr, {}};
    branchingSet.basis = XPRBsavebasis(masterSolver);
//    nodeNum++;

    for (auto &ctr : targetSet.branchingCons) {
        branchingSet.branchingCons.push_back(ctr);
    }

    int bound = (int) XPRBgetsol(masterLocations[branchingLocationId]);
    if (ctrType == 0) {
        string ctrName = "y_";
        ctrName += to_string(branchingLocationId);
        ctrName += " <= ";
        ctrName += to_string(bound);
        BranchingConstraint left = {ctrName, branchingLocationId, 0, bound};
        branchingSet.branchingCons.push_back(left);
    } else {
        string ctrName = "y_";
        ctrName += to_string(branchingLocationId);
        ctrName += " >= ";
        ctrName += to_string(bound + 1);
        BranchingConstraint right = {ctrName, branchingLocationId, 1, bound + 1};
        branchingSet.branchingCons.push_back(right);
    }

    branchingNodes.insert(branchingNodes.begin(), branchingSet);

}


void branching(BranchingConstraintSet &set) {
    int targetBranchingLocationId = 0;
    double gapToHalf = DBL_MAX;
    for (int i = 1; i <= numFacility; i++) {
        if (abs(XPRBgetsol(masterLocations[i]) - round(XPRBgetsol(masterLocations[i]))) <= INT_GAP) {
            continue;
        }
        double fractional = XPRBgetsol(masterLocations[i]) - (int) XPRBgetsol(masterLocations[i]);
        if (abs(fractional - 0.5) < gapToHalf) {
            gapToHalf = abs(fractional - 0.5);
            targetBranchingLocationId = i;
        }
    }

    if (targetBranchingLocationId == 0) {
        return;
    }

    buildBranchingConstraintSet(set, targetBranchingLocationId, 0);
    buildBranchingConstraintSet(set, targetBranchingLocationId, 1);
}
