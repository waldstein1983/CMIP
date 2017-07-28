//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BRANCHING_H
#define CMIP_BRANCHING_H

#include <xprb.h>
#include <vector>
#include "model.h"
#include "decision_var.h"

using namespace std;

int LOCAL_BRANCHING_K = 10;

int nodeNum = 0;

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

void buildBranchingConstraintSet(BranchingConstraintSet &targetSet, int branchingLocationId, int ctrType) {
    BranchingConstraintSet branchingSet = {nullptr, {}};
    branchingSet.basis = XPRBsavebasis(masterSolver);
    nodeNum++;

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

void localBranching(Solution &best, LocalBranchingConstraintSet &set, int k) {
    buildLocalBranchingConstraintSet(best, set, 0, k);
    buildLocalBranchingConstraintSet(best, set, 1, k);
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

#endif //CMIP_BRANCHING_H
