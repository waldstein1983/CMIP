//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_MODEL_H
#define CMIP_MODEL_H

#include <xprb.h>
#include <map>

using namespace std;

map<int, XPRBprob> subSolvers;
XPRBprob masterSolver = XPRBnewprob("master");
map<int, map<int, XPRBctr>> subBoundingCtrs;

//facility -> customer -> dual
map<int, map<int, double>> boundingVarSubDuals;

struct Solution {
    map<int, int> selectedLocations;
    double totalCost;
};


#endif //CMIP_MODEL_H
