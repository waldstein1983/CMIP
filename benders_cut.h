//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BENDERS_CUT_H
#define CMIP_BENDERS_CUT_H

//#define MAX DBL_MAX

#include "ufl_model.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include "rounding_heu.h"

int globalBendersCutId = 1;



bool isMasterSolutionInteger();

void addOptimalityCut();

bool addBendersCutForEachSubProblemToMaster();

//bool addBendersCutForEachSubProblemToMaster(Solution &solution) {
//    bool newCut = false;
//    for (int j = 1; j <= numCustomer; j++) {
//        double sumDual = 0;
//        for (int i = 1; i <= numFacility; i++) {
////            cout << boundingVarSubDuals[i][j] << " * " << XPRBgetsol(masterLocations[i]) << endl;
//            sumDual += boundingVarSubDuals[i][j] * solution.selectedLocations[i];
//        }
//
//        XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("benders cut %d", globalBendersCutId), XPRB_G);
//
//        for (int i = 1; i <= numFacility; i++) {
//            XPRBaddterm(cut, masterLocations[i], -boundingVarSubDuals[i][j]);
////            cout <<  XPRBgetsol(masterLocations[i]) << endl;
////            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
//        }
//        XPRBaddterm(cut, masterAlphas[j], 1);
//        XPRBaddterm(cut, nullptr, XPRBgetobjval(subSolvers[j]) - sumDual);
//        newCut = true;
//    }
//    return newCut;
//}


Solution bendersDecomposition(bool useOptimalityCut, bool useDualSimplex, bool useRoundingHeuristic);




#endif //CMIP_BENDERS_CUT_H
