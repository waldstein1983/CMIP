//
// Created by baohuaw on 7/28/17.
//


#include "benders_cut.h"



bool isMasterSolutionInteger() {
    for (int i = 1; i <= numFacility; i++) {
        if (abs(XPRBgetsol(masterLocations[i]) - round(XPRBgetsol(masterLocations[i]))) <= INT_GAP)
            continue;
        return false;
    }
    return true;
}

void addOptimalityCut() {
    XPRBctr cut = XPRBnewctr(masterSolver, "Optimality cut", XPRB_G);
    double sumServingCost = 0;
    for (int j = 1; j <= numCustomer; j++) {
        sumServingCost += XPRBgetobjval(subSolvers[j]);
        XPRBaddterm(cut, masterAlphas[j], 1);

    }
    for (int i = 1; i < numFacility; i++) {
        if (abs(XPRBgetsol(masterLocations[i]) - 0) <= INT_GAP) {
            XPRBaddterm(cut, masterLocations[i], sumServingCost);
        }
    }

    XPRBaddterm(cut, nullptr, sumServingCost);
}

bool addBendersCutForEachSubProblemToMaster() {
    bool newCut = false;
    for (int j = 1; j <= numCustomer; j++) {
        double sumDual = 0;
        for (int i = 1; i <= numFacility; i++) {
//            cout << boundingVarSubDuals[i][j] << " * " << XPRBgetsol(masterLocations[i]) << endl;
            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
        }

        XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("benders cut %d", globalBendersCutId), XPRB_G);

        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(cut, masterLocations[i], -boundingVarSubDuals[i][j]);
//            cout <<  XPRBgetsol(masterLocations[i]) << endl;
//            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
        }
        XPRBaddterm(cut, masterAlphas[j], 1);
        XPRBaddterm(cut, nullptr, XPRBgetobjval(subSolvers[j]) - sumDual);
//        XPRBprintctr(cut);
        newCut = true;
    }
    return newCut;
}


Solution bendersDecomposition(bool useOptimalityCut, bool useDualSimplex, bool useRoundingHeuristic) {
    Solution best = {{}, MAX};
    XPRBbasis basis = XPRBsavebasis(masterSolver);
    double nodeLB = 0, nodeUB = MAX;
    while (abs(nodeUB - nodeLB) > 1) {
        if (XPRBgetobjval(masterSolver) > ub) {
            break;
        }
        nodeLB = XPRBgetobjval(masterSolver);
        solveSubModel();
        nodeUB = computeTotalCost();

        bool useHeuristicSolutionToAddCut = false;
        Solution solution;
        if (isMasterSolutionInteger()) {
            if (useOptimalityCut)
                addOptimalityCut();
            if (nodeUB < ub) {
                ub = nodeUB;

            }

            if(nodeUB < best.totalCost){
                best.totalCost = nodeUB;
                for (int i = 1; i <= numFacility; i++) {
                    if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
                        best.selectedLocations[i] = 1;
                    } else {
                        best.selectedLocations[i] = 0;
                    }
                }
            }
        } else {
            if (useRoundingHeuristic) {
                solution = roundingHeuristic();
                if (solution.totalCost < ub) {
                    ub = solution.totalCost;

//                    useHeuristicSolutionToAddCut = true;

//                    best.totalCost = ub;
//                    for (int i = 1; i <= numFacility; i++) {
//                        if (solution.selectedLocations[i] == 1) {
//                            best.selectedLocations[i] = 1;
//                        } else {
//                            best.selectedLocations[i] = 0;
//                        }
//                    }
                }

                if(solution.totalCost < best.totalCost){
                    best = solution;
//                    best.totalCost = solution.totalCost;
//                    for (int i = 1; i <= numFacility; i++) {
//                        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
//                            best.selectedLocations[i] = 1;
//                        } else {
//                            best.selectedLocations[i] = 0;
//                        }
//                    }
                }
            }
        }
        if (abs(nodeUB - nodeLB) <= 1)
            break;

        if (useHeuristicSolutionToAddCut) {
//            solveSubModel(solution);
//            addBendersCutForEachSubProblemToMaster(solution);
        } else {
            addBendersCutForEachSubProblemToMaster();
        }
        if (useDualSimplex) {
            XPRBloadmat(masterSolver);
            XPRBloadbasis(basis);
            basis = nullptr;
            XPRBlpoptimise(masterSolver, "d");
            basis = XPRBsavebasis(masterSolver);
        } else
            XPRBlpoptimise(masterSolver, "");
    }
    return best;
}

