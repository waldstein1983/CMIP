//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BENDERS_CUT_H
#define CMIP_BENDERS_CUT_H

//#define MAX DBL_MAX

#include "problem.h"
#include "model.h"
#include <cmath>
#include <cfloat>
#include "rounding_heuristic.h"

int globalBendersCutId = 1;


double LAMDA = 0.2;
double DELTA = 2 * INT_GAP;

map<int, double> yy;


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

bool addBendersCutForEachSubProblemToMaster(Solution &solution) {
    bool newCut = false;
    for (int j = 1; j <= numCustomer; j++) {
        double sumDual = 0;
        for (int i = 1; i <= numFacility; i++) {
//            cout << boundingVarSubDuals[i][j] << " * " << XPRBgetsol(masterLocations[i]) << endl;
            sumDual += boundingVarSubDuals[i][j] * solution.selectedLocations[i];
        }

        XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("benders cut %d", globalBendersCutId), XPRB_G);

        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(cut, masterLocations[i], -boundingVarSubDuals[i][j]);
//            cout <<  XPRBgetsol(masterLocations[i]) << endl;
//            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
        }
        XPRBaddterm(cut, masterAlphas[j], 1);
        XPRBaddterm(cut, nullptr, XPRBgetobjval(subSolvers[j]) - sumDual);
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
            solveSubModel(solution);
            addBendersCutForEachSubProblemToMaster(solution);
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


map<int, double> stabilize() {
    for (int i = 1; i <= numFacility; i++) {
        double temp = (XPRBgetsol(masterLocations[i]) + yy[i]) * 0.5;
        yy[i] = temp;
    }
    map<int, double> separator;
    for (int i = 1; i <= numFacility; i++) {
        double temp = LAMDA * XPRBgetsol(masterLocations[i]) + (1 - LAMDA) * yy[i] + DELTA;
        separator[i] = temp;
    }
    return separator;
}

void addBendersCutWithSeparator(map<int, double> &separator) {
    for (int j = 1; j <= numCustomer; j++) {
        double sumDual = 0;
        for (int i = 1; i <= numFacility; i++) {
//            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
            sumDual += boundingVarSubDuals[i][j] * separator[i];
        }

        XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("benders cut %d", globalBendersCutId), XPRB_G);

        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(cut, masterLocations[i], -boundingVarSubDuals[i][j]);
        }
        XPRBaddterm(cut, masterAlphas[j], 1);
        XPRBaddterm(cut, nullptr, XPRBgetobjval(subSolvers[j]) - sumDual);
    }
}






Solution bendersDecompositionWithInOut(bool useDualSimplex, bool useRoundingHeuristic) {
    Solution best = {{}, 0};
    XPRBbasis basis = XPRBsavebasis(masterSolver);
    for (int i = 1; i <= numFacility; i++) {
        yy[i] = 1;
    }

    int LBnoImproveStep = 0;
    double nodeLB = 0, nodeUB = MAX;
//    while (abs(nodeUB - nodeLB) > 1) {
//    while ((nodeLB < nodeUB) || abs(nodeUB - nodeLB) > BD_GAP) {
    while (abs(nodeUB - nodeLB) > BD_GAP) {
        if (nodeLB <= nodeUB) {
            cout << nodeLB << "    " << nodeUB << endl;
        }
        if (XPRBgetobjval(masterSolver) > ub) {
            break;
        }

//        bool useHeuristicSolutionToAddCut = false;
//        Solution solution;
////        if (useRoundingHeuristic) {
////            solution = roundingHeuristic2();
////            if (solution.totalCost < ub) {
////                ub = solution.totalCost;
////                solveSubModel(solution);
////                useHeuristicSolutionToAddCut = true;
////
////                best.totalCost = ub;
////                for (int i = 1; i <= numFacility; i++) {
////                    if (solution.selectedLocations[i] == 1) {
////                        best.selectedLocations[i] = 1;
////                    } else {
////                        best.selectedLocations[i] = 0;
////                    }
////                }
////            }
////        }
//
////        solveSubModel();
////        nodeUB = computeTotalCost();
//        if (isMasterSolutionInteger()) {
//            solveSubModel();
//            double cost = computeTotalCost();
//            if (cost < ub){
//                ub = cost;
//                best.totalCost = ub;
//                for (int i = 1; i <= numFacility; i++) {
//                    if (solution.selectedLocations[i] == 1) {
//                        best.selectedLocations[i] = 1;
//                    } else {
//                        best.selectedLocations[i] = 0;
//                    }
//                }
//            }
//
//        }else{
//            if (useRoundingHeuristic) {
//                solution = roundingHeuristic2();
//                if (solution.totalCost < ub) {
//                    ub = solution.totalCost;
//                    solveSubModel(solution);
//                    useHeuristicSolutionToAddCut = true;
//
//                    best.totalCost = ub;
//                    for (int i = 1; i <= numFacility; i++) {
//                        if (solution.selectedLocations[i] == 1) {
//                            best.selectedLocations[i] = 1;
//                        } else {
//                            best.selectedLocations[i] = 0;
//                        }
//                    }
//                }
//            }
//        }
//
//        map<int, double> separator = stabilize();
//        solveSubModelWithSeprator(separator);
//        nodeUB = computeTotalCost();
//        addBendersCutWithSeparator(separator);


        /////////////////////////////////////////////////////////////////////////////


        bool useHeuristicSolutionToAddCut = false;
        Solution solution;
//        if (useRoundingHeuristic) {
//            solution = roundingHeuristic2();
//            if (solution.totalCost < ub) {
//                ub = solution.totalCost;
//                solveSubModel(solution);
//                useHeuristicSolutionToAddCut = true;
//
//                best.totalCost = ub;
//                for (int i = 1; i <= numFacility; i++) {
//                    if (solution.selectedLocations[i] == 1) {
//                        best.selectedLocations[i] = 1;
//                    } else {
//                        best.selectedLocations[i] = 0;
//                    }
//                }
//            }
//        }

        solveSubModel();
        nodeUB = computeTotalCost();
        if (isMasterSolutionInteger()) {
//            solveSubModel();
//            double cost = computeTotalCost();
            if (nodeUB < ub) {
                ub = nodeUB;
                best.totalCost = ub;
                for (int i = 1; i <= numFacility; i++) {
                    if (solution.selectedLocations[i] == 1) {
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
//                    solveSubModel(solution);
//                    useHeuristicSolutionToAddCut = true;

                    best.totalCost = ub;
                    for (int i = 1; i <= numFacility; i++) {
                        if (solution.selectedLocations[i] == 1) {
                            best.selectedLocations[i] = 1;
                        } else {
                            best.selectedLocations[i] = 0;
                        }
                    }
                }
            }
        }

        map<int, double> separator = stabilize();
        solveSubModelWithSeprator(separator);
//        nodeUB = computeTotalCost();
        addBendersCutWithSeparator(separator);


        ////////////////////////////////////////////////////////////////

//        if (useHeuristicSolutionToAddCut) {
//            addBendersCutForEachSubProblemToMaster(solution);
//        } else {
//            map<int, double> separator = stabilize();
//            solveSubModelWithSeprator(separator);
//            nodeUB = computeTotalCost();
//            addBendersCutWithSeparator(separator);
////            addBendersCutForEachSubProblemToMaster();
//        }

        if (abs(XPRBgetobjval(masterSolver) - nodeLB) <= 1) {
            LBnoImproveStep++;

            if (LBnoImproveStep == 5) {
                if (LAMDA != 1)
                    LAMDA = 1;
                else {
                    if (DELTA == 0.2) {
                        DELTA = 0;
                    } else {
                        break;
                    }
                }

                LBnoImproveStep = 0;
            }
        } else {
            LBnoImproveStep = 0;
        }
        nodeLB = XPRBgetobjval(masterSolver);


        if (abs(nodeUB - nodeLB) <= BD_GAP)
            break;


//        XPRBlpoptimise(masterSolver, "");
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

#endif //CMIP_BENDERS_CUT_H
