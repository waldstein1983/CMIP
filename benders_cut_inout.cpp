//
// Created by baohuaw on 7/28/17.
//

#include "benders_cut_inout.h"
#include "benders_cut.h"

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