#include <iostream>
#include <chrono>
#include <xprb.h>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <cfloat>

#define MAX DBL_MAX
#define BD_GAP 1

using namespace std;
using namespace std::chrono;

int nodeNum = 0;

#define INT_GAP 0.00001

double LAMDA = 0.2;
double DELTA = 2 * INT_GAP;

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
    int bound;
};

struct BranchingConstraintSet {
    XPRBbasis basis;
    vector<BranchingConstraint> branchingCons;
};

struct LocalBranchingConstraintSet {
    XPRBbasis basis;
    vector<LocalBranchingConstraint> branchingCons;
};

milliseconds start;

vector<BranchingConstraintSet> branchingNodes;

vector<LocalBranchingConstraintSet> localBranchingNodes;

int numFacility;
int numCustomer;
map<int, XPRBprob> subSolvers;
XPRBprob masterSolver = XPRBnewprob("master");
map<int, XPRBvar> masterLocations;
map<int, XPRBvar> masterAlphas;
map<int, map<int, XPRBvar>> subCovers;

map<int, map<int, XPRBctr>> subBoundingCtrs;

double ub = MAX;

map<int, double> openingCosts;
map<int, map<int, double>> servingCosts;

int globalBendersCutId = 1;
//map<int, vector<int>> customerCriticals;

map<int, double> yy;

//facility -> customer -> dual
map<int, map<int, double>> boundingVarSubDuals;

struct Solution {
    map<int, int> selectedLocations;
    double totalCost;
};

Solution best = {{}, 0};

void print();

Solution roundingHeuristic2();

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

void readFromFile(const string &fileName) {
    ifstream in(fileName);
    if (!in) {
        cout << "Cannot open input file.\n";
        return;
    }

    string str;
    int lineId = 0;
    while (getline(in, str)) {
        istringstream iss(str);
        if (lineId == 0) {
            cout << str << endl;
        } else if (lineId == 1) {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));

            numFacility = stoi(tokens[0]);
            numCustomer = stoi(tokens[1]);
        } else {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));

            int facilityId = stoi(tokens[0]);
            double openingCost = stof(tokens[1]);
            openingCosts[facilityId] = openingCost;

            for (int j = 1; j <= tokens.size() - 2; j++) {
                servingCosts[facilityId][j] = stof(tokens[j - 1 + 2]);
            }
        }

        lineId++;
    }
    in.close();
}

void initMaster() {
    for (int i = 1; i <= numFacility; i++) {
        masterLocations[i] = XPRBnewvar(masterSolver, XPRB_UI, XPRBnewname("y_%d", i), 0, 1);
    }

    for (int j = 1; j <= numCustomer; j++) {
        masterAlphas[j] = XPRBnewvar(masterSolver, XPRB_PL, XPRBnewname("alpha_%d", j), 0, MAX);
    }

    XPRBctr obj = XPRBnewctr(masterSolver, "Obj", XPRB_N);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(obj, masterLocations[i], openingCosts[i]);
    }

    for (int j = 1; j <= numCustomer; j++) {
        XPRBaddterm(obj, masterAlphas[j], 1.0);
    }
    XPRBsetobj(masterSolver, obj);

    XPRBctr facilityExistence = XPRBnewctr(masterSolver, "Facility_existence", XPRB_G);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(facilityExistence, masterLocations[i], 1.0);
    }

    XPRBaddterm(facilityExistence, NULL, 1.0);
    XPRBsetsense(masterSolver, XPRB_MINIM);

}

void initSubModel() {
    for (int j = 1; j <= numCustomer; j++) {
        XPRBprob customer = XPRBnewprob("Sub");
        XPRBsetmsglevel(customer, 1);
        for (int i = 1; i <= numFacility; i++) {
            subCovers[i][j] = XPRBnewvar(customer, XPRB_PL, XPRBnewname("x_%d_%d", i, j), 0, MAX);
        }

        XPRBctr obj = XPRBnewctr(customer, XPRBnewname("Obj of Sub problem %d", j), XPRB_N);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(obj, subCovers[i][j], servingCosts[i][j]);
        }
        XPRBsetobj(customer, obj);

        XPRBctr fulfill = XPRBnewctr(customer, XPRBnewname("Fullfill of Sub problem %d", j), XPRB_E);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(fulfill, subCovers[i][j], 1.0);
        }
        XPRBaddterm(fulfill, nullptr, 1.0);

        for (int i = 1; i <= numFacility; i++) {
            XPRBctr ctr = XPRBnewctr(customer, XPRBnewname("Bounding with facility %d", i), XPRB_L);
            XPRBaddterm(ctr, subCovers[i][j], 1);
            XPRBaddterm(ctr, nullptr, 0);
            subBoundingCtrs[j][i] = ctr;
        }

        XPRBsetsense(customer, XPRB_MINIM);
        subSolvers[j] = customer;
    }
}

bool solveSubModel() {
    for (auto it = subSolvers.begin(); it != subSolvers.end(); ++it) {
        for (int i = 1; i <= numFacility; i++) {
            XPRBsetrange(subBoundingCtrs[it->first][i], -MAX, XPRBgetsol(masterLocations[i]));
        }

        XPRBlpoptimise(subSolvers[it->first], "");
        for (int i = 1; i <= numFacility; i++) {
            string ctrName = "Bounding with y_" + i;
            boundingVarSubDuals[i][it->first] = XPRBgetdual(subBoundingCtrs[it->first][i]);
        }
    }
    return true;
}

bool solveSubModel(Solution &solution) {
    for (auto it = subSolvers.begin(); it != subSolvers.end(); ++it) {
        for (int i = 1; i <= numFacility; i++) {
            XPRBsetrange(subBoundingCtrs[it->first][i], -MAX, solution.selectedLocations[i]);
        }

        XPRBlpoptimise(subSolvers[it->first], "");
        for (int i = 1; i <= numFacility; i++) {
            string ctrName = "Bounding with y_" + i;
            boundingVarSubDuals[i][it->first] = XPRBgetdual(subBoundingCtrs[it->first][i]);
        }
    }
    return true;
}

bool solveSubModelWithSeprator(map<int, double> &seprator) {
    for (auto it = subSolvers.begin(); it != subSolvers.end(); ++it) {
        for (int i = 1; i <= numFacility; i++) {
            XPRBsetrange(subBoundingCtrs[it->first][i], -MAX, seprator[i]);
        }

        XPRBlpoptimise(subSolvers[it->first], "");
        for (int i = 1; i <= numFacility; i++) {
            string ctrName = "Bounding with y_" + i;
            boundingVarSubDuals[i][it->first] = XPRBgetdual(subBoundingCtrs[it->first][i]);
        }
    }
    return true;
}

void print() {
    for (int i = 1; i <= numFacility; i++) {
        cout << "Facility " << i << " -> " << XPRBgetsol(masterLocations[i]) << endl;
        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {

        }
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            if (abs(XPRBgetsol(subCovers[i][j]) - 1) <= INT_GAP) {
                cout << "Customer " << j << " -> " << "Facility " << i << endl;
            }
        }
    }
}

void print(Solution solution) {
    for (int i = 1; i <= numFacility; i++) {
        cout << "Facility " << i << " -> " << solution.selectedLocations[i] << endl;
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            if (abs(XPRBgetsol(subCovers[i][j]) - 1) <= INT_GAP) {
                cout << "Customer " << j << " -> " << "Facility " << i << endl;
            }
        }
    }
}

double computeTotalCost() {
    double currentUb = 0;

    for (int i = 1; i <= numFacility; i++) {
        currentUb += XPRBgetsol(masterLocations[i]) * openingCosts[i];
//        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
//
//        }
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            currentUb += XPRBgetsol(subCovers[i][j]) * servingCosts[i][j];
//            if (abs(XPRBgetsol(subCovers[i][j]) - 1) <= INT_GAP) {
//                currentUb += XPRBgetsol(subCovers[i][j]) * servingCosts[i][j];
//            }
        }
    }
    return currentUb;
}
//
//void updateUB() {
//    double currentUb = 0;
//
//    for (int i = 1; i <= numFacility; i++) {
//        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
//            currentUb += XPRBgetsol(masterLocations[i]) * openingCosts[i];
//        }
//    }
//
//    for (int i = 1; i <= numFacility; i++) {
//        for (int j = 1; j <= numCustomer; j++) {
//            if (abs(XPRBgetsol(subCovers[i][j]) - 1) <= INT_GAP) {
//                currentUb += XPRBgetsol(subCovers[i][j]) * servingCosts[i][j];
//            }
//        }
//    }
//
////    print();
//    if (currentUb < ub)
//        ub = currentUb;
//}

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

void destroy() {
    XPRBdelprob(masterSolver);
    for (auto &subSolver : subSolvers) {
        XPRBdelprob(subSolver.second);
    }
    XPRBfinish();

}

bool isMasterSolutionInteger() {
    for (int i = 1; i <= numFacility; i++) {
//        if (abs(XPRBgetsol(masterLocations[i]) - 0) <= INT_GAP)
        if (abs(XPRBgetsol(masterLocations[i]) - round(XPRBgetsol(masterLocations[i]))) <= INT_GAP)
            continue;
        return false;
    }
    return true;
}

void bendersDecomposition(bool useOptimalityCut, bool useDualSimplex, bool useRoundingHeuristic) {
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
                best.totalCost = ub;
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
                solution = roundingHeuristic2();
                if (solution.totalCost < ub) {
                    ub = solution.totalCost;
                    solveSubModel(solution);
                    useHeuristicSolutionToAddCut = true;

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
        if (abs(nodeUB - nodeLB) <= 1)
            break;

        if (useHeuristicSolutionToAddCut) {
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
}

void bendersDecompositionWithInOut(bool useDualSimplex, bool useRoundingHeuristic) {
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
                solution = roundingHeuristic2();
                if (solution.totalCost < ub) {
                    ub = solution.totalCost;
                    solveSubModel(solution);
                    useHeuristicSolutionToAddCut = true;

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
}

void buildLocalBranchingConstraintSet(LocalBranchingConstraintSet &targetSet, int ctrType, int k) {
    LocalBranchingConstraintSet branchingSet = {nullptr, {}};

    for (auto &ctr : targetSet.branchingCons) {
        branchingSet.branchingCons.push_back(ctr);
    }

    vector<int> nonZeros;
    vector<int> zeros;
    for (int i = 1; i <= numFacility; i++) {
        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
            nonZeros.push_back(i);
        } else {
            zeros.push_back(i);
        }
    }

    if (ctrType == 0) {
        LocalBranchingConstraint left = {"Left", nonZeros, zeros, 0, k};
        branchingSet.branchingCons.push_back(left);
    } else {
        LocalBranchingConstraint right = {"Right", nonZeros, zeros, 0, k + 1};
        branchingSet.branchingCons.push_back(right);
    }
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

void localBranching(LocalBranchingConstraintSet &set, int k) {
    buildLocalBranchingConstraintSet(set, 0, k);
    buildLocalBranchingConstraintSet(set, 1, k);
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

//Solution roundingHeuristic() {
//    Solution solution = {{}, 0};
//    int tempLocation[numFacility + 1];
//    for (int i = 1; i <= numFacility; i++) {
//        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
//            tempLocation[i] = 1;
//            solution.totalCost += openingCosts[i];
//        } else {
//            if (XPRBgetsol(masterLocations[i]) <= 0.5) {
//                tempLocation[i] = 0;
//            } else {
//                solution.totalCost += openingCosts[i];
//                tempLocation[i] = 1;
//            }
//        }
//    }
//
//    bool locationSelected[numFacility + 1];
//
//    for (int i = 1; i <= numFacility; i++) {
//        locationSelected[i] = false;
//    }
//
//    for (int j = 1; j <= numCustomer; j++) {
//        int neareastFacility = 0;
//        double nearestDistance = MAX;
//        for (int i = 1; i <= numFacility; i++) {
//            if (tempLocation[i] == 1) {
//                if (servingCosts[i][j] < nearestDistance) {
//                    nearestDistance = servingCosts[i][j];
//                    neareastFacility = i;
//                }
//            }
//        }
//        solution.totalCost += servingCosts[neareastFacility][j];
//        locationSelected[neareastFacility] = true;
//    }
//
//    for (int i = 1; i <= numFacility; i++) {
//        solution.selectedLocations[i] = 0;
//        if (!locationSelected[i] && tempLocation[i] == 1) {
//            solution.totalCost -= openingCosts[i];
//        } else if (locationSelected[i] && tempLocation[i] == 1) {
//            solution.selectedLocations[i] = 1;
//        }
//    }
//    return solution;
//}

Solution roundingHeuristic2() {
    double sortedLocation[numFacility + 1];
    for (int i = 1; i <= numFacility; i++) {
        sortedLocation[i] = XPRBgetsol(masterLocations[i]);
    }

    sort(sortedLocation + 1, sortedLocation + numFacility + 1);

    vector<Solution> solutions;
    double curThreshold = -1;
    for (int i = 1; i <= numFacility; i++) {
        if (sortedLocation[i] != curThreshold) {
            Solution solution = {{}, 0};
            curThreshold = sortedLocation[i];
            int roundingLocation[numFacility + 1];
            for (int ii = 1; ii <= numFacility; ii++) {
                if (abs(XPRBgetsol(masterLocations[ii]) - 1) <= INT_GAP) {
                    roundingLocation[ii] = 1;
                    solution.totalCost += openingCosts[ii];
                } else {
                    if (XPRBgetsol(masterLocations[ii]) < curThreshold) {
                        roundingLocation[ii] = 0;
                    } else {
                        solution.totalCost += openingCosts[ii];
                        roundingLocation[ii] = 1;
                    }
                }
            }

            bool locationSelected[numFacility + 1];

            for (int ii = 1; ii <= numFacility; ii++) {
                locationSelected[ii] = false;
            }

            for (int j = 1; j <= numCustomer; j++) {
                int neareastFacility = 0;
                double nearestDistance = MAX;
                for (int ii = 1; ii <= numFacility; ii++) {
                    if (roundingLocation[ii] == 1) {
                        if (servingCosts[ii][j] < nearestDistance) {
                            nearestDistance = servingCosts[ii][j];
                            neareastFacility = ii;
                        }
                    }
                }
                solution.totalCost += servingCosts[neareastFacility][j];
                locationSelected[neareastFacility] = true;
            }

            for (int ii = 1; ii <= numFacility; ii++) {
                solution.selectedLocations[ii] = 0;
                if (!locationSelected[ii] && roundingLocation[ii] == 1) {
                    solution.totalCost -= openingCosts[ii];
                } else if (locationSelected[ii] && roundingLocation[ii] == 1) {
                    solution.selectedLocations[ii] = 1;
                }
            }
//            cout << i << "    " << solution.totalCost << endl;
            solutions.push_back(solution);
        }
    }


    sort(solutions.begin(), solutions.begin() + solutions.size(),
         [](Solution const &a, Solution const &b) { return a.totalCost < b.totalCost; });

    return solutions[0];
}

void branchAndCut(bool useOptimalityCutInBendersDecomposition, bool useDualSimplexInBendersDecomposition,
                  bool useRoundingHeuristicInBendersDecomposition,
                  bool useDualSimplexInBranch, bool useInOutStrategyInBendersDecomposition) {
    XPRBsetmsglevel(masterSolver, 1);
    initMaster();
    initSubModel();

    XPRBlpoptimise(masterSolver, "");
    if (!useInOutStrategyInBendersDecomposition) {
        bendersDecomposition(useOptimalityCutInBendersDecomposition, useDualSimplexInBendersDecomposition,
                             useRoundingHeuristicInBendersDecomposition);
    } else {
        bendersDecompositionWithInOut(useDualSimplexInBendersDecomposition, useRoundingHeuristicInBendersDecomposition);
    }


    if (isMasterSolutionInteger()) {
        cout << "UB = " << ub << " at the root node " << endl;
        print();
        return;
    }

    BranchingConstraintSet target = {nodeNum, {}};
    branching(target);

    int step = 1;
    while (!branchingNodes.empty()) {
        milliseconds end = duration_cast<milliseconds>(
                system_clock::now().time_since_epoch()
        );
        if (step % 1 == 0) {
            cout << ub << "   node size " << branchingNodes.size() << " step = " << step << "   Time "
                 << (end.count() - start.count()) << endl;

        }
        if ((end.count() - start.count()) / 1000 > 720000) {
            cout << "Terminate due to time limit of 3600 sec" << endl;
            break;
        }
        target = branchingNodes[0];
        branchingNodes.erase(branchingNodes.begin());

        for (auto &branching : target.branchingCons) {
            if (branching.ctrType == 0) {
                XPRBctr ctr = XPRBnewctr(masterSolver, branching.name.c_str(), XPRB_L);
                XPRBaddterm(ctr, masterLocations[branching.branchingLocationId], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);

            } else if (branching.ctrType == 1) {
                XPRBctr ctr = XPRBnewctr(masterSolver, branching.name.c_str(), XPRB_G);
                XPRBaddterm(ctr, masterLocations[branching.branchingLocationId], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);
            }
        }

        if (useDualSimplexInBranch) {
            XPRBloadmat(masterSolver);
            XPRBloadbasis(target.basis);
//            basis = nullptr;
            XPRBlpoptimise(masterSolver, "d");
//            basis = XPRBsavebasis(masterSolver);
        } else {
            XPRBlpoptimise(masterSolver, "");
        }


        if (XPRBgetlpstat(masterSolver) == XPRB_LP_OPTIMAL) {
            if (XPRBgetobjval(masterSolver) < ub && abs(XPRBgetobjval(masterSolver) - ub) >= INT_GAP) {
                bendersDecomposition(useOptimalityCutInBendersDecomposition, useDualSimplexInBendersDecomposition,
                                     useRoundingHeuristicInBendersDecomposition);
                branching(target);
            }
        }

        for (auto &branching : target.branchingCons) {
            auto branchingCtr = (XPRBctr) XPRBgetbyname(masterSolver, branching.name.c_str(), XPRB_CTR);
            XPRBdelctr(branchingCtr);
        }
        step++;
    }
    cout << "UB = " << ub << endl;

    cout << "Total Cost " << best.totalCost << endl;
    for (int i = 1; i <= numFacility; i++) {
        cout << "facility " << i << " -> " << best.selectedLocations[i] << endl;
    }
    for (int j = 1; j <= numCustomer; j++) {
        int neareastFacility = 0;
        double nearestDistance = MAX;
        for (int i = 1; i < numFacility; i++) {
            if (abs(best.selectedLocations[i] - 1) <= INT_GAP) {
                if (servingCosts[i][j] < nearestDistance) {
                    nearestDistance = servingCosts[i][j];
                    neareastFacility = i;
                }
            }
        }
        cout << "Customer " << j << " -> " << neareastFacility << endl;
//        solution.totalCost += servingCosts[neareastFacility][j];
//        locationSelected[neareastFacility] = true;
    }

}

int main() {
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/simpleExample2.txt";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/50/50.1";
//        string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/70/70.1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/150/150.5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/200/200.10";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/a/gs250a-5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/b/gs250b-3";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/500/a/gs500a-2";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/750/a/gs750a-1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/kmedian/1000-10";
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/kmedian/2500-10";
    readFromFile(fileName);
    start = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );
//    bendersDecomposition();

//    branchAndCut(false,true);
    branchAndCut(true, true, true, false, true);
//    branchAndCutWithInOut(false);
//    branchAndCutWithDualSimplex(false);
    milliseconds end = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    cout << "Time " << (end.count() - start.count()) << endl;
    destroy();
    return 0;
}