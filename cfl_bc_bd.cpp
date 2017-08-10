//
// Created by baohuaw on 8/1/17.
// This implementation should ref. Fesichetti's paper, where xij is the ratio that demand j is assigned to location i
//

#include <iostream>
#include <chrono>
#include <map>
#include <iterator>
#include <algorithm>
#include <cfloat>
#include <fstream>
#include <sstream>

#include <xprb.h>

#define BD_TERMINATE_GAP 100
#define INT_GAP 0.00001
#define MAX DBL_MAX

using namespace std;
using namespace std::chrono;

int numFacility;
int numCustomer;

double ub = MAX;

int globalBendersCutId = 1;

map<int, double> openingCosts;
map<int, double> capacities;
map<int, map<int, double>> servingCosts;

map<int, double> demands;


map<int, XPRBvar> y;
XPRBvar masterOmega;
map<int, map<int, XPRBvar>> x;

XPRBprob originalSolver = XPRBnewprob("original");

XPRBprob subSolver = XPRBnewprob("sub");
map<int, map<int, XPRBctr>> boundingCtrs;
map<int, XPRBctr> capacityCtrs;
XPRBprob masterSolver = XPRBnewprob("master");
XPRBprob kp = XPRBnewprob("knapSack");
XPRBctr kpObj;
XPRBctr kpCtr;
map<int, XPRBvar> kpVars;
map<int, XPRBctr> fracCnsrvtnCtrs;
map<int, double> fracCnsrvtnDuals;

struct Solution {
    map<int, int> selectedLocations;
    double totalCost;
};

////////////////////////////////////////////
//branching
///////////////////////////////////////////
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

void printMasterVar();

vector<BranchingConstraintSet> branchingNodes;

/////////////////////////////////////////////////////
//Build master model
void initMaster() {
    for (int i = 1; i <= numFacility; i++) {
        y[i] = XPRBnewvar(masterSolver, XPRB_UI, XPRBnewname("y_%d", i), 0, 1);
    }

    masterOmega = XPRBnewvar(masterSolver, XPRB_PL, XPRBnewname("omega"), 0, MAX);


    XPRBctr obj = XPRBnewctr(masterSolver, "Obj", XPRB_N);
//    for (int i = 1; i <= numFacility; i++) {
//        XPRBaddterm(obj, y[i], openingCosts[i]);
//    }
    XPRBaddterm(obj, masterOmega, 1.0);

    XPRBsetobj(masterSolver, obj);
//    XPRBprintctr(obj);

    XPRBctr totalCapacity = XPRBnewctr(masterSolver, "Sum Capacity >= Sum Demand", XPRB_G);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(totalCapacity, y[i], capacities[i]);
    }


    int totalDemand = 0;
    for (int j = 1; j <= numCustomer; j++) {
        totalDemand += demands[j];
    }
    XPRBaddterm(totalCapacity, NULL, totalDemand);

//    XPRBprintctr(totalCapacity);

    XPRBsetsense(masterSolver, XPRB_MINIM);

}


void printProblem() {
    for (int i = 1; i <= numFacility; i++) {
        cout << "facility " << i << " opening cost " << openingCosts[i] << "  capacity " << capacities[i] << endl;
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            cout << "facility " << i << " Customer " << j << "   " << servingCosts[i][j] << endl;
        }
    }

    for (int j = 1; j <= numCustomer; j++) {
        cout << "Customer " << j << "   " << demands[j] << endl;
    }
}

/////////////////////////////////////////////////////
//Build original model
void initOriginal() {
    for (int i = 1; i <= numFacility; i++) {
        y[i] = XPRBnewvar(originalSolver, XPRB_BV, XPRBnewname("y_%d", i), 0, 1);
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            x[i][j] = XPRBnewvar(originalSolver, XPRB_PL, XPRBnewname("x_%d_%d", i, j), 0, 1);
        }

    }

    XPRBctr obj = XPRBnewctr(originalSolver, "Obj", XPRB_N);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(obj, y[i], openingCosts[i]);
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            XPRBaddterm(obj, x[i][j], servingCosts[i][j] * demands[j]);
//            XPRBaddterm(obj, x[i][j], servingCosts[i][j]);
        }

    }

//    XPRBprintctr(obj);
    XPRBsetobj(originalSolver, obj);

    XPRBctr totalCapacity = XPRBnewctr(originalSolver, "Sum Capacity >= Sum Demand", XPRB_G);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(totalCapacity, y[i], capacities[i]);
    }


    int totalDemand = 0;
    for (int j = 1; j <= numCustomer; j++) {
        totalDemand += demands[j];
    }
    XPRBaddterm(totalCapacity, NULL, totalDemand);

//    XPRBprintctr(totalCapacity);

    for (int j = 1; j <= numCustomer; j++) {
        XPRBctr demandConserv = XPRBnewctr(originalSolver, XPRBnewname("demand_%d", j), XPRB_E);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(demandConserv, x[i][j], 1.0);

        }
        XPRBaddterm(demandConserv, nullptr, 1);
//        XPRBprintctr(demandConserv);
    }


    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            XPRBctr ctr = XPRBnewctr(originalSolver, XPRBnewname("Bouding %d %d ", i, j), XPRB_L);
            XPRBaddterm(ctr, x[i][j], 1);
            XPRBaddterm(ctr, y[i], -1);
            XPRBaddterm(ctr, nullptr, 0);

//            XPRBprintctr(ctr);
//            boundingCtrs[i][j] = ctr;
        }
    }

    for (int i = 1; i <= numFacility; i++) {
        XPRBctr ctr = XPRBnewctr(originalSolver, XPRBnewname("Capacity %d", i), XPRB_L);
        for (int j = 1; j <= numCustomer; j++) {
            XPRBaddterm(ctr, x[i][j], demands[j]);
//            XPRBaddterm(ctr, x[i][j], 1);
//            cout << -capacities[i] << endl;

        }
        XPRBaddterm(ctr, y[i], -capacities[i]);
        XPRBaddterm(ctr, nullptr, 0);

//        XPRBprintctr(ctr);
//        capacityCtrs[i] = ctr;
    }

//    XPRBprintctr(totalCapacity);

    XPRBsetsense(originalSolver, XPRB_MINIM);

}

void printSubVar() {
    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            if (abs(XPRBgetsol(x[i][j]) - 0) > INT_GAP) {
                cout << "facility " << i << " -> Customer " << j << "   " << XPRBgetsol(x[i][j]) << "    "
                     << XPRBgetsol(x[i][j]) * demands[j] << endl;
            }
        }
    }
}

void verifySolution() {
    for (int i = 1; i <= numFacility; i++) {
        double total = 0;
        for (int j = 1; j <= numCustomer; j++) {
            total += XPRBgetsol(x[i][j]) * demands[j];
        }

        if (total > capacities[i]) {
            cout << "Capacity violation of facility " << i << endl;
        }
    }

    for (int j = 1; j <= numCustomer; j++) {
        double total = 0;
        for (int i = 1; i <= numFacility; i++) {
            total += XPRBgetsol(x[i][j]) * demands[j];
        }

        if (total != demands[j]) {
            cout << "Demand violation of customer " << j << endl;
        }
    }
}

void solveOriginalModel() {
    initOriginal();
    XPRBmipoptimise(originalSolver, "");
    cout << "Optimal -> " << XPRBgetobjval(originalSolver) << endl;

//
    verifySolution();

    XPRBdelprob(originalSolver);
    XPRBfinish;

}

///////////////////////////////////////
//build sub model
void initSubModel() {
    XPRBsetmsglevel(subSolver, 1);
    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            x[i][j] = XPRBnewvar(subSolver, XPRB_PL, XPRBnewname("x_%d_%d", i, j), 0, 1);
        }

    }

    XPRBctr obj = XPRBnewctr(subSolver, XPRBnewname("SubObj"), XPRB_N);
    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            double coeff = servingCosts[i][j] * demands[j];
            XPRBaddterm(obj, x[i][j], servingCosts[i][j] * demands[j]);
        }

    }
    XPRBsetobj(subSolver, obj);
//    XPRBprintctr(obj);


    for (int j = 1; j <= numCustomer; j++) {
        XPRBctr fracConserv = XPRBnewctr(subSolver, XPRBnewname("Sum_i x_i_d% == 1", j), XPRB_E);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(fracConserv, x[i][j], 1.0);
        }
        XPRBaddterm(fracConserv, nullptr, 1.0);

        fracCnsrvtnCtrs[j] = fracConserv;
    }


    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            XPRBctr ctr = XPRBnewctr(subSolver, XPRBnewname("x_%d_d% <= y_d%", i, j, i), XPRB_L);
            XPRBaddterm(ctr, x[i][j], 1);
            XPRBaddterm(ctr, nullptr, 0);

            boundingCtrs[i][j] = ctr;
        }
    }

    for (int i = 1; i <= numFacility; i++) {
        XPRBctr ctr = XPRBnewctr(subSolver, XPRBnewname("Sum_j x_d%_j <= Capacity d%", i, i), XPRB_L);
        for (int j = 1; j <= numCustomer; j++) {
            XPRBaddterm(ctr, x[i][j], demands[j]);
        }
        XPRBaddterm(ctr, nullptr, 0);
        capacityCtrs[i] = ctr;
    }
    XPRBsetsense(subSolver, XPRB_MINIM);

}

void initKnapSack() {
//    map<int, XPRBvar> z;
    for (int j = 1; j <= numCustomer; j++) {
        kpVars[j] = XPRBnewvar(kp, XPRB_PL, XPRBnewname("z_%d", j), 0, 1);
    }
    kpObj = XPRBnewctr(kp, "Obj", XPRB_N);

    for (int j = 1; j <= numCustomer; j++) {
        XPRBaddterm(kpObj, kpVars[j], 1);
    }

    kpCtr = XPRBnewctr(kp, XPRBnewname("capacity"), XPRB_L);
    for (int j = 1; j <= numCustomer; j++) {
        XPRBaddterm(kpCtr, kpVars[j], demands[j]);
    }
    XPRBaddterm(kpCtr, nullptr, MAX);

    XPRBsetobj(kp, kpObj);
}

void solveSubModel() {

//    cout << "Solving sub model.." << endl;
    for (int i = 1; i <= numFacility; i++) {
//        cout << XPRBgetsol(y[i]) << endl;
        for (int j = 1; j <= numCustomer; j++) {
            XPRBsetterm(boundingCtrs[i][j], nullptr, XPRBgetsol(y[i]));
//            XPRBprintctr(boundingCtrs[i][j]);
        }


//        XPRBsetrange(subBoundingCtrs[it->first][i], -MAX, XPRBgetsol(y[i]));
    }

    for (int i = 1; i <= numFacility; i++) {
        XPRBsetterm(capacityCtrs[i], nullptr, XPRBgetsol(y[i]) * capacities[i]);
    }

    XPRBlpoptimise(subSolver, "");

    if (XPRBgetlpstat(subSolver) != XPRB_LP_OPTIMAL) {
        cout << "Sub model -> " << XPRBgetlpstat(subSolver) << endl;
        return;
    }

    for (int j = 1; j <= numCustomer; j++) {
        fracCnsrvtnDuals[j] = XPRBgetdual(fracCnsrvtnCtrs[j]);
    }
}

map<int, double> solveKnapSacks() {

    map<int, double> kpObjVals;
    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            double coeff = 0;
            coeff += demands[j] * servingCosts[i][j];
            coeff -= fracCnsrvtnDuals[j];
            XPRBsetterm(kpObj, kpVars[j], coeff);
        }


        XPRBsetterm(kpCtr, nullptr, capacities[i]);

        XPRBlpoptimise(kp, "");

        if (XPRBgetlpstat(kp) != XPRB_LP_OPTIMAL) {
            cout << "Knapsack status -> " << XPRBgetlpstat(kp) << endl;
            return kpObjVals;
        }

        kpObjVals[i] = XPRBgetobjval(kp);
    }


    return kpObjVals;
}


void addBendersCut() {
    XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("Benders Cut"), XPRB_G);


    double totalDual = 0;
    for (int j = 1; j <= numCustomer; j++) {
        totalDual += fracCnsrvtnDuals[j];
    }

    map<int, double> kpObjVals = solveKnapSacks();

    XPRBaddterm(cut, masterOmega, 1);
    for (int i = 1; i <= numFacility; i++) {
        double coeff = 0;
        coeff += openingCosts[i];
        coeff += kpObjVals[i];
        XPRBaddterm(cut, y[i], -coeff);
    }

    XPRBaddterm(cut, nullptr, totalDual);
}

double computeTotalCost() {
    double currentUb = 0;

    for (int i = 1; i <= numFacility; i++) {
        currentUb += XPRBgetsol(y[i]) * openingCosts[i];
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
//            cout << XPRBgetsol(x[i][j]) << endl;
            currentUb += XPRBgetsol(x[i][j]) * servingCosts[i][j] * demands[j];
        }
    }
    return currentUb;
}

bool isMasterSolutionInteger() {
    for (int i = 1; i <= numFacility; i++) {
        if (abs(XPRBgetsol(y[i]) - round(XPRBgetsol(y[i]))) <= INT_GAP)
            continue;
        return false;
    }
    return true;
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
//        printSubVar();
        nodeUB = computeTotalCost();
        cout << nodeLB << " , " << nodeUB << endl;

        bool useHeuristicSolutionToAddCut = false;
        Solution solution;
//        printMasterVar();
        if (isMasterSolutionInteger()) {
//            if (useOptimalityCut)
//                addOptimalityCut();
            if (nodeUB < ub) {
                ub = nodeUB;

            }

            if (nodeUB < best.totalCost) {
                best.totalCost = nodeUB;
                for (int i = 1; i <= numFacility; i++) {
                    if (abs(XPRBgetsol(y[i]) - 1) <= INT_GAP) {
                        best.selectedLocations[i] = 1;
                    } else {
                        best.selectedLocations[i] = 0;
                    }
                }
            }
        } else {
//            if (useRoundingHeuristic) {
//                solution = roundingHeuristic();
//                if (solution.totalCost < ub) {
//                    ub = solution.totalCost;
//
////                    useHeuristicSolutionToAddCut = true;
//
////                    best.totalCost = ub;
////                    for (int i = 1; i <= numFacility; i++) {
////                        if (solution.selectedLocations[i] == 1) {
////                            best.selectedLocations[i] = 1;
////                        } else {
////                            best.selectedLocations[i] = 0;
////                        }
////                    }
//                }
//
//                if (solution.totalCost < best.totalCost) {
//                    best = solution;
////                    best.totalCost = solution.totalCost;
////                    for (int i = 1; i <= numFacility; i++) {
////                        if (abs(XPRBgetsol(y[i]) - 1) <= INT_GAP) {
////                            best.selectedLocations[i] = 1;
////                        } else {
////                            best.selectedLocations[i] = 0;
////                        }
////                    }
//                }
//            }
        }
        if (abs(nodeUB - nodeLB) <= BD_TERMINATE_GAP)
            break;
        addBendersCut();
//
//        if (useHeuristicSolutionToAddCut) {
////            solveSubModel(solution);
////            addBendersCutForEachSubProblemToMaster(solution);
//        } else {
//            addBendersCutForEachSubProblemToMaster();
//        }
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

void buildBranchingConstraintSet(BranchingConstraintSet &targetSet, int branchingLocationId, int ctrType) {
    BranchingConstraintSet branchingSet = {nullptr, {}};
    branchingSet.basis = XPRBsavebasis(masterSolver);
    for (auto &ctr : targetSet.branchingCons) {
        branchingSet.branchingCons.push_back(ctr);
    }

    int bound = (int) XPRBgetsol(y[branchingLocationId]);
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
        if (abs(XPRBgetsol(y[i]) - round(XPRBgetsol(y[i]))) <= INT_GAP) {
            continue;
        }
        double fractional = XPRBgetsol(y[i]) - (int) XPRBgetsol(y[i]);
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


milliseconds start;

void printMasterVar() {
    for (int i = 1; i <= numFacility; i++) {
        cout << "facility " << i << " -> " << XPRBgetsol(y[i]) << endl;
    }
}

void branchAndCut(bool useOptimalityCutInBendersDecomposition, bool useDualSimplexInBendersDecomposition,
                  bool useRoundingHeuristicInBendersDecomposition,
                  bool useDualSimplexInBranch, bool useInOutStrategyInBendersDecomposition) {
    XPRBsetmsglevel(masterSolver, 1);
    XPRBsetmsglevel(subSolver, 1);
    XPRBsetmsglevel(kp, 1);


    initMaster();
    initSubModel();
    initKnapSack();

    XPRBlpoptimise(masterSolver, "");
//    printMasterVar();
    if (!useInOutStrategyInBendersDecomposition) {
        bendersDecomposition(useOptimalityCutInBendersDecomposition, useDualSimplexInBendersDecomposition,
                             useRoundingHeuristicInBendersDecomposition);
    } else {
//        bendersDecompositionWithInOut(useDualSimplexInBendersDecomposition, useRoundingHeuristicInBendersDecomposition);
    }


    if (isMasterSolutionInteger()) {
        cout << fixed;
        cout << "UB = " << ub << " at the root node " << endl;
//        print();
        return;
    }

    BranchingConstraintSet target = {nullptr, {}};
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
                XPRBaddterm(ctr, y[branching.branchingLocationId], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);

            } else if (branching.ctrType == 1) {
                XPRBctr ctr = XPRBnewctr(masterSolver, branching.name.c_str(), XPRB_G);
                XPRBaddterm(ctr, y[branching.branchingLocationId], 1);
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
}

void readFromFile(const string &fileName) {
    ifstream in(fileName);
    if (!in) {
        cout << "Cannot open input file.\n";
        return;
    }

    string str;
    int lineId = 1;
//    int demandContent =
    while (getline(in, str)) {
        istringstream iss(str);
        if (lineId == 1) {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));

            numFacility = stoi(tokens[0]);
            numCustomer = stoi(tokens[1]);
        } else {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));
            if (tokens.size() == 2) {
                openingCosts[lineId - 1] = stod(tokens[1]);
                capacities[lineId - 1] = stod(tokens[0]);
            } else if (tokens.size() == numCustomer) {
                if (lineId == numFacility + 1 + 1) {
                    for (int id = 0; id < tokens.size(); id++) {
                        demands[id + 1] = stod(tokens[id]);
                    }
                } else {
                    for (int id = 0; id < tokens.size(); id++) {
//                        double t = stod(tokens[id]);
                        servingCosts[lineId - numFacility - 2][id + 1] = stod(tokens[id]);
                    }
                }

            }
        }
        lineId++;
    }
    in.close();
}

void destroy() {
    XPRBdelprob(masterSolver);
    XPRBdelprob(subSolver);
    XPRBdelprob(kp);
    XPRBfinish();
}


int main() {
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/simpleExample2.txt";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/50/50.1";
//        string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/70/70.1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/150/150.5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/200/200.10";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/a/gs250a-5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/b/gs250b-3";
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/cfl/Beasley/capa1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/750/a/gs750a-1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/kmedian/1000-10";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/kmedian/2500-10";
    readFromFile(fileName);
//    printProblem();
    start = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );
    solveOriginalModel();

//    bendersDecomposition();

//    branchAndCut(false,true);
    bool useOptimalityCutInBendersDecomposition = true;
    bool useDualSimplexInBendersDecomposition = false;
    bool useRoundingHeuristicInBendersDecomposition = true;
    bool useDualSimplexInBranch = false;
    bool useInOutStrategyInBendersDecomposition = false;

    bool useLocalBranching = false;

    if (!useLocalBranching) {

        branchAndCut(useOptimalityCutInBendersDecomposition,
                     useDualSimplexInBendersDecomposition,
                     useRoundingHeuristicInBendersDecomposition,
                     useDualSimplexInBranch,
                     useInOutStrategyInBendersDecomposition);

//        branchAndCut2(useOptimalityCutInBendersDecomposition,
//                      useDualSimplexInBendersDecomposition,
//                      useRoundingHeuristicInBendersDecomposition,
//                      useDualSimplexInBranch,
//                      useInOutStrategyInBendersDecomposition);
    } else {
//        branchAndCutWithLocalBranching(useOptimalityCutInBendersDecomposition,
//                                       useDualSimplexInBendersDecomposition,
//                                       useRoundingHeuristicInBendersDecomposition,
//                                       useDualSimplexInBranch,
//                                       useInOutStrategyInBendersDecomposition);
    }


//    branchAndCutWithLocalBranching(true, true, true, false, false);
//    branchAndCutWithInOut(false);
//    branchAndCutWithDualSimplex(false);
    milliseconds end = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    cout << "Time " << (end.count() - start.count()) << endl;
    destroy();
    return 0;
}


