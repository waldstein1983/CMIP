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

using namespace std;
using namespace std::chrono;

int nodeNum = 0;

struct BranchingConstraint {
    string name;
    int branchingLocationId;
    int ctrType;
    int bound;
};

struct BranchingConstraintSet {
    int nodeId;
    vector<BranchingConstraint> branchingCons;
};

vector<BranchingConstraintSet> branchingNodes;

int numFacility;
int numCustomer;
map<int, XPRBprob> subSolvers;
XPRBprob masterSolver = XPRBnewprob("master");
map<int, XPRBvar> masterLocations;
map<int, XPRBvar> masterAlphas;
map<int, map<int, XPRBvar>> subCovers;

map<int, map<int, XPRBctr>> subBoundingCtrs;

//const double MAX = 10000000000;
//double lb = -MAX;
double ub = MAX;

map<int, double> openingCosts;
map<int, map<int, double>> servingCosts;

int globalBendersCutId = 1;
map<int, vector<int>> customerCriticals;
const double INT_GAP = 0.00001;

//facility -> customer -> dual
map<int, map<int, double>> boundingVarSubDuals;

void print();

void readFromFile(string fileName) {
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
        XPRBaddterm(fulfill, NULL, 1.0);

        for (int i = 1; i <= numFacility; i++) {
            XPRBctr ctr = XPRBnewctr(customer, XPRBnewname("Bounding with facility %d", i), XPRB_L);
            XPRBaddterm(ctr, subCovers[i][j], 1);
            XPRBaddterm(ctr, NULL, 0);
            subBoundingCtrs[j][i] = ctr;
        }

        XPRBsetsense(customer, XPRB_MINIM);
        subSolvers[j] = customer;
    }
}

void solveMaster() {
    XPRBmipoptimise(masterSolver, "");
//    lb = XPRBgetobjval(masterSolver);
//    print();
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

void print() {
    for (int i = 1; i <= numFacility; i++) {
        cout << "Facility " << i << " -> " << XPRBgetsol(masterLocations[i]) << endl;
        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {

        }
    }

    for (int j = 1; j <= numCustomer; j++) {
        cout << "alpha " << j << XPRBgetsol(masterAlphas[j]) << endl;
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            if (abs(XPRBgetsol(subCovers[i][j]) - 1) <= INT_GAP) {
//                cout << "Customer " << j << " -> " << "Facility " << i << endl;
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

void updateUB() {
    double currentUb = 0;

    for (int i = 1; i <= numFacility; i++) {
        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
            currentUb += XPRBgetsol(masterLocations[i]) * openingCosts[i];
        }
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            if (abs(XPRBgetsol(subCovers[i][j]) - 1) <= INT_GAP) {
                currentUb += XPRBgetsol(subCovers[i][j]) * servingCosts[i][j];
            }
        }
    }

//    print();
    if (currentUb < ub)
        ub = currentUb;
}

void addOptimalityCut() {
    XPRBctr cut = XPRBnewctr(masterSolver, "Optimality cut", XPRB_G);
    double sumServingCost = 0;
    for (int j = 1; j <= numCustomer; j++) {
        sumServingCost += XPRBgetobjval(subSolvers[j]);
        XPRBaddterm(cut, masterAlphas[j], 1);

    }
    for (int i = 1; i < numFacility; i++) {
        if(abs(XPRBgetsol(masterLocations[i]) - 0) <= INT_GAP){
            XPRBaddterm(cut, masterLocations[i], sumServingCost);
        }
    }

    XPRBaddterm(cut, nullptr, sumServingCost);
//    XPRBprintctr(cut);
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
//        cout << XPRBgetsol(masterLocations[i]) << endl;
//        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP)
//            continue;
        return false;
    }
    return true;
}

void bendersDecomposition() {

    double nodeLB = 0;
    double nodeUB = MAX;
    while (abs(nodeUB - nodeLB) > 1) {
        if (XPRBgetobjval(masterSolver) > ub) {
            break;
        }
        nodeLB = XPRBgetobjval(masterSolver);
        solveSubModel();
        nodeUB = computeTotalCost();
        if (isMasterSolutionInteger() ) {
            addOptimalityCut();
            if(nodeUB < ub)
                ub = nodeUB;
        }
        if (abs(nodeUB - nodeLB) <= 1)
            break;
        addBendersCutForEachSubProblemToMaster();
        XPRBlpoptimise(masterSolver, "");
    }
}

void buildBranchingConstraintSet(BranchingConstraintSet targetSet, int branchingLocationId, int ctrType) {
    BranchingConstraintSet branchingSet = {0, {}};
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

void branching(BranchingConstraintSet set) {
    int targetBranchingLocationId = 0;


//    for (int i = 1; i <= numFacility; i++) {
//        if (abs(XPRBgetsol(masterLocations[i]) - round(XPRBgetsol(masterLocations[i]))) <= INT_GAP)
//            continue;
//        targetBranchingLocationId = i;
//        break;
//    }

    double gapToHalf = DBL_MAX;
    for(int i = 1;i <= numFacility;i++){
        if (abs(XPRBgetsol(masterLocations[i]) - round(XPRBgetsol(masterLocations[i]))) <= INT_GAP){
            continue;
        }
        double fractional = XPRBgetsol(masterLocations[i]) - (int)XPRBgetsol(masterLocations[i]);
        if (abs(fractional - 0.5) < gapToHalf){
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

void branchAndCut() {
    XPRBsetmsglevel(masterSolver, 1);
    initMaster();
    initSubModel();

    XPRBlpoptimise(masterSolver, "");
    bendersDecomposition();

    if (isMasterSolutionInteger()) {
        cout << "UB = " << ub << " at the root node " << endl;
        return;
    }

//    addOptimalityCut();
//    XPRBlpoptimise(masterSolver, "");
//
//    if (isMasterSolutionInteger()) {
//        double cost = computeTotalCost();
//        if(cost < ub)
//            ub = cost;
//        cout << "UB = " << ub << " at the root node " << endl;
//        return;
//    }

    BranchingConstraintSet target = {nodeNum, {}};
    branching(target);

    int step = 1;
    while (!branchingNodes.empty()) {
        if (step % 1 == 0) {
            cout << ub << "   node size " << branchingNodes.size() << " step = " << step << endl;
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

        XPRBlpoptimise(masterSolver, "");

        if (XPRBgetlpstat(masterSolver) == XPRB_LP_OPTIMAL) {
//            print();
//            cout << XPRBgetobjval(masterSolver) << endl;
            if (XPRBgetobjval(masterSolver) < ub && abs(XPRBgetobjval(masterSolver) - ub) >= INT_GAP) {
                bendersDecomposition();
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


int main() {
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/simpleExample2.txt";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/50/50.1";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/150/150.5";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/GalvaoRaggi/200/200.10";
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/a/gs250a-1";
    readFromFile(fileName);
    milliseconds start = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );
//    bendersDecomposition();

    branchAndCut();
    milliseconds end = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    cout << "Time " << (end.count() - start.count()) << endl;
    destroy();
    return 0;
}