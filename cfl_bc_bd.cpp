//
// Created by baohuaw on 8/1/17.
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

#define BD_GAP 1
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

map<int, int> demands;


map<int, XPRBvar> masterLocations;
XPRBvar masterOmega;
map<int, map<int, XPRBvar>> fcltCstmFracs;

XPRBprob subSolver = XPRBnewprob("sub");
XPRBprob masterSolver = XPRBnewprob("master");
XPRBprob kp = XPRBnewprob("knapSack");
XPRBctr kpObj;
XPRBctr kpCtr;
map<int, XPRBvar > kpVars;
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

vector<BranchingConstraintSet> branchingNodes;

/////////////////////////////////////////////////////
//Build master model
void initMaster() {
    for (int i = 1; i <= numFacility; i++) {
        masterLocations[i] = XPRBnewvar(masterSolver, XPRB_UI, XPRBnewname("y_%d", i), 0, 1);
    }

    masterOmega = XPRBnewvar(masterSolver, XPRB_PL, XPRBnewname("omega"), 0, MAX);


    XPRBctr obj = XPRBnewctr(masterSolver, "Obj", XPRB_N);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(obj, masterLocations[i], openingCosts[i]);
    }
    XPRBaddterm(obj, masterOmega, 1.0);

    XPRBsetobj(masterSolver, obj);

    XPRBctr totalCapacity = XPRBnewctr(masterSolver, "Sum Capacity >= Sum Demand", XPRB_G);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(totalCapacity, masterLocations[i], capacities[i]);
    }

    int totalDemand = 0;
    for (int j = 1; j <= numCustomer; j++) {
        totalDemand += demands[j];
    }
    XPRBaddterm(totalCapacity, NULL, totalDemand);

    XPRBsetsense(masterSolver, XPRB_MINIM);

}

///////////////////////////////////////
//build sub model
void initSubModel() {
    XPRBsetmsglevel(subSolver, 1);
    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            fcltCstmFracs[i][j] = XPRBnewvar(subSolver, XPRB_PL, XPRBnewname("x_%d_%d", i, j), 0, 1);
        }

    }

    XPRBctr obj = XPRBnewctr(subSolver, XPRBnewname("SubObj"), XPRB_N);
    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            XPRBaddterm(obj, fcltCstmFracs[i][j], servingCosts[i][j] * demands[j]);
        }

    }
    XPRBsetobj(subSolver, obj);


    for (int j = 1; j <= numCustomer; j++) {
        XPRBctr fracConserv = XPRBnewctr(subSolver, XPRBnewname("Sum_i x_i_d% == 1", j), XPRB_E);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(fracConserv, fcltCstmFracs[i][j], 1.0);
        }
        XPRBaddterm(fracConserv, nullptr, 1.0);

        fracCnsrvtnCtrs[j] = fracConserv;
    }


    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            XPRBctr ctr = XPRBnewctr(subSolver, XPRBnewname("x_%d_d% <= y_d%", i, j, i), XPRB_L);
            XPRBaddterm(ctr, fcltCstmFracs[i][j], 1);
            XPRBaddterm(ctr, masterLocations[i], -1);
            XPRBaddterm(ctr, nullptr, 0);
        }
    }

    for (int i = 1; i <= numFacility; i++) {
        XPRBctr ctr = XPRBnewctr(subSolver, XPRBnewname("Sum_j x_d%_j <= Capacity d%", i, i), XPRB_L);
        for (int j = 1; j <= numCustomer; j++) {
            XPRBaddterm(ctr, fcltCstmFracs[i][j], demands[j]);
        }
        XPRBaddterm(ctr, masterLocations[i], capacities[i]);
    }
    XPRBsetsense(subSolver, XPRB_MINIM);

}

void initKnapSack(){
//    map<int, XPRBvar> z;
    for (int j = 1; j <= numCustomer; j++) {
        kpVars[j] = XPRBnewvar(kp, XPRB_PL, XPRBnewname("z_%d", j), 0, 1);
    }
    kpObj = XPRBnewctr(kp, "Obj", XPRB_N);

    for (int j = 1; j <= numCustomer; j++) {
        XPRBaddterm(kpObj, kpVars[j], 1);
    }

    kpCtr = XPRBnewctr(subSolver, XPRBnewname("capacity"), XPRB_L);
    for (int j = 1; j <= numCustomer; j++) {
        XPRBaddterm(kpCtr, kpVars[j], demands[j]);
    }
    XPRBaddterm(kpCtr, nullptr, MAX);

    XPRBsetobj(kp, kpObj);
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

        kpObjVals[i] = XPRBgetobjval(kp);
    }


    return kpObjVals;
}


void addBendersCut() {
    XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("Benders Cut d%", globalBendersCutId), XPRB_G);
    double totalDual = 0;
    for (int j = 1; j <= numCustomer; j++) {
        totalDual += fracCnsrvtnDuals[j];
    }

    map<int ,double> kpObjVals = solveKnapSacks();


}

