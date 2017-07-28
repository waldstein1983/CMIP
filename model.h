//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_MODEL_H
#define CMIP_MODEL_H

#define BD_GAP 1
#define INT_GAP 0.00001

#include <xprb.h>
#include <map>
#include "problem.h"


using namespace std;

map<int, XPRBvar> masterLocations;
map<int, XPRBvar> masterAlphas;
map<int, map<int, XPRBvar>> subCovers;

map<int, XPRBprob> subSolvers;
XPRBprob masterSolver = XPRBnewprob("master");
map<int, map<int, XPRBctr>> subBoundingCtrs;

//facility -> customer -> dual

map<int, map<int, double>> boundingVarSubDuals;

struct Solution {
    map<int, int> selectedLocations;
    double totalCost;
};



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

void destroy() {
    XPRBdelprob(masterSolver);
    for (auto &subSolver : subSolvers) {
        XPRBdelprob(subSolver.second);
    }
    XPRBfinish();

}


double computeTotalCost() {
    double currentUb = 0;

    for (int i = 1; i <= numFacility; i++) {
        currentUb += XPRBgetsol(masterLocations[i]) * openingCosts[i];
    }

    for (int i = 1; i <= numFacility; i++) {
        for (int j = 1; j <= numCustomer; j++) {
            currentUb += XPRBgetsol(subCovers[i][j]) * servingCosts[i][j];
        }
    }
    return currentUb;
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

#endif //CMIP_MODEL_H
