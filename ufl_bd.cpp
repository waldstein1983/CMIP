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

using namespace std;
using namespace std::chrono;

int numFacility;
int numCustomer;
map<int, XPRBprob> subSolvers;
XPRBprob masterSolver = XPRBnewprob("master");
map<int, XPRBvar> masterLocations;
map<int, XPRBvar> masterAlphas;
map<int, map<int, XPRBvar>> subCovers;

map<int, map<int, XPRBctr>> subBoundingCtrs;

const double MAX = 1000000;
double lb = -MAX;
double ub = MAX;

map<int, double> openingCosts;
map<int, map<int, double>> servingCosts;

int globalBendersCutId = 1;
map<int, vector<int>> customerCriticals;
const double INT_GAP = 0.00001;


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
//    XPRBsetmsglevel(masterSolver, 1);
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
    lb = XPRBgetobjval(masterSolver);
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
        if (abs(XPRBgetsol(masterLocations[i]) - 1) <= INT_GAP) {
            cout << "Facility " << i << " -> " << XPRBgetsol(masterLocations[i]) << endl;
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
    ub = currentUb;
}


bool addBendersCutForEachSubProblemToMaster() {
    bool newCut = false;
    for (int j = 1; j <= numCustomer; j++) {
        map<int, double> targetServingCosts;
        for (int i = 1; i <= numFacility; i++) {
            targetServingCosts[i] = servingCosts[i][j];
        }

        vector<pair<int, double>> pairs;
        for (auto &targetServingCost : targetServingCosts)
            pairs.push_back(targetServingCost);

        sort(pairs.begin(), pairs.end(), [=](std::pair<int, double> &a, std::pair<int, double> &b) {
                 return a.second < b.second;
             }
        );

        int temp = 0;

        int critical = -1;
        for (auto &pair : pairs) {
            temp += XPRBgetsol(masterLocations[pair.first]);
            if (temp >= 1 && temp - XPRBgetsol(masterLocations[pair.first]) < 1) {
                critical = pair.first;
                break;
            }
        }

        if (customerCriticals.count(j)) {
            if (std::find(customerCriticals[j].begin(), customerCriticals[j].end(), critical) !=
                customerCriticals[j].end()) {
                continue;
            }
        }
        customerCriticals[j].push_back(critical);
        XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("benders cut %d", globalBendersCutId), XPRB_G);
        globalBendersCutId++;
        for (int i = 0; i < pairs.size(); i++) {
            if (pairs[i].first == critical) {
                break;
            } else {
                if (servingCosts[critical][j] - pairs[i].second != 0) {
                    XPRBaddterm(cut, masterLocations[pairs[i].first], servingCosts[critical][j] - pairs[i].second);
                }
            }
        }

        XPRBaddterm(cut, masterAlphas[j], 1.0);
        XPRBaddterm(cut, NULL, servingCosts[critical][j]);
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

bool addBendersCutForEachSubProblemToMaster2() {
    bool newCut = false;
    for (int j = 1; j <= numCustomer; j++) {
        double sumDual  = 0;
        for(int i = 1;i <= numFacility;i++){
//            cout << boundingVarSubDuals[i][j] << " * " << XPRBgetsol(masterLocations[i]) << endl;
            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
        }

        XPRBctr cut = XPRBnewctr(masterSolver, XPRBnewname("benders cut %d", globalBendersCutId), XPRB_G);

        for(int i = 1;i <= numFacility;i++){
            XPRBaddterm(cut, masterLocations[i], -boundingVarSubDuals[i][j]);
//            cout <<  XPRBgetsol(masterLocations[i]) << endl;
//            sumDual += boundingVarSubDuals[i][j] * XPRBgetsol(masterLocations[i]);
        }
        XPRBaddterm(cut, masterAlphas[j],1);
        XPRBaddterm(cut, nullptr, XPRBgetobjval(subSolvers[j]) - sumDual);
//        XPRBprintctr(cut);
        newCut = true;
    }
    return newCut;
}

void bendersDecomposition() {
    XPRBsetmsglevel(masterSolver, 1);
    initMaster();
    initSubModel();
    solveMaster();
    solveSubModel();
    updateUB();
    bool newCutAdded = addBendersCutForEachSubProblemToMaster2();
    if (!newCutAdded)
        return;
    while (abs(ub - lb) >= 1) {
        cout << "$$$$$$$$$$$$$$   " << lb << " , " << ub << endl;
        solveMaster();
        solveSubModel();
        updateUB();
        newCutAdded = addBendersCutForEachSubProblemToMaster2();
        if (!newCutAdded)
            break;
    }
    cout << "UB = " << ub << endl;
    print();


}


int main() {
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/simpleExample2.txt";
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/KoerkelGhosh-sym/250/a/gs250a-2";
    readFromFile(fileName);
    milliseconds start = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );
    bendersDecomposition();

    milliseconds end = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    cout << "Time " << (end.count() - start.count()) << endl;
    destroy();
    return 0;
}