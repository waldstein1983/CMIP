#include <iostream>


#include <string>
#include <xprb.h>
#include <map>
#include <cfloat>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <xprb_cpp.h>
#include <cmath>
#include <algorithm>

using namespace std;
//using namespace ::dashoptimization;

int numFacility;
int numCustomer;
map<int, XPRBprob> subSolvers;
XPRBprob masterSolver = XPRBnewprob("master");
map<int, XPRBvar> masterLocations;
map<int, XPRBvar> masterAlphas;
map<int, map<int, XPRBvar>> subCovers;

map<int, map<int, XPRBctr>> subBoundingCtrs;

double lb = -DBL_MAX;
double ub = DBL_MAX;
//vector<string> complicatingVarNames;
map<int, double> openingCosts;
map<int, map<int, double>> servingCosts;

int globalBendersCutId = 1;
map<int, vector<int>> customerCriticals;
int INT_GAP = 0.00001;


map<string, map<int, double>> boundingVarSubDuals;

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
        masterAlphas[j] = XPRBnewvar(masterSolver, XPRB_PL, XPRBnewname("alpha_%d", j), 0, DBL_MAX);
    }

    XPRBctr obj = XPRBnewctr(masterSolver, "Obj", XPRB_N);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(obj, masterLocations[i], openingCosts[i]);
    }
    XPRBsetobj(masterSolver, obj);

    XPRBctr facilityExistence = XPRBnewctr(masterSolver, "Facility_existence", XPRB_E);
    for (int i = 1; i <= numFacility; i++) {
        XPRBaddterm(facilityExistence, masterLocations[i], 1.0);
    }
    XPRBsetrange(facilityExistence, 1, 1);
//    XPRBprintctr(facilityExistence);
//    masterSolver.newCtr("Facility_exist", facilityExistence >= 1);
//    masterSolver.setSense(XPRB_MINIM);
    XPRBsetsense(masterSolver, XPRB_MINIM);
    XPRBsetmsglevel(masterSolver, 1);
}

void initSubModel() {
    for (int j = 1; j <= numCustomer; j++) {
        XPRBprob customer = XPRBnewprob("Sub");
        XPRBsetmsglevel(customer, 1);
        for (int i = 1; i <= numFacility; i++) {
            string varName = "x_" + i;
            varName.append("_" + j);
            subCovers[i][j] = XPRBnewvar(customer, XPRB_PL, varName.c_str(), 0, DBL_MAX);
        }

        XPRBctr obj = XPRBnewctr(customer, "Obj of Sub problem " + j, XPRB_N);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(obj, subCovers[i][j], servingCosts[i][j]);
        }
        XPRBsetobj(customer, obj);

        XPRBctr fulfill = XPRBnewctr(customer, "Fullfill of Sub problem " + j, XPRB_E);
        for (int i = 1; i <= numFacility; i++) {
            XPRBaddterm(fulfill, subCovers[i][j], 1.0);
        }
        XPRBsetrange(fulfill, 1, 1);

        for (int i = 1; i <= numFacility; i++) {
            XPRBctr ctr = XPRBnewctr(customer, "Bounding with facility " + i, XPRB_L);
            XPRBaddterm(ctr, subCovers[i][j], 1);
            XPRBsetrange(ctr, -DBL_MAX, 0);
            subBoundingCtrs[j][i] = ctr;
        }

        XPRBsetsense(customer, XPRB_MINIM);
        subSolvers[j] = customer;
    }
}

void solveMaster() {
    XPRBmipoptimise(masterSolver, "");
//    masterSolver.mipOptimize();
    lb = XPRBgetobjval(masterSolver);
//    lb = masterSolver.getObjVal();
//    for (int i = 1; i <= numFacility; i++) {
//        cout << "facility " << i << " : " << XPRBgetsol(masterLocations[i]) << endl;
//    }
}

bool solveSubModel() {
    for (map<int, XPRBprob>::iterator it = subSolvers.begin(); it != subSolvers.end(); ++it) {
        for (int i = 1; i <= numFacility; i++) {
            XPRBsetrange(subBoundingCtrs[it->first][i], -DBL_MAX, XPRBgetsol(masterLocations[i]));
        }

        XPRBlpoptimise(subSolvers[it->first], "");
        for (int i = 1; i <= numFacility; i++) {
            string ctrName = "Bounding with y_" + i;
            boundingVarSubDuals["y_" + i][it->first] = XPRBgetdual(subBoundingCtrs[it->first][i]);
        }
    }
    return true;
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

    ub = currentUb;
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

bool addBendersCutForEachSubProblemToMaster() {

    XPRBctr ctr;
    ctr = NULL;
    while ((ctr = XPRBgetnextctr(masterSolver, ctr)) != NULL) {
        XPRBprintctr(ctr);
    }

    XPRBcut cut[numCustomer];
    bool newCut = false;
    for (int j = 1; j <= numCustomer; j++) {
        map<int, double> targetServingCosts;
        for (int i = 1; i <= numFacility; i++) {
            targetServingCosts[i] = servingCosts[i][j];
        }

        vector<pair<int, double>> pairs;
        for (auto itr = targetServingCosts.begin(); itr != targetServingCosts.end(); ++itr)
            pairs.push_back(*itr);

        sort(pairs.begin(), pairs.end(), [=](std::pair<int, double> &a, std::pair<int, double> &b) {
                 return a.second < b.second;
             }
        );

        int temp = 0;

        int critical = -1;
//        map<string, Double> cutTerms = new LinkedHashMap<>();
        for (int i = 0; i < pairs.size(); i++) {
            temp += XPRBgetsol(masterLocations[pairs[i].first]);
            if (temp >= 1 && temp - XPRBgetsol(masterLocations[pairs[i].first]) < 1) {
                critical = pairs[i].first;
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

//        string cutId("Benders Cut " + j);
//        cutId += " ";
//        cutId += globalBendersCutId;
        XPRBctr cut = XPRBnewctr(masterSolver, "benders cut " + globalBendersCutId, XPRB_G);

//        XPRBcut cc = masterSolver.newCut(0);
//        XPRBcut cut = XPRBnewcut(masterSolver, XPRB_G, globalBendersCutId);

//        XPRBexpr cut;
//        cut += masterLocations[1];
        globalBendersCutId++;
        for (int i = 0; i < pairs.size(); i++) {
            if (pairs[i].first == critical) {
                break;
            } else {
                if (servingCosts[critical][j] - pairs[i].second != 0) {
//                    cut += masterLocations[1];

                    XPRBaddterm(cut, masterLocations[pairs[i].first], servingCosts[critical][j] - pairs[i].second);
//                    XPRBaddcutterm(cut, masterLocations[pairs[i].first], servingCosts[critical][j] - pairs[i].second);
//                    cut.addTerm();
                }
            }
        }

        XPRBaddterm(cut, masterAlphas[j], 1.0);

        XPRBsetrange(cut, servingCosts[critical][j], DBL_MAX);

        XPRBprintctr(cut);
//        XPRBaddcutterm(cut, NULL, servingCosts[critical][j]);
//        XPRBset
//        XPRBsetrange(cut, servingCosts[critical][j], DBL_MAX);
//        string varName = "alpha_" + j;
//        cut.addTerm(masterAlphas[j], 1.0);

//        string cutName = "Benders cut " + j;
//        cutName.append(" " + globalBendersCutId);
//        masterSolver.newCtr(cutName.c_str(), cut >= servingCosts[critical][j]);
        newCut = true;
    }

//    XPRBctr ctr;
//    ctr = NULL;
//    while((ctr = XPRBgetnextctr(masterSolver, ctr)) != NULL){
//        XPRBprintctr(ctr);
//    }
//    XPRBgetc
    return newCut;
}

void bendersDecomposition() {
    initMaster();
    initSubModel();
    solveMaster();
    solveSubModel();
    updateUB();
    bool newCutAdded = addBendersCutForEachSubProblemToMaster();
    if (newCutAdded == false)
        return;
    while (abs(ub - lb) >= 1) {
        cout << lb << " , " << ub << endl;
        solveMaster();
        solveSubModel();
        updateUB();
        newCutAdded = addBendersCutForEachSubProblemToMaster();
        if (newCutAdded == false)
            break;
    }
    cout << "UB = " << ub << endl;
    print();
}


int main() {
//    std::cout << "Hello, World!" << std::endl;
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/ufl/simpleExample2.txt";
    readFromFile(fileName);
    bendersDecomposition();
    return 0;
}