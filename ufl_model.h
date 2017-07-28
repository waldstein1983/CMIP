//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_MODEL_H
#define CMIP_MODEL_H


#include <map>
#include <xprb.h>
#include <cfloat>

#define BD_GAP 1
#define INT_GAP 0.00001
#define MAX DBL_MAX

using namespace std;

int numFacility;
int numCustomer;

double ub = MAX;

map<int, double> openingCosts;
map<int, map<int, double>> servingCosts;


void readFromFile(const string &fileName);


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


void initMaster();

void initSubModel();

void destroy();


double computeTotalCost();

bool solveSubModel();

bool solveSubModel(Solution &solution);

bool solveSubModelWithSeprator(map<int, double> &seprator);

void print();

void print(Solution solution);


#endif //CMIP_MODEL_H
