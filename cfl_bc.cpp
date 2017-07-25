//
// Created by baohuaw on 7/25/17.
//
//#include <iostream>
//#include <chrono>
//#include <xprb.h>
//#include <map>
//#include <vector>
//#include <fstream>
//#include <sstream>
//#include <iterator>
//#include <cmath>
//#include <algorithm>
//#include <cfloat>
//
//#define MAX DBL_MAX
//
//using namespace std;
//using namespace std::chrono;
//
//int nodeNum = 0;
//
//#define INT_GAP 0.00001
//
//double LAMDA = 0.2;
//double DELTA = 2 * INT_GAP;
//
//struct BranchingConstraint {
//    string name;
//    int branchingLocationId;
//    int ctrType;
//    int bound;
//};
//
//struct BranchingConstraintSet {
//    int nodeId;
//    XPRBbasis basis;
//    vector<BranchingConstraint> branchingCons;
//};
//
//milliseconds start;
//
//vector<BranchingConstraintSet> branchingNodes;
//
//int numFacility;
//int numCustomer;
//map<int, XPRBprob> subSolvers;
//XPRBprob masterSolver = XPRBnewprob("master");
//map<int, XPRBvar> masterLocations;
//map<int, XPRBvar> masterAlphas;
//map<int, map<int, XPRBvar>> subCovers;
//
//map<int, map<int, XPRBctr>> subBoundingCtrs;
//
//double ub = MAX;
//
//map<int, double> openingCosts;
//map<int, map<int, double>> servingCosts;
//
//int globalBendersCutId = 1;
//map<int, vector<int>> customerCriticals;
//
//map<int, double> yy;
//
////facility -> customer -> dual
//map<int, map<int, double>> boundingVarSubDuals;
//
//struct Solution {
//    map<int, int> selectedLocations;
//    double totalCost;
//};
//
//Solution best = {{},0};
