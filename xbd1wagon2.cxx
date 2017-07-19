/********************************************************
Xpress-BCL C++ Example Problems
===============================

file d1wagon2.cpp
````````````````
Load balancing of train wagons
(second version, using heuristic solution as 
 start solution for MIP)
   
(c) 2014 Fair Isaac Corporation
author: L.Bertacco, September 2014
********************************************************/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <vector>
#include <algorithm>
#include "xprb_cpp.h"
#include "xprs.h"

using namespace std;
using namespace ::dashoptimization;

#define NBOXES (sizeof(WEIGHT)/sizeof(*WEIGHT)) /* Number of boxes                       */
#define NWAGONS 3                               /* Number of wagons                      */

/* Box weights                           */
int WEIGHT[] = {34, 6, 8, 17, 16, 5, 13, 21, 25, 31, 14, 13, 33, 9, 25, 25};
int WMAX = 100;                                 /* Weight limit of the wagons            */

int HeurSol[NBOXES];                            /* Heuristic solution: for each box      */
/* specifies in which wagon it is loaded */

/****VARIABLES****/
XPRBvar load[NBOXES][NWAGONS];                  /* 1 if box loaded on wagon, 0 otherwise */
XPRBvar maxweight;                              /* Weight of the heaviest wagon load     */

XPRBprob prob;

void XPRS_CC solnotify(XPRSprob my_prob, void *my_object, const char *solname, int status);

void d1w2_model(XPRBprob &prob) {
    /****VARIABLES****/

    /* Create load[box,wagon] (binary) */
    for (int b = 0; b < NBOXES; b++)
        for (int w = 0; w < NWAGONS; w++)
            load[b][w] = prob.newVar(XPRBnewname("load_%d_%d", b + 1, w + 1), XPRB_BV);

    /* Create maxweight (continuous with lb=ceil((sum(b in BOXES) WEIGHT(b))/NBOXES) */
    double sum_weights = 0;
    for (int b = 0; b < NBOXES; b++) sum_weights += WEIGHT[b];
    maxweight = prob.newVar("maxweight", XPRB_PL, ceil(sum_weights / NBOXES), XPRB_INFINITY);

    /****CONSTRAINTS****/

    /* Every box into one wagon: forall(b in BOXES) sum(w in WAGONS) load(b,w) = 1 */
    for (int b = 0; b < NBOXES; b++) {
        XPRBexpr eq;
        for (int w = 0; w < NWAGONS; w++) eq += load[b][w];
        prob.newCtr(eq == 1);
    }

    /* Limit the weight loaded into every wagon: forall(w in WAGONS) sum(b in BOXES) WEIGHT(b)*load(b,w) <= maxweight */
    for (int w = 0; w < NWAGONS; w++) {
        XPRBexpr le;
        for (int b = 0; b < NBOXES; b++) le += WEIGHT[b] * load[b][w];
        prob.newCtr(le <= maxweight);
    }

    /****OBJECTIVE****/

    prob.setObj(maxweight);
    prob.setSense(XPRB_MINIM);
}

void d1w2_solve(XPRBprob &prob) {
    int b, w;

    /* Alternative to lower bound on maxweight: adapt the optimizer cutoff value  */
    /* XPRSsetdblcontrol(XPRBgetXPRSprob(prob), XPRS_MIPADDCUTOFF, -0.99999); */

    /* Comment out the following line to enable the optimizer log */
    XPRSsetintcontrol(prob.getXPRSprob(), XPRS_OUTPUTLOG, 0);

    /* Create a BCL solution from the heuristic solution we have found */
    XPRBsol sol = prob.newSol();
    /* Set the solution values for all discrete variables that are non-zero */
    for (b = 0; b < NBOXES; b++) sol.setVar(load[b][HeurSol[b]], 1);

    /* It is possible, but not necessary, to set values for ALL discrete vars  */
    /* by uncommenting the following line. In this case, the usersolnotify     */
    /* callback would return status equal to 2 (instead of 3), as the solution */
    /* would be feasible without the need of a local search.                   */
    /* for (b=0; b<NBOXES; b++) for (w=0; w<NWAGONS; w++) XPRBsetsolvar(sol, load[b][w], w==HeurSol[b]); */

    prob.addMIPSol(sol, "heurSol");      /* Send the solution to the optimizer */
    /* Request notification of solution status after processing */
    XPRSaddcbusersolnotify(prob.getXPRSprob(), solnotify, NULL, 0);

    /* Parameter settings to make use of loaded solution */
    XPRSsetdblcontrol(prob.getXPRSprob(), XPRS_HEURSEARCHEFFORT, 2);
    XPRSsetintcontrol(prob.getXPRSprob(), XPRS_HEURSEARCHROOTSELECT, 31);
    XPRSsetintcontrol(prob.getXPRSprob(), XPRS_HEURSEARCHTREESELECT, 19);

    prob.mipOptimize();              /* Solve the LP-problem */
    int statmip = prob.getMIPStat(); /* Get the problem status */
    if (statmip == XPRB_MIP_SOLUTION || statmip == XPRB_MIP_OPTIMAL) { /* An integer solution has been found */
        cout << "Optimal solution:\n Max weight: " << prob.getObjVal() << endl;
        for (w = 0; w < NWAGONS; w++) {
            int tot_weight = 0;
            cout << " " << (w + 1) << ":";
            for (b = 0; b < NBOXES; b++)
                if (load[b][w].getSol() > .5) {
                    cout << " " << (b + 1);
                    tot_weight += WEIGHT[b];
                }
            cout << " (total weight: " << tot_weight << ")" << endl;
        }
    }
}

/***********************************************************************/

/* LPT (Longest processing time) heuristic:     */
/* One at a time, place the heaviest unassigned */
/* box onto the wagon with the least load       */
bool weight_greater(int i, int j) { return WEIGHT[i] > WEIGHT[j]; }

double solve_heur() {
    vector<int> ORDERW(
            NBOXES);           /* Box indices sorted in decreasing weight order                                              */
    int CurNum[NWAGONS] = {
            0};          /* For each wagon w, this is the number of boxes currently loaded                             */
    int CurWeight[NWAGONS] = {
            0};       /* For each wagon w, this is the current weight, i.e. the sum of weights of loaded boxes      */
    int Load[NWAGONS][NBOXES] = {
            0};    /* Load[w][i] (for i=0..CurNum[w]-1) contains the box index of the i-th box loaded on wagon w */

    /* Copy the box indices into array ORDERW and sort them in decreasing     */
    /* order of box weights (the sorted indices are returned in array ORDERW) */
    for (int b = 0; b < NBOXES; b++) ORDERW[b] = b;
    sort(ORDERW.begin(), ORDERW.end(), weight_greater);

    /* Distribute the loads to the wagons using the LPT heuristic  */
    for (int b = 0; b < NBOXES; b++) {
        int v = 0;                          /* Find wagon v with the smallest load */
        for (int w = 0; w < NWAGONS; w++) if (CurWeight[w] <= CurWeight[v]) v = w;
        Load[v][CurNum[v]] = ORDERW[b];     /* Add current box to wagon v */
        CurNum[v]++;                        /* Increase the counter of boxes on v */
        CurWeight[v] += WEIGHT[ORDERW[b]];  /* Update current weight of the wagon */
    }

    /* Calculate the solution value */
    double heurobj = 0;                   /* heuristic solution objective value (max wagon weight) */
    for (int w = 0; w < NWAGONS; w++) if (CurWeight[w] > heurobj) heurobj = CurWeight[w];

    /* Solution printing */
    cout << "Heuristic solution:\n Max weight: " << heurobj << endl;
    for (int w = 0; w < NWAGONS; w++) {
        cout << " " << (w + 1) << ":";
        for (int i = 0; i < CurNum[w]; i++) cout << " " << (Load[w][i] + 1);
        cout << " (total weight: " << CurWeight[w] << ")" << endl;
    }

    /* Save the heuristic solution into the HeurSol array */
    for (int w = 0; w < NWAGONS; w++) for (int i = 0; i < CurNum[w]; i++) HeurSol[Load[w][i]] = w;

    return heurobj;
}

/* Callback function reporting loaded solution status */
void XPRS_CC solnotify(XPRSprob my_prob, void *my_object, const char *solname, int status) {
    cout << "Optimizer loaded solution " << (solname ? solname : "(null)") << " with status=" << status << endl;
}

/***********************************************************************/

int main(int argc, char **argv) {
    XPRBprob prob("d1wagon2"); /* Initialize a new problem in BCL */;

    if (solve_heur() <= WMAX) {
        cout << "Heuristic solution fits capacity limits" << endl;
        exit(0);
    }

    d1w2_model(prob);             /* Model the problem */
    d1w2_solve(prob);             /* Solve the problem */

    return 0;
}

