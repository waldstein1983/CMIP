/********************************************************
  Xpress-BCL C++ Example Problems
  ===============================

  file xbelsc.cxx 
  ```````````````
  Economic lot sizing, ELS, problem, solved by adding
  (l,S)-inequalities) in a branch-and-cut heuristic 
  (using the cut manager).
  
  ELS considers production planning over a horizon
  of T periods. In period t, t=1,...,T, there is a
  given demand DEMAND[t] that must be satisfied by
  production prod[t] in period t and by inventory
  carried over from previous periods. There is a
  set-up up cost SETUPCOST[t] associated with
  production in period t. The unit production cost
  in period t is PRODCOST[t]. There is no inventory
  or stock-holding cost.

  (c) 2008 Fair Isaac Corporation
      author: S.Heipcke, 2005, rev. Mar. 2011
********************************************************/

#include <iostream>
#include "xprb_cpp.h"
#include "xprs.h"

using namespace std;
using namespace ::dashoptimization;

#define T 6                             /* Number of time periods */

/****DATA****/
int DEMAND[] = {1, 3, 5, 3, 4, 2};  /* Demand per period */
int SETUPCOST[] = {17, 16, 11, 6, 9, 6};  /* Setup cost per period */
int PRODCOST[] = {5, 3, 2, 1, 3, 1};  /* Production cost per period */
int D[T][T];                            /* Total demand in periods t1 - t2 */

XPRBvar prod[T];                        /* Production in period t */
XPRBvar setup[T];                       /* Setup in period t */

struct myobj {
    XPRBprob *prob;
    double tol;
};

XPRBprob p("Els");                      /* Initialize a new problem in BCL */

/***********************************************************************/

void modEls() {
    int s, t, k;
    XPRBexpr cobj, le;

    for (s = 0; s < T; s++)
        for (t = 0; t < T; t++)
            for (k = s; k <= t; k++)
                D[s][t] += DEMAND[k];

/****VARIABLES****/
    for (t = 0; t < T; t++) {
        prod[t] = p.newVar(XPRBnewname("prod%d", t + 1));
        setup[t] = p.newVar(XPRBnewname("setup%d", t + 1), XPRB_BV);
    }

/****OBJECTIVE****/
    for (t = 0; t < T; t++)                       /* Minimize total cost */
        cobj += SETUPCOST[t] * setup[t] + PRODCOST[t] * prod[t];
    p.setObj(cobj);

/****CONSTRAINTS****/
    /* Production in period t must not exceed the total demand for the
       remaining periods; if there is production during t then there
       is a setup in t */
    for (t = 0; t < T; t++)
        p.newCtr("Production", prod[t] <= D[t][T - 1] * setup[t]);

    /* Production in periods 0 to t must satisfy the total demand
       during this period of time */
    for (t = 0; t < T; t++) {
        le = 0;
        for (s = 0; s <= t; s++) le += prod[s];
        p.newCtr("Demand", le >= D[0][t]);
    }

}

/**************************************************************************/
/*  Cut generation loop at the tree node:                                 */
/*    get the solution values                                             */
/*    identify and set up violated constraints                            */
/*    add cuts to the matrix                                              */
/**************************************************************************/
int XPRS_CC cbNode(XPRSprob oprob, void *mobj) {
    struct myobj *mo;
    double objval;                  /* Objective value */
    int t, l;
    int ncut;                       /* Counters for cuts */
    double solprod[T], solsetup[T]; /* Solution values for var.s prod & setup */
    double ds;
    int depth, node;
    XPRBcut cut[T];
    XPRBexpr le;

    mo = (struct myobj *) mobj;
    mo->prob->beginCB(oprob);

    ncut = 0;
    XPRSgetintattrib(oprob, XPRS_NODEDEPTH, &depth);
    XPRSgetintattrib(oprob, XPRS_NODES, &node);

    /* Get the solution values */
    mo->prob->sync(XPRB_XPRS_SOL);
    for (t = 0; t < T; t++) {
        solprod[t] = prod[t].getSol();
        solsetup[t] = setup[t].getSol();
    }

    /* Search for violated constraints: */
    for (l = 0; l < T; l++) {
        for (ds = 0.0, t = 0; t <= l; t++) {
            if (solprod[t] < D[t][l] * solsetup[t] + mo->tol) ds += solprod[t];
            else ds += D[t][l] * solsetup[t];
        }

        /* Add the violated inequality: the minimum of the actual production
           prod[t] and the maximum potential production D[t][l]*setup[t]
           in periods 0 to l must at least equal the total demand in periods
           0 to l.
           sum(t=1:l) min(prod[t], D[t][l]*setup[t]) >= D[0][l]
         */
        if (ds < D[0][l] - mo->tol) {
            le = 0;
            for (t = 0; t <= l; t++) {
                if (solprod[t] < D[t][l] * solsetup[t] + mo->tol)
                    le += prod[t];
                else
                    le += D[t][l] * setup[t];
            }
            cut[ncut] = mo->prob->newCut(le >= D[0][l]);
            ncut++;
        }
    }

/* Add cuts to the problem */
    if (ncut > 0) {
        mo->prob->addCuts(cut, ncut);
        XPRSgetdblattrib(oprob, XPRS_LPOBJVAL, &objval);
        cout << "Cuts added : " << ncut << " (depth " << depth << ", node ";
        cout << node << "), obj. " << objval << endl;
    }
    mo->prob->endCB();

    return 0;
}

/***********************************************************************/
void treeCutGen() {
    XPRSprob oprob;
    struct myobj mo;
    double feastol;
    int starttime, t;

    starttime = XPRB::getTime();

    oprob = p.getXPRSprob();                         /* Get Optimizer problem */

    XPRSsetintcontrol(oprob, XPRS_LPLOG, 0);
    XPRSsetintcontrol(oprob, XPRS_MIPLOG, 3);

    XPRSsetintcontrol(oprob, XPRS_CUTSTRATEGY, 0);   /* Disable automatic cuts */
    XPRSsetintcontrol(oprob, XPRS_PRESOLVE, 0);      /* Switch presolve off */
    XPRSsetintcontrol(oprob, XPRS_EXTRAROWS, 5000);  /* Reserve extra rows */

    XPRSgetdblcontrol(oprob, XPRS_FEASTOL, &feastol);  /* Get zero tolerance */
    feastol *= 10;

    mo.prob = &p;
    mo.tol = feastol;
    p.setCutMode(1);
    XPRSsetcbcutmgr(oprob, cbNode, &mo);

    p.mipOptimize("");                               /* Solve the MIP */
    cout << "(" << (XPRB::getTime() - starttime) / 1000.0 << " sec) Global status ";
    cout << p.getMIPStat() << ", best solution: " << p.getObjVal() << endl;
    for (t = 0; t < T; t++) {
        cout << "Period " << t + 1 << ": prod " << prod[t].getSol() << " (demand: ";
        cout << DEMAND[t] << ", cost: " << PRODCOST[t] << "), setup ";
        cout << setup[t].getSol() << " (cost: " << SETUPCOST[t] << ")" << endl;
    }
}

/***********************************************************************/

int main(int argc, char **argv) {
    modEls();                    /* Model the problem */
    treeCutGen();                /* Solve the problem */

    return 0;
} 
