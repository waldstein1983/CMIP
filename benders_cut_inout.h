//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BENDERS_CUT_INOUT_H
#define CMIP_BENDERS_CUT_INOUT_H

#include "ufl_model.h"


double LAMDA = 0.2;
double DELTA = 2 * INT_GAP;

map<int, double> yy;


map<int, double> stabilize();

void addBendersCutWithSeparator(map<int, double> &separator);


Solution bendersDecompositionWithInOut(bool useDualSimplex, bool useRoundingHeuristic);

#endif //CMIP_BENDERS_CUT_INOUT_H
