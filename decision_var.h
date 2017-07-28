//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_DECISION_VAR_H
#define CMIP_DECISION_VAR_H

#include <xprb.h>
#include <map>

using namespace std;

map<int, XPRBvar> masterLocations;
map<int, XPRBvar> masterAlphas;
map<int, map<int, XPRBvar>> subCovers;


#endif //CMIP_DECISION_VAR_H
