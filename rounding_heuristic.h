//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_ROUNDING_HEURISTIC_H
#define CMIP_ROUNDING_HEURISTIC_H

#include "model.h"
#include "problem.h"
#include <algorithm>

Solution roundingHeuristic() {
    double sortedLocation[numFacility + 1];
    for (int i = 1; i <= numFacility; i++) {
        sortedLocation[i] = XPRBgetsol(masterLocations[i]);
    }

    sort(sortedLocation + 1, sortedLocation + numFacility + 1);

    vector<Solution> solutions;
    double curThreshold = -1;
    for (int i = 1; i <= numFacility; i++) {
        if (sortedLocation[i] != curThreshold) {
            Solution solution = {{}, 0};
            curThreshold = sortedLocation[i];
            int roundingLocation[numFacility + 1];
            for (int ii = 1; ii <= numFacility; ii++) {
                if (abs(XPRBgetsol(masterLocations[ii]) - 1) <= INT_GAP) {
                    roundingLocation[ii] = 1;
                    solution.totalCost += openingCosts[ii];
                } else {
                    if (XPRBgetsol(masterLocations[ii]) < curThreshold) {
                        roundingLocation[ii] = 0;
                    } else {
                        solution.totalCost += openingCosts[ii];
                        roundingLocation[ii] = 1;
                    }
                }
            }

            bool locationSelected[numFacility + 1];

            for (int ii = 1; ii <= numFacility; ii++) {
                locationSelected[ii] = false;
            }

            for (int j = 1; j <= numCustomer; j++) {
                int neareastFacility = 0;
                double nearestDistance = MAX;
                for (int ii = 1; ii <= numFacility; ii++) {
                    if (roundingLocation[ii] == 1) {
                        if (servingCosts[ii][j] < nearestDistance) {
                            nearestDistance = servingCosts[ii][j];
                            neareastFacility = ii;
                        }
                    }
                }
                solution.totalCost += servingCosts[neareastFacility][j];
                locationSelected[neareastFacility] = true;
            }

            for (int ii = 1; ii <= numFacility; ii++) {
                solution.selectedLocations[ii] = 0;
                if (!locationSelected[ii] && roundingLocation[ii] == 1) {
                    solution.totalCost -= openingCosts[ii];
                } else if (locationSelected[ii] && roundingLocation[ii] == 1) {
                    solution.selectedLocations[ii] = 1;
                }
            }
//            cout << i << "    " << solution.totalCost << endl;
            solutions.push_back(solution);
        }
    }


    sort(solutions.begin(), solutions.begin() + solutions.size(),
         [](Solution const &a, Solution const &b) { return a.totalCost < b.totalCost; });

    return solutions[0];
}

#endif //CMIP_ROUNDING_HEURISTIC_H
