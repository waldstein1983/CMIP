//
// Created by baohuaw on 7/28/17.
//
//#include <xprb.h>
#include <chrono>
#include "branch_cut.h"
#include "ufl_model.h"
#include "benders_cut.h"
#include "branching.h"
#include "benders_cut_inout.h"

using namespace std::chrono;

milliseconds start;

void branchAndCut(bool useOptimalityCutInBendersDecomposition, bool useDualSimplexInBendersDecomposition,
                  bool useRoundingHeuristicInBendersDecomposition,
                  bool useDualSimplexInBranch, bool useInOutStrategyInBendersDecomposition) {
    XPRBsetmsglevel(masterSolver, 1);
    initMaster();
    initSubModel();

    XPRBlpoptimise(masterSolver, "");
    if (!useInOutStrategyInBendersDecomposition) {
        bendersDecomposition(useOptimalityCutInBendersDecomposition, useDualSimplexInBendersDecomposition,
                             useRoundingHeuristicInBendersDecomposition);
    } else {
        bendersDecompositionWithInOut(useDualSimplexInBendersDecomposition, useRoundingHeuristicInBendersDecomposition);
    }


    if (isMasterSolutionInteger()) {
        cout << "UB = " << ub << " at the root node " << endl;
        print();
        return;
    }

    BranchingConstraintSet target = {nullptr, {}};
    branching(target);

    int step = 1;
    while (!branchingNodes.empty()) {
        milliseconds end = duration_cast<milliseconds>(
                system_clock::now().time_since_epoch()
        );
        if (step % 1 == 0) {
            cout << ub << "   node size " << branchingNodes.size() << " step = " << step << "   Time "
                 << (end.count() - start.count()) << endl;

        }
        if ((end.count() - start.count()) / 1000 > 720000) {
            cout << "Terminate due to time limit of 3600 sec" << endl;
            break;
        }
        target = branchingNodes[0];
        branchingNodes.erase(branchingNodes.begin());

        for (auto &branching : target.branchingCons) {
            if (branching.ctrType == 0) {
                XPRBctr ctr = XPRBnewctr(masterSolver, branching.name.c_str(), XPRB_L);
                XPRBaddterm(ctr, masterLocations[branching.branchingLocationId], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);

            } else if (branching.ctrType == 1) {
                XPRBctr ctr = XPRBnewctr(masterSolver, branching.name.c_str(), XPRB_G);
                XPRBaddterm(ctr, masterLocations[branching.branchingLocationId], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);
            }
        }

        if (useDualSimplexInBranch) {
            XPRBloadmat(masterSolver);
            XPRBloadbasis(target.basis);
//            basis = nullptr;
            XPRBlpoptimise(masterSolver, "d");
//            basis = XPRBsavebasis(masterSolver);
        } else {
            XPRBlpoptimise(masterSolver, "");
        }


        if (XPRBgetlpstat(masterSolver) == XPRB_LP_OPTIMAL) {
            if (XPRBgetobjval(masterSolver) < ub && abs(XPRBgetobjval(masterSolver) - ub) >= INT_GAP) {
                bendersDecomposition(useOptimalityCutInBendersDecomposition, useDualSimplexInBendersDecomposition,
                                     useRoundingHeuristicInBendersDecomposition);
                branching(target);
            }
        }

        for (auto &branching : target.branchingCons) {
            auto branchingCtr = (XPRBctr) XPRBgetbyname(masterSolver, branching.name.c_str(), XPRB_CTR);
            XPRBdelctr(branchingCtr);
        }
        step++;
    }
    cout << "UB = " << ub << endl;

//    cout << "Total Cost " << best.totalCost << endl;
//    for (int i = 1; i <= numFacility; i++) {
//        cout << "facility " << i << " -> " << best.selectedLocations[i] << endl;
//    }
//    for (int j = 1; j <= numCustomer; j++) {
//        int neareastFacility = 0;
//        double nearestDistance = MAX;
//        for (int i = 1; i < numFacility; i++) {
//            if (abs(best.selectedLocations[i] - 1) <= INT_GAP) {
//                if (servingCosts[i][j] < nearestDistance) {
//                    nearestDistance = servingCosts[i][j];
//                    neareastFacility = i;
//                }
//            }
//        }
//        cout << "Customer " << j << " -> " << neareastFacility << endl;
////        solution.totalCost += servingCosts[neareastFacility][j];
////        locationSelected[neareastFacility] = true;
//    }

}