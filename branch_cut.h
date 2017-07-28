//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_BRANCH_CUT_H
#define CMIP_BRANCH_CUT_H

void branchAndCut(bool useOptimalityCutInBendersDecomposition, bool useDualSimplexInBendersDecomposition,
                  bool useRoundingHeuristicInBendersDecomposition,
                  bool useDualSimplexInBranch, bool useInOutStrategyInBendersDecomposition);

#endif //CMIP_BRANCH_CUT_H
