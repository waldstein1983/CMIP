//
// Created by baohuaw on 8/8/17.
// Elementary Shortest Path Problem, ESPP
//
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

struct Node {
    int id;
    int earliestTime;
    int latestTime;
    int demand;
};

struct Arc {
    Node source;
    Node target;
    double cost;
    int time;
};

struct Label {
    Node node;
    Label preLabel;
    double cost;
    int arrivalTime;
    int curDemand;

    bool operator<(const Label &other) const {
        return arrivalTime < other.arrivalTime;
    }
};

vector<Label> NPS;

const int capacity = 10;

map<Node, vector<Label>> nodeLabels;

map<Node, vector<Arc>> outArcs;

map<int, Node> nodes;

void init() {
    nodes[0] = {0, 0, 0, 0};
    nodes[1] = {1, 6, 14, 0};
    nodes[2] = {2, 9, 12, 0};
    nodes[3] = {3, 8, 12, 0};
    nodes[4] = {4, 9, 35, 0};
    nodes[5] = {5, 15, 25, 0};

    outArcs[nodes[0]].push_back({nodes[0], nodes[1], 3, 8});
    outArcs[nodes[0]].push_back({nodes[0], nodes[2], 5, 5});
    outArcs[nodes[0]].push_back({nodes[0], nodes[3], 2, 12});

    outArcs[nodes[1]].push_back({nodes[1], nodes[4], 13, 4});
    outArcs[nodes[1]].push_back({nodes[1], nodes[5], 6, 6});

    outArcs[nodes[2]].push_back({nodes[2], nodes[4], 8, 2});

    outArcs[nodes[3]].push_back({nodes[3], nodes[4], 16, 4});

    outArcs[nodes[5]].push_back({nodes[5], nodes[4], 3, 7});
}

void shortestPath(Node &source) {
    Label sourceLabel = {source, nullptr, 0, 0, 0};
    nodeLabels[source].push_back(sourceLabel);
    NPS.push_back(sourceLabel);

    while (!NPS.empty()) {
        sort(NPS.begin(), NPS.end());
        Label targetLabel = NPS[0];
        NPS.erase(NPS.begin());

        for (auto &arc : outArcs[targetLabel.node]) {
            if (arc.target.id != targetLabel.node.id &&
                targetLabel.arrivalTime + arc.time <= arc.target.latestTime &&
                targetLabel.curDemand + arc.target.demand <= capacity) {
                if (nodeLabels.count(arc.target) == 0) {
                    nodeLabels[arc.target].push_back({arc.target, targetLabel, targetLabel.cost + arc.cost,
                                                      max(targetLabel.arrivalTime + arc.time, arc.target.earliestTime),
                                                      targetLabel.curDemand + arc.target.demand});
                } else {
                    bool addNewLabel = true;
                    for (auto it = nodeLabels[arc.target].begin();
                         it != nodeLabels[arc.target].end(); /*it++*/) {

                        if (it->cost >= targetLabel.cost + arc.cost &&
                                it->arrivalTime >= targetLabel.arrivalTime + arc.time &&
                                it->curDemand >= targetLabel.curDemand + arc.target.demand){

                            NPS.erase(remove_if(NPS.begin(), NPS.end(), it),NPS.end());
                            it = nodeLabels[arc.target].erase(it);

                        } else{
                            if(it->cost <= targetLabel.cost + arc.cost &&
                                    it->arrivalTime <= targetLabel.arrivalTime + arc.time &&
                                    it->curDemand <= targetLabel.curDemand + arc.target.demand){
                                addNewLabel = false;
                                break;
                            }
                            ++it;
                        }

                    }

                    if(addNewLabel){
                        Label newLabel = {arc.target, targetLabel, targetLabel.cost + arc.cost,
                                          max(targetLabel.arrivalTime + arc.time, arc.target.earliestTime),
                                          targetLabel.curDemand + arc.target.demand};
                        nodeLabels[arc.target].push_back(newLabel);
                        NPS.push_back(newLabel);
                    }


                }
            }
        }
    }
}
