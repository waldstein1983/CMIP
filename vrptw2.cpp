//
// Created by baohuaw on 8/14/17.
//

#define MAX 10000000;
#define INT_GAP 0.00001

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <xprb.h>
#include <iterator>
#include <sstream>
#include <fstream>
#include <chrono>



using namespace std;
using namespace std::chrono;

double ub = MAX;

struct BranchingConstraint {
    string name;
//    int branchingId;
    Path *branchingPath;
    int ctrType;
    int bound;
};


struct BranchingConstraintSet {
//    XPRBbasis basis;
    vector<BranchingConstraint> branchingCons;
};

vector<BranchingConstraintSet> branchingNodes;

struct Node {
    int id;
    double earliestTime;
    double latestTime;
    double serviceTime;
    int demand;
    double cost;
    double x;
    double y;

    bool operator<(const Node &other) const {
        return id < other.id;
    }

    bool operator==(const Node &other) {
        return this == &other;

    }
};

struct Arc {
    Node *source;
    Node *target;
    double cost;
    double time;
//    double originalCost;

    Arc(Node *source, Node *target, double cost, double time) : source(source), target(target), cost(cost),
                                                                time(time) {}
};

struct Label {
    Node *node;
    Label *preLabel;
    double cost;
    double arrivalTime;
    int curDemand;

    bool operator<(const Label *other) const {
        return arrivalTime < other->arrivalTime;
    }

    bool operator==(const Label &other) {
        return this == &other;

    }

    Label(Node *node, Label *preLabel, double cost, double arrivalTime, int curDemand) : node(node), preLabel(preLabel),
                                                                                         cost(cost),
                                                                                         arrivalTime(arrivalTime),
                                                                                         curDemand(curDemand) {}
};

struct Path {
    vector<Node *> content;
    double cost;

    Path(vector<Node *> &content, double cost) : content(content), cost(cost) {}

    bool operator<(const Path *other) const {
        return cost < other->cost;
    }
};

struct InsertingTuple {
    Path *path;
    int pos;
    Node *insertion;
    double score;

    InsertingTuple(Path *path, int pos, Node *insertion, double score) : path(path), pos(pos),
                                                                         insertion(insertion), score(score) {}

    bool operator<(const InsertingTuple &other) const {
        return score < other.score;
    }
};

vector<Path *> pathSet;

vector<Path *> solution;
vector<Label *> NPS;

//vector<Vehicle> vehicles;

const int capacity = 200;

map<Node *, vector<Label *>> nodeLabels;

map<Node *, vector<Arc *>> outArcs;

vector<Arc *> arcs;

Node virtualSource;
const int sourceNodeId = 0;
const int targetNodeId = 101;
Node virtualTarget;
map<int, Node> nodes;
vector<Node *> unhandledTasks;

void computePathCost(Path *path);

string pathToString(Path *path);

//XPRBprob model = XPRBnewprob("vrptw");
//map<Vehicle *, map<Path *, XPRBvar>> x_v_p;
//map<Node *, double> duals;
//map<Node *, XPRBctr> ctrs;

//bool isPathCycle(Path *path) {
//    for (int i = 1; i <= path->content.size() - 1; i++) {
//        for (int j = i + 1; j <= path->content.size() - 1; j++) {
//            if (path->content[j]->id == path->content[i]->id) {
//                return true;
//            }
//        }
//    }
//    return false;
//}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

bool isPathCycle(Label *targetLabel, Node *arcTarget) {
    Label *target = targetLabel;
    while (target->preLabel != nullptr) {
        if (target->node->id == arcTarget->id) {
            return true;
        }
        target = target->preLabel;
    }
    return false;
}

bool checkVehicleTimeWindow(Path *path) {
    double departTime = 0;
    double arriveTime = 0;

    arriveTime = departTime + distance(virtualSource.x, virtualSource.y, path->content[0]->x, path->content[0]->y);
    if (arriveTime < path->content[0]->earliestTime) {
        departTime = path->content[0]->earliestTime + path->content[0]->serviceTime;
    } else if (arriveTime > path->content[0]->latestTime) {
        return false;
    } else {
        departTime = arriveTime + path->content[0]->serviceTime;
    }

    for (int i = 1; i < path->content.size(); i++) {
        arriveTime = departTime + distance(path->content[i - 1]->x, path->content[i - 1]->y, path->content[i]->x,
                                           path->content[i]->y);
        if (arriveTime < path->content[i]->earliestTime) {
            departTime = path->content[i]->earliestTime + path->content[i]->serviceTime;
        } else if (arriveTime > path->content[i]->latestTime) {
            return false;
        } else {
            departTime = arriveTime + path->content[i]->serviceTime;
        }
    }


    arriveTime =
            departTime +
            distance(path->content[path->content.size() - 1]->x, path->content[path->content.size() - 1]->y,
                     virtualTarget.x, virtualTarget.y);
//    if (arriveTime < virtualTarget.earliestTime) {
//        arriveTime = virtualTarget.earliestTime;
//    }

    if (arriveTime > virtualTarget.latestTime) {
        return false;
    }

    return true;
}

bool isPathExisting(Path *p) {
    for (auto &path : pathSet) {
        if (path->content.size() == p->content.size()) {
            bool samePath = true;
            for (int i = 0; i < path->content.size(); i++) {
                if (path->content[i]->id != p->content[i]->id) {
                    samePath = false;
                    break;
                }
            }
            if (samePath) {
                return true;
            }
        }
    }

    return false;
}

Path *shortestPath(Node *source, bool includeCyclePath) {
    auto *sourceLabel = new Label(source, nullptr, 0, 0, 0);
    nodeLabels[source].push_back(sourceLabel);
    NPS.push_back(sourceLabel);

    while (!NPS.empty()) {

        sort(NPS.begin(), NPS.end(),
             [](const Label *a, const Label *b) -> bool {
                 return a->arrivalTime < b->arrivalTime;
             });

//        sort(NPS.begin(), NPS.end(), comparePtrToLabel);
        Label *targetLabel = NPS[0];
        NPS.erase(NPS.begin());
//
//        Label *targetLabel = nullptr;
//        double minArriveTime = MAX;
//        for (auto &label : NPS) {
//            if (label->arrivalTime < minArriveTime) {
//                minArriveTime = label->arrivalTime;
//                targetLabel = label;
//            }
//        }
//        for (auto it = NPS.begin(); it != NPS.end();) {
//            if (*it == targetLabel) {
//                it = NPS.erase(it);
//            } else {
//                ++it;
//            }
//        }


        for (auto &arc : outArcs[targetLabel->node]) {
            if (!includeCyclePath && isPathCycle(targetLabel, arc->target))
                continue;

            if (arc->target->id != targetLabel->node->id &&
                targetLabel->arrivalTime + targetLabel->node->serviceTime + arc->time <= arc->target->latestTime &&
                targetLabel->curDemand + arc->target->demand <= capacity) {
                if (nodeLabels.count(arc->target) == 0 ||
                    (nodeLabels.count(arc->target) != 0 && nodeLabels[arc->target].empty())) {
                    Label *newLabel = new Label(arc->target, targetLabel,
                                                targetLabel->cost + arc->cost + arc->target->cost,
                                                max(targetLabel->arrivalTime + targetLabel->node->serviceTime +
                                                    arc->time, arc->target->earliestTime),
                                                targetLabel->curDemand + arc->target->demand);
                    nodeLabels[arc->target].push_back(newLabel);
                    NPS.push_back(newLabel);
                } else {
                    bool addNewLabel = true;
                    for (auto it = nodeLabels[arc->target].begin(); it != nodeLabels[arc->target].end();) {
                        if ((*it)->cost >= targetLabel->cost + arc->cost + arc->target->cost &&
                            (*it)->arrivalTime >=
                            targetLabel->arrivalTime + targetLabel->node->serviceTime + arc->time &&
                            (*it)->curDemand >= targetLabel->curDemand + arc->target->demand) {
                            NPS.erase(std::remove(NPS.begin(), NPS.end(), *it), NPS.end());

//                            for (auto it2 = NPS.begin(); it2 != NPS.end();) {
//                                if ((*it2)->node->id == (*it)->node->id &&
//                                    (*it2)->preLabel == (*it)->preLabel &&
//                                    (*it2)->cost == (*it)->cost &&
//                                    (*it2)->arrivalTime == (*it)->arrivalTime &&
//                                    (*it2)->curDemand == (*it)->curDemand) {
//                                    it2 = NPS.erase(it2);
//                                } else {
//                                    ++it2;
//                                }
//                            }
                            it = nodeLabels[arc->target].erase(it);
                        } else {
                            if ((*it)->cost <= targetLabel->cost + arc->cost + arc->target->cost &&
                                (*it)->arrivalTime <=
                                targetLabel->arrivalTime + targetLabel->node->serviceTime + arc->time &&
                                (*it)->curDemand <= targetLabel->curDemand + arc->target->demand) {
                                addNewLabel = false;
                                break;
                            }
                            ++it;
                        }
                    }

                    if (addNewLabel) {
                        Label *newLabel = new Label(arc->target, targetLabel,
                                                    targetLabel->cost + arc->cost + arc->target->cost,
                                                    max(targetLabel->arrivalTime + targetLabel->node->serviceTime +
                                                        arc->time, arc->target->earliestTime),
                                                    targetLabel->curDemand + arc->target->demand);
                        nodeLabels[arc->target].push_back(newLabel);
                        NPS.push_back(newLabel);
                    }
                }
            }
        }
    }

    vector<Path *> paths;
    int pathId = 0;
    for (auto &node : nodeLabels) {
        if (node.first->id != targetNodeId)
            continue;
        for (auto &label : node.second) {
            if (label->cost >= 0)
                continue;

            vector<Node *> content;
            Label *targetLabel = label;
            if (targetLabel->arrivalTime > targetLabel->node->latestTime) {
                cout << "Node " << targetLabel->node->id << " violate time window" << endl;
                return nullptr;
            }
            string path;
            while (targetLabel->preLabel != nullptr) {
                path = " - " + path;
                path = to_string(targetLabel->preLabel->node->id) + path;
                if (targetLabel->preLabel->node->id == sourceNodeId) {
//                    content.insert(content.begin(), &virtualSource);
                } else {
                    content.insert(content.begin(), &nodes[targetLabel->preLabel->node->id]);
                }

                targetLabel = targetLabel->preLabel;
                if (targetLabel->arrivalTime > targetLabel->node->latestTime) {
                    cout << "Node " << targetLabel->node->id << " violate time window" << endl;
                    return nullptr;
                }
            }
            auto *p = new Path(content, (int) label->cost);

            if (!p->content.empty() && !checkVehicleTimeWindow(p)) {
                cout << "Violate time window " << endl;
            }
            pathId++;
            paths.push_back(p);


        }
    }


    for (auto itr = nodeLabels.begin(); itr != nodeLabels.end(); itr++) {
        for (auto it = itr->second.begin(); it != itr->second.end(); ++it) {
            delete (*it);
        }
        itr->second.clear();
    }
    nodeLabels.clear();


    Path *shortest = nullptr;


    sort(paths.begin(), paths.end(),
         [](const Path *a, const Path *b) -> bool {
             return a->cost < b->cost;
         });


    shortest = paths[0];
//    double minCost = MAX;
//    for (auto &path : paths) {
//        if (path->cost < minCost && !isPathExisting(path)) {
//            minCost = path->cost;
//            shortest = path;
//        }
//    }

    if (shortest == nullptr)
        return nullptr;
//    sort(paths.begin(), paths.end());
//    string path;
//    for (auto &pathNode : paths[0]->content) {
//        path += to_string(pathNode->id) + "->";
//    }
////    path = path + " -> " + to_string(targetNodeId);
//    cout << "least reduced cost path: " << path << endl;
//    Path* shortest = paths[0];

    auto *target = new Path(shortest->content, shortest->cost);
    for (auto it = paths.begin() + 1; it != paths.end(); ++it) {
        delete (*it);
    }
    paths.clear();
    return target;
}


Path *shortestPath(Node *source, bool includeCyclePath, vector<Node *> blockNodes) {
    auto *sourceLabel = new Label(source, nullptr, 0, 0, 0);
    nodeLabels[source].push_back(sourceLabel);
    NPS.push_back(sourceLabel);

    while (!NPS.empty()) {
        Label *targetLabel = nullptr;
        double minArriveTime = MAX;
        for (auto &label : NPS) {
            if (label->arrivalTime < minArriveTime) {
                minArriveTime = label->arrivalTime;
                targetLabel = label;
            }
        }
        for (auto it = NPS.begin(); it != NPS.end();) {
            if (*it == targetLabel) {
                it = NPS.erase(it);
            } else {
                ++it;
            }
        }

        for (auto &arc : outArcs[targetLabel->node]) {

            if (!(find(blockNodes.begin(), blockNodes.end(), arc->target) == blockNodes.end()))
                continue;

            if (!includeCyclePath && isPathCycle(targetLabel, arc->target))
                continue;

            if (arc->target->id != targetLabel->node->id &&
                targetLabel->arrivalTime + targetLabel->node->serviceTime + arc->time <= arc->target->latestTime &&
                targetLabel->curDemand + arc->target->demand <= capacity) {
                if (nodeLabels.count(arc->target) == 0 ||
                    (nodeLabels.count(arc->target) != 0 && nodeLabels[arc->target].empty())) {
                    Label *newLabel = new Label(arc->target, targetLabel,
                                                targetLabel->cost + arc->cost,
                                                max(targetLabel->arrivalTime + targetLabel->node->serviceTime +
                                                    arc->time, arc->target->earliestTime),
                                                targetLabel->curDemand + arc->target->demand);
                    nodeLabels[arc->target].push_back(newLabel);
                    NPS.push_back(newLabel);
                } else {
                    bool addNewLabel = true;
                    for (auto it = nodeLabels[arc->target].begin(); it != nodeLabels[arc->target].end();) {
                        if ((*it)->cost >= targetLabel->cost + arc->cost &&
                            (*it)->arrivalTime >=
                            targetLabel->arrivalTime + targetLabel->node->serviceTime + arc->time &&
                            (*it)->curDemand >= targetLabel->curDemand + arc->target->demand) {

                            for (auto it2 = NPS.begin(); it2 != NPS.end();) {
                                if ((*it2)->node->id == (*it)->node->id &&
                                    (*it2)->preLabel == (*it)->preLabel &&
                                    (*it2)->cost == (*it)->cost &&
                                    (*it2)->arrivalTime == (*it)->arrivalTime &&
                                    (*it2)->curDemand == (*it)->curDemand) {
                                    it2 = NPS.erase(it2);
                                } else {
                                    ++it2;
                                }
                            }
                            it = nodeLabels[arc->target].erase(it);
                        } else {
                            if ((*it)->cost <= targetLabel->cost + arc->cost &&
                                (*it)->arrivalTime <=
                                targetLabel->arrivalTime + targetLabel->node->serviceTime + arc->time &&
                                (*it)->curDemand <= targetLabel->curDemand + arc->target->demand) {
                                addNewLabel = false;
                                break;
                            }
                            ++it;
                        }
                    }

                    if (addNewLabel) {
                        Label *newLabel = new Label(arc->target, targetLabel,
                                                    targetLabel->cost + arc->cost,
                                                    max(targetLabel->arrivalTime + targetLabel->node->serviceTime +
                                                        arc->time, arc->target->earliestTime),
                                                    targetLabel->curDemand + arc->target->demand);
                        nodeLabels[arc->target].push_back(newLabel);
                        NPS.push_back(newLabel);
                    }
                }
            }
        }
    }

//    Path sp;
    vector<Path *> paths;
//    Path* shortest = nullptr;
//    Path shortest(0, );
//    int minCost = MAX;
    int pathId = 0;
    for (auto &node : nodeLabels) {
        if (node.first->id != targetNodeId)
            continue;
        for (auto &label : node.second) {
            vector<Node *> content;
            Label *targetLabel = label;
            if (targetLabel->arrivalTime > targetLabel->node->latestTime) {
                cout << "Node " << targetLabel->node->id << " violate time window" << endl;
                return nullptr;
            }
//            string path = to_string(targetLabel->node->id);
            string path;
//            content.insert(content.begin(), &nodes[targetLabel->node->id]);
//            content.insert(content.begin(), &virtualTarget);

            while (targetLabel->preLabel != nullptr) {
                path = " - " + path;
                path = to_string(targetLabel->preLabel->node->id) + path;
                if (targetLabel->preLabel->node->id == sourceNodeId) {
//                    content.insert(content.begin(), &virtualSource);
                } else {
                    content.insert(content.begin(), &nodes[targetLabel->preLabel->node->id]);
                }

                targetLabel = targetLabel->preLabel;
                if (targetLabel->arrivalTime > targetLabel->node->latestTime) {
                    cout << "Node " << targetLabel->node->id << " violate time window" << endl;
                    return nullptr;
                }
            }
//            if (label->cost < 0)
            cout << "Path from " << source->id << " to " << label->node->id << "   -> " << path << "   cost: "
                 << label->cost << endl;

//            if(label->cost < minCost){
//                Path p(pathId, content, (int)label->cost);
//                shortest = &p;
////
//            }
            Path *p = new Path(content, label->cost);

            if (!p->content.empty() && !checkVehicleTimeWindow(p)) {
                cout << "Violate time window " << endl;
            }
            pathId++;
            paths.push_back(p);


        }
    }


    for (auto itr = nodeLabels.begin(); itr != nodeLabels.end(); itr++) {
        for (auto it = itr->second.begin(); it != itr->second.end(); ++it) {
            delete (*it);
        }
        itr->second.clear();
    }
    nodeLabels.clear();


    Path *shortest = nullptr;
    double minCost = MAX;
    for (auto &path : paths) {
        if (path->cost < minCost) {
            minCost = path->cost;
            shortest = path;
        }
    }
    Path *target = new Path(shortest->content, shortest->cost);
    for (auto it = paths.begin() + 1; it != paths.end(); ++it) {
        delete (*it);
    }
    paths.clear();
    return target;
}


void readFromFile(const string &fileName) {
    ifstream in(fileName);
    if (!in) {
        cout << "Cannot open input file.\n";
        return;
    }

    string str;
    int lineId = 1;
    while (getline(in, str)) {
        istringstream iss(str);

//        if(str == " ")
//            continue;
        if (lineId == 1) {
            cout << str << endl;
        } else if (lineId < 10) {
            if (lineId < 7 && !(str == " "))
                cout << str << endl;
        } else if (lineId == 10) {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));
//            int customerId = stoi(tokens[0]);
            double x = stod(tokens[1]);
            double y = stod(tokens[2]);
            int demand = stoi(tokens[3]);
            double readyTime = stod(tokens[4]);
            double dueTime = stod(tokens[5]);
            double serviceTime = stod(tokens[6]);
            virtualSource = {sourceNodeId, readyTime, dueTime, serviceTime, demand, 0, x, y};

            virtualTarget = {targetNodeId, readyTime, dueTime, serviceTime, demand, 0, x, y};

        } else {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));
            int customerId = stoi(tokens[0]);
            double x = stod(tokens[1]);
            double y = stod(tokens[2]);
            int demand = stoi(tokens[3]);
            double readyTime = stod(tokens[4]);
            double dueTime = stod(tokens[5]);
            double serviceTime = stod(tokens[6]);
            nodes[customerId] = {customerId, readyTime, dueTime, serviceTime, demand, 0, x, y};
        }
        lineId++;
    }
//    targetNodeId = nodes.rbegin()->first + 1;
    in.close();
}


void initArcs() {
    for (auto &node : nodes) {
        Arc *arc = new Arc(&virtualSource, &node.second,
                           distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y),
                           distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y));
        outArcs[&virtualSource].push_back(arc);
        arcs.push_back(arc);
    }
//
//    for (auto &node : nodes) {
//        Arc *arc = new Arc(&node.second, &virtualSource,
//                           distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y),
//                           distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y));
//        outArcs[&node.second].push_back(arc);
//        arcs.push_back(arc);
//    }

//    outArcs[&virtualSource].push_back(
//            {&virtualSource, &virtualTarget, distance(virtualSource.x, virtualSource.y, virtualTarget.x, virtualTarget.y),
//             distance(virtualSource.x, virtualSource.y, virtualTarget.x, virtualTarget.y)});

    for (auto &node : nodes) {
        Arc *arc = new Arc(&node.second, &virtualTarget,
                           distance(virtualTarget.x, virtualTarget.y, node.second.x, node.second.y),
                           distance(virtualTarget.x, virtualTarget.y, node.second.x, node.second.y));
        outArcs[&node.second].push_back(arc);
        arcs.push_back(arc);
    }

    for (auto &source : nodes) {
        for (auto &target : nodes) {
            Arc *arc = new Arc(&source.second, &target.second,
                               distance(source.second.x, source.second.y, target.second.x,
                                        target.second.y),
                               distance(source.second.x, source.second.y, target.second.x,
                                        target.second.y));
            outArcs[&source.second].push_back(arc);
            arcs.push_back(arc);
        }
    }
}

double detourDistance(Path *path, int pos, Node *insertion) {
    double prev_x, prev_y, next_x, next_y;
    if (pos == 0) {
        prev_x = virtualSource.x;
        prev_y = virtualSource.y;
    } else {
        prev_x = path->content[pos - 1]->x;
        prev_y = path->content[pos - 1]->y;
    }

    if (pos == path->content.size()) {
        next_x = virtualTarget.x;
        next_y = virtualTarget.y;
    } else {
        next_x = path->content[pos]->x;
        next_y = path->content[pos]->y;
    }

    double detour = distance(prev_x, prev_y, insertion->x, insertion->y);
    detour += distance(insertion->x, insertion->y, next_x, next_y);
    detour -= distance(prev_x, prev_y, next_x, next_y);
    return detour;
}

bool checkVehicleCapacity(Path *path, Node *insertion) {
    int total = 0;
    for (auto &node : path->content) {
        total += node->demand;
    }
    return total + insertion->demand <= capacity;
}

bool checkVehicleCapacity(Path *path) {
    int total = 0;
    for (auto &node : path->content) {
        total += node->demand;
    }
    return total <= capacity;
}

bool checkTimeWindow(Path *path, Node *insertion, int pos) {
    vector<Node *> temp;
    for (auto &node : path->content) {
        temp.push_back(node);
    }

    temp.insert(temp.begin() + pos, insertion);

    double departTime = 0;
    double arriveTime;

    arriveTime = departTime + distance(virtualSource.x, virtualSource.y, temp[0]->x, temp[0]->y);
    if (arriveTime < temp[0]->earliestTime) {
        arriveTime = temp[0]->earliestTime;
    }

    if (arriveTime > temp[0]->latestTime) {
        return false;
    }

    departTime = arriveTime + temp[0]->serviceTime;

    for (int i = 1; i < temp.size(); i++) {
        arriveTime = departTime + distance(temp[i - 1]->x, temp[i - 1]->y, temp[i]->x, temp[i]->y);
        if (arriveTime < temp[i]->earliestTime) {
            arriveTime = temp[i]->earliestTime;
        }
        if (arriveTime > temp[i]->latestTime) {
            return false;
        }
        departTime = arriveTime + temp[i]->serviceTime;
    }

    arriveTime =
            departTime + distance(temp[temp.size() - 1]->x, temp[temp.size() - 1]->y, virtualTarget.x, virtualTarget.y);
    if (arriveTime < virtualTarget.earliestTime) {
        arriveTime = virtualTarget.earliestTime;
    }

    if (arriveTime > virtualTarget.latestTime) {
        return false;
    }
    return true;
}

//Vehicle *findEmptyVehicle(Node *insertion) {
//    for (auto &vehicle : vehicles) {
//        if (!vehicle.path->content.empty()) {
//            continue;
//        }
//        if (!checkVehicleCapacity(&vehicle, insertion)) {
//            continue;
//        }
//
//        if (!checkVehicleTimeWindow(&vehicle, insertion, 0)) {
//            continue;
//        }
//        return &vehicle;
//    }
//    return nullptr;
//}

InsertingTuple *findNextInsertingTuple(Path *path) {
    vector<InsertingTuple *> tuples;
    for (auto &task : unhandledTasks) {
        if (!checkVehicleCapacity(path, task)) {
            continue;
        }

        for (int pos = 0; pos < path->content.size() + 1; pos++) {
            double detour = detourDistance(path, pos, task);
            double score;
            score = -detour;

//            score -= task->latestTime - task->earliestTime;
            score += -task->cost * 1;

//            double critical = 1 - ((double) (task->latestTime - task->earliestTime) / task->earliestTime);
//
//            score += critical * 7;

//            score += distance(virtualSource.x, virtualSource.y, task->x, task->y);
            if (!checkTimeWindow(path, task, pos)) {
                continue;
            }

            tuples.push_back(new InsertingTuple(path, pos, task, score));
        }
    }

    if (!tuples.empty()) {
//        sort(tuples.begin(), tuples.end());
        InsertingTuple *maxScoreTuple = nullptr;

        sort(tuples.begin(), tuples.end(),
             [](const InsertingTuple *a, const InsertingTuple *b) -> bool {
                 return a->score > b->score;
             });

        maxScoreTuple = tuples[0];

//        double maxScore = -MAX;
//        for (auto &tuple : tuples) {
//            if (tuple->score > maxScore) {
//                maxScore = tuple->score;
//                maxScoreTuple = tuple;
//            }
//        }

        if (maxScoreTuple == nullptr)
            return nullptr;

        auto *target = new InsertingTuple(maxScoreTuple->path, maxScoreTuple->pos,
                                          maxScoreTuple->insertion,
                                          maxScoreTuple->score);
        for (auto it = tuples.begin(); it != tuples.end(); ++it) {
            delete *it;
        }
        tuples.clear();
//        auto maxScoreIndex = static_cast<int>(tuples.size() - 1);
        return target;
    }
    return nullptr;
}

struct RemovingNode {
    Path *path;
    int pos;
    double distance;

    RemovingNode(Path *path, int pos, double distance) : path(path), pos(pos), distance(distance) {}

};

RemovingNode *findMaxDetourNode(vector<Path *> &paths) {
    RemovingNode *targetRemoving = nullptr;
    vector<RemovingNode *> removingNodes;
    double maxDetour = 0;
    for (auto &path : paths) {
        if (path->content.size() == 1) {
            Node *prev = &virtualSource;
            Node *target = path->content[0];
            Node *next = &virtualTarget;

            double detour = distance(prev->x, prev->y, target->x, target->y);
            detour += distance(target->x, target->y, next->x, next->y);
            detour -= distance(prev->x, prev->y, next->x, next->y);

            removingNodes.push_back(new RemovingNode(path, 0, detour));
        } else {
            Node *prev = &virtualSource;
            Node *target = path->content[0];
            Node *next = path->content[1];

            double detour = distance(prev->x, prev->y, target->x, target->y);
            detour += distance(target->x, target->y, next->x, next->y);
            detour -= distance(prev->x, prev->y, next->x, next->y);

            removingNodes.push_back(new RemovingNode(path, 0, detour));

            for (int i = 1; i < path->content.size() - 1; i++) {
                prev = path->content[i - 1];
                target = path->content[i];
                next = path->content[1 + 1];

                double detour = distance(prev->x, prev->y, target->x, target->y);
                detour += distance(target->x, target->y, next->x, next->y);
                detour -= distance(prev->x, prev->y, next->x, next->y);

                removingNodes.push_back(new RemovingNode(path, i, detour));
            }

            prev = path->content[path->content.size() - 2];
            target = path->content[path->content.size() - 1];
            next = &virtualTarget;

            detour = distance(prev->x, prev->y, target->x, target->y);
            detour += distance(target->x, target->y, next->x, next->y);
            detour -= distance(prev->x, prev->y, next->x, next->y);

            removingNodes.push_back(new RemovingNode(path, static_cast<int>(path->content.size() - 1), detour));
        }

    }

//    Path *shortest = nullptr;
//    double minCost = MAX;

    sort(removingNodes.begin(), removingNodes.end(),
         [](const RemovingNode *a, const RemovingNode *b) -> bool {
             return a->distance > b->distance;
         });

    targetRemoving = removingNodes[0];
//    for (auto &removingNode : removingNodes) {
//        if (removingNode->distance > maxDetour) {
//            maxDetour = removingNode->distance;
//            targetRemoving = removingNode;
//        }
//    }

    auto *maxDetourNode = new RemovingNode(targetRemoving->path, targetRemoving->pos, targetRemoving->distance);
//    Path *target = new Path(shortest->content, shortest->cost);
    for (auto it = removingNodes.begin() + 1; it != removingNodes.end(); ++it) {
        delete (*it);
    }
    removingNodes.clear();

    return maxDetourNode;
}

void randomRemovingSearch(vector<Path *> &paths, int removingSize) {
    vector<Node *> removingNodes;
    for (auto &path : paths) {
        for (auto &node : path->content) {

        }
    }
}

void maxDetourRemovingSearch(vector<Path *> &paths) {
    for (int i = 0; i < 30; i++) {
        double totalCost = 0;
        for (auto &path : paths) {
            computePathCost(path);
            totalCost += path->cost;
        }

        RemovingNode *maxDetour = findMaxDetourNode(paths);

        Node *node = maxDetour->path->content[maxDetour->pos];
        maxDetour->path->content.erase(maxDetour->path->content.begin() + maxDetour->pos);

        vector<InsertingTuple *> tuples;
        for (auto &path : paths) {
            if (!checkVehicleCapacity(path, node)) {
                continue;
            }

            for (int pos = 0; pos < path->content.size() + 1; pos++) {
                double detour = detourDistance(path, pos, node);
                double score;
                score = -detour;

//            score -= task->latestTime - task->earliestTime;
//                score += -node->cost * 1;

//            double critical = 1 - ((double) (task->latestTime - task->earliestTime) / task->earliestTime);
//
//            score += critical * 7;

//                score += distance(virtualSource.x, virtualSource.y, node->x, node->y);
                if (!checkTimeWindow(path, node, pos)) {
                    continue;
                }

                tuples.push_back(new InsertingTuple(path, pos, node, score));
            }
        }


        if (!tuples.empty()) {
            InsertingTuple *maxScoreTuple = nullptr;
            sort(tuples.begin(), tuples.end(),
                 [](const InsertingTuple *a, const InsertingTuple *b) -> bool {
                     return a->score > b->score;
                 });
            maxScoreTuple = tuples[0];
//            double maxScore = -MAX;
//            for (auto &tuple : tuples) {
//                if (tuple->score > maxScore) {
//                    maxScore = tuple->score;
//                    maxScoreTuple = tuple;
//                }
//            }

            if (maxScoreTuple == nullptr)
                return;

            maxScoreTuple->path->content.insert(maxScoreTuple->path->content.begin() + maxScoreTuple->pos,
                                                maxScoreTuple->insertion);
            for (auto it = tuples.begin(); it != tuples.end(); ++it) {
                delete *it;
            }
            tuples.clear();


            double totalCost = 0;
            for (auto &path : paths) {
                computePathCost(path);
                totalCost += path->cost;
            }
//            cout << endl;
        }
    }
}

//void remove()
void computePathCost(Path *path) {
    double pathCost = 0;
    pathCost += distance(virtualSource.x, virtualSource.y, path->content[0]->x, path->content[0]->y);
    for (int i = 1; i < path->content.size(); i++) {
        pathCost += distance(path->content[i - 1]->x, path->content[i - 1]->y,
                             path->content[i]->x, path->content[i]->y);
    }
    pathCost += distance(virtualTarget.x, virtualTarget.y, path->content[path->content.size() - 1]->x,
                         path->content[path->content.size() - 1]->y);
    path->cost = pathCost;
}

void
buildPath(vector<Path *> &newPaths, bool considerDetourDistance, bool considerTimeWindowWidth, bool considerNodeCost,
          bool considerDistanceFromDepot) {
    double maxScore = -MAX;
    Node *maxScoreTask = nullptr;
    for (auto &task : unhandledTasks) {
        double detourDuration = distance(virtualSource.x, virtualSource.y, task->x, task->y);
        double score = 0;
        score -= detourDuration;
//        score -= task->latestTime - task->earliestTime;
        score += -task->cost * 1;
//        score -= task->latestTime - task->earliestTime;
//        double critical = 1 - ((double) (task->latestTime - task->earliestTime) / task->earliestTime);
//        score += critical * 7;
        if (maxScore < score) {
            maxScore = score;
            maxScoreTask = task;
        }
    }


//    Vehicle *vehicle = findEmptyVehicle(maxScoreTask);
//    if (vehicle == nullptr)
//        return;
    vector<Node *> content;
    content.push_back(maxScoreTask);
    auto *path = new Path(content, 0);

    for (auto it = unhandledTasks.begin(); it != unhandledTasks.end();) {
        if ((*it)->id == maxScoreTask->id) {
            it = unhandledTasks.erase(it);
        } else {
            ++it;
        }
    }

    InsertingTuple *insertingTuple = findNextInsertingTuple(path);

    while (insertingTuple != nullptr) {
        insertingTuple->path->content.insert(
                insertingTuple->path->content.begin() + insertingTuple->pos,
                insertingTuple->insertion);

        for (auto it = unhandledTasks.begin(); it != unhandledTasks.end();) {
            if ((*it)->id == insertingTuple->insertion->id) {
                it = unhandledTasks.erase(it);
            } else {
                ++it;
            }
        }

        insertingTuple = findNextInsertingTuple(path);
    }

    if (isPathExisting(path)) {
        delete path;
        return;
    }
    computePathCost(path);
    newPaths.push_back(path);
}

void buildPathByHeuristic(vector<Path *> &newPaths, bool considerDetourDistance, bool considerTimeWindowWidth,
                          bool considerNodeCost,
                          bool considerDistanceFromDepot) {
    for (auto &node : nodes) {
        unhandledTasks.push_back(&node.second);
    }
    while (!unhandledTasks.empty()) {
        buildPath(newPaths, considerDetourDistance, considerTimeWindowWidth, considerNodeCost,
                  considerDistanceFromDepot);
    }
}

XPRBprob model = XPRBnewprob("vrptw");
map<Path *, XPRBvar> x;
map<Node *, XPRBctr> ctrs;
XPRBctr obj;
XPRBbasis basis;

void initModel() {
    obj = XPRBnewctr(model, "Obj", XPRB_N);
    XPRBsetobj(model, obj);
    XPRBsetsense(model, XPRB_MINIM);
}

void buildMathModel(vector<Path *> &newPaths) {

    if (x.size() != 0) {
        XPRBloadmat(model);
        XPRBloadbasis(basis);
        basis = nullptr;
    }

    for (int i = 0; i < newPaths.size(); i++) {
        x[newPaths[i]] = XPRBnewvar(model, XPRB_PL, XPRBnewname("x_%d", pathSet.size() + i), 0, 1);
    }

    for (auto &path : pathSet) {
        XPRBsetterm(obj, x[path], path->cost);
    }

    for (auto &path : newPaths) {
        XPRBaddterm(obj, x[path], path->cost);
    }

    XPRBsetobj(model, obj);
//    XPRBprintctr(obj);


    for (auto &node : nodes) {
//        if()
        if (ctrs.count(&node.second) == 0) {
            XPRBctr fulfill = XPRBnewctr(model, "Demand Fulfillment Once", XPRB_E);
            for (auto &path : newPaths) {

                if (find(path->content.begin(), path->content.end(), &(node.second)) != path->content.end()) {
                    XPRBaddterm(fulfill, x[path], 1);
                }
            }
            XPRBaddterm(fulfill, nullptr, 1);
            ctrs[&node.second] = fulfill;
        } else {
            for (auto &path : newPaths) {

                if (find(path->content.begin(), path->content.end(), &(node.second)) != path->content.end()) {
                    XPRBaddterm(ctrs[&(node.second)], x[path], 1);
                }
            }
        }
    }

//    cout << "Solve Math Model.... " << endl;
    XPRBsetmsglevel(model, 1);
    XPRBlpoptimise(model, "");

    basis = XPRBsavebasis(model);

    solution.clear();
    for (auto &var : x) {
        if (abs(XPRBgetsol(var.second) - 1) <= 0.00001) {
            solution.push_back(var.first);
        }
//        cout << XPRBgetvarname(var.second) << "     " << XPRBgetsol(var.second) << "    " << pathToString(var.first) <<  "   cost " << var.first->cost << endl;
//
    }

    cout << "Optimal Cost " << XPRBgetobjval(model) << endl;

//    for(auto &arc : arcs){
//        if(arc->target == &virtualTarget)
//            continue;
//        arc->originalCost = arc->cost;
//        arc->cost = arc->cost - (int) XPRBgetdual(ctrs[(arc->target)]);
//    }

    for (auto &node : nodes) {
        node.second.cost = 0;

        node.second.cost = -XPRBgetdual(ctrs[&node.second]);
    }
}

void buildSimplePath(vector<Path *> &newPaths) {
    for (auto &node : nodes) {
        vector<Node *> content;
        content.push_back(&node.second);
        auto *p = new Path(content, 0);
        computePathCost(p);
        newPaths.push_back(p);
    }

}

string pathToString(Path *path) {
    string pathStr = to_string(virtualSource.id);
    for (auto &node : path->content) {
        pathStr += " -> " + to_string(node->id);
    }
    pathStr += " -> " + to_string(virtualTarget.id);
    return pathStr;
}

double computePathReducedCost(Path *path) {
    double cost = 0;
    cost += distance(virtualSource.x, virtualSource.y, path->content[0]->x, path->content[0]->y);
    cost += path->content[0]->cost;
    for (int i = 1; i < path->content.size(); i++) {
        cost += distance(path->content[i - 1]->x, path->content[i - 1]->y, path->content[i]->x, path->content[i]->y);
        cost += path->content[i]->cost;
    }

    cost += distance(virtualTarget.x, virtualTarget.y, path->content[path->content.size() - 1]->x,
                     path->content[path->content.size() - 1]->y);

    return cost;
}

bool pathInclude(Path *path1, Path *path2) {
    int target = -1;
    for (int i = 0; i < path1->content.size(); i++) {
        if (find(path2->content.begin(), path2->content.end(), path1->content[i]) != path2->content.end()) {
            target = i;
        } else {
            return false;
        }
    }
}

//void removePositiveReducedCostPaths() {
//    for (auto it = pathSet.begin(); it != pathSet.end();) {
//
////        if (computePathReducedCost(*it) >= 0 && (*it)->content.size() != 1) {
////            it = pathSet.erase(it);
////        } else {
////            ++it;
////        }
//    }
//}

//void buildPathBasedOnMinReducedCostPath()

double computeSolutionCost(vector<Path *> &paths) {
    double totalCost = 0;
    for (auto &path : paths) {
        computePathCost(path);
        totalCost += path->cost;
    }
    return totalCost;
}

void solveRelaxModel() {


}

bool isSolutionInteger(vector<Path *> &paths) {
    for (auto &var : x) {
        if (abs(XPRBgetsol(var.second) - round(XPRBgetsol(var.second))) <= INT_GAP) {
            continue;
        }
        return false;
    }
    return true;
}



void buildBranchingConstraintSet(BranchingConstraintSet &targetSet, Path* path, int ctrType) {
    BranchingConstraintSet branchingSet = {{}};
//    branchingSet.basis = XPRBsavebasis(masterSolver);
    for (auto &ctr : targetSet.branchingCons) {
        branchingSet.branchingCons.push_back(ctr);
    }

    int bound = (int) XPRBgetsol(x[path]);
    if (ctrType == 0) {
        string ctrName = "x(";
        ctrName += pathToString(path);
        ctrName += ") <= ";
        ctrName += to_string(bound);
        BranchingConstraint left = {ctrName, path, 0, bound};
        branchingSet.branchingCons.push_back(left);
    } else {
        string ctrName = "x(";
        ctrName += pathToString(path);
        ctrName += ") >= ";
        ctrName += to_string(bound + 1);
        BranchingConstraint right = {ctrName, path, 1, bound + 1};
        branchingSet.branchingCons.push_back(right);
    }

    branchingNodes.insert(branchingNodes.begin(), branchingSet);

}

void branching(BranchingConstraintSet &set) {
    Path* targetBranchingPath = nullptr;
    double gapToHalf = MAX;
    for (auto &var : x) {
        if (abs(XPRBgetsol(var.second) - round(XPRBgetsol(var.second))) <= INT_GAP) {
            continue;
        }
        double fractional = XPRBgetsol(var.second) - (int) XPRBgetsol(var.second);
        if (abs(fractional - 0.5) < gapToHalf) {
            gapToHalf = abs(fractional - 0.5);
            targetBranchingPath = var.first;
        }
    }

    if (targetBranchingPath == nullptr) {
        return;
    }

    buildBranchingConstraintSet(set, targetBranchingPath, 0);
    buildBranchingConstraintSet(set, targetBranchingPath, 1);
}

void pricing() {
    vector<Path *> newPaths;

    buildPathByHeuristic(newPaths, true, true, true, true);

    buildSimplePath(newPaths);
    buildMathModel(newPaths);

    maxDetourRemovingSearch(solution);

    for (auto &path : newPaths) {
        pathSet.push_back(path);
    }

    newPaths.clear();

    if (!isSolutionInteger(solution)) {
        BranchingConstraintSet target = {{}};
        branching(target);
        return;
    }

    int step = 1;
    while (step < 10000) {
        Path *minReducePath = shortestPath(&virtualSource, false);
        if (minReducePath == nullptr)
            break;
        computePathCost(minReducePath);
        cout << "Step " << step << " , " << pathToString(minReducePath) << " , reduced cost "
             << computePathReducedCost(minReducePath) << "  , numOfPath " << pathSet.size() << endl;

        newPaths.push_back(minReducePath);
        unhandledTasks.clear();
        for (auto &node : nodes) {
            if (find(minReducePath->content.begin(), minReducePath->content.end(), &node.second) ==
                minReducePath->content.end()) {
                unhandledTasks.push_back(&node.second);
            }
        }
        while (!unhandledTasks.empty()) {
            buildPath(newPaths, true, true, false, true);
        }

        buildMathModel(newPaths);

        maxDetourRemovingSearch(solution);

        for (auto &path : newPaths) {
            pathSet.push_back(path);
        }
//
        newPaths.clear();

        if (!isSolutionInteger(solution)) {
            BranchingConstraintSet target = {{}};
            branching(target);
            break;
        }

        step++;
    }
}

void branchAndPrice() {
    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/vrptw/solomon_100/C103.txt";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/vrptw/s-cvrptw/C1_2_1.TXT";
    readFromFile(fileName);
    initArcs();
    initModel();

    pricing();


    while (!branchingNodes.empty()) {
//        milliseconds end = duration_cast<milliseconds>(
//                system_clock::now().time_since_epoch()
//        );
//        if (step % 1 == 0) {
//            cout << ub << "   node size " << branchingNodes.size() << " step = " << step << "   Time "
//                 << (end.count() - start.count()) << endl;
//
//        }
//        if ((end.count() - start.count()) / 1000 > 720000) {
//            cout << "Terminate due to time limit of 3600 sec" << endl;
//            break;
//        }
        BranchingConstraintSet target = branchingNodes[0];
        branchingNodes.erase(branchingNodes.begin());

        for (auto &branching : target.branchingCons) {
            if (branching.ctrType == 0) {
                XPRBctr ctr = XPRBnewctr(model, branching.name.c_str(), XPRB_L);
                XPRBaddterm(ctr, x[branching.branchingPath], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);

            } else if (branching.ctrType == 1) {
                XPRBctr ctr = XPRBnewctr(model, branching.name.c_str(), XPRB_G);
                XPRBaddterm(ctr, x[branching.branchingPath], 1);
                XPRBaddterm(ctr, nullptr, branching.bound);
            }
        }

        XPRBlpoptimise(model, "");

//        if (useDualSimplexInBranch) {
//            XPRBloadmat(masterSolver);
//            XPRBloadbasis(target.basis);
////            basis = nullptr;
//            XPRBlpoptimise(masterSolver, "d");
////            basis = XPRBsavebasis(masterSolver);
//        } else {
//
//        }


        if (XPRBgetlpstat(model) == XPRB_LP_OPTIMAL) {
            if (XPRBgetobjval(model) < ub && abs(XPRBgetobjval(model) - ub) >= INT_GAP) {
//                bendersDecomposition(useOptimalityCutInBendersDecomposition, useDualSimplexInBendersDecomposition,
//                                     useRoundingHeuristicInBendersDecomposition);
                branching(target);
            }
        }

        for (auto &branching : target.branchingCons) {
            auto branchingCtr = (XPRBctr) XPRBgetbyname(model, branching.name.c_str(), XPRB_CTR);
            XPRBdelctr(branchingCtr);
        }
//        step++;
    }


}


int main() {

    milliseconds start = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/vrptw/solomon_100/C103.txt";
//    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/vrptw/s-cvrptw/C1_2_1.TXT";
    readFromFile(fileName);
    initArcs();
    initModel();

    vector<Path *> newPaths;

//    buildPathByHeuristic(newPaths);

    buildPathByHeuristic(newPaths, true, true, true, true);



//    double totalCost = 0;
//    for (auto &path : newPaths) {
//        computePathCost(path);
//        totalCost += path->cost;
//    }

    buildSimplePath(newPaths);


    buildMathModel(newPaths);

    maxDetourRemovingSearch(solution);

//    double solutionCost = computeSolutionCost(solution);

//    for(auto &path : solution){
//        cout << pathToString(path) <<  "   cost " << path->cost << endl;
//    }
    for (auto &path : newPaths) {
        pathSet.push_back(path);
    }

    newPaths.clear();


    int step = 1;
    while (step < 10000) {


//        vector<Node *> blockNodes;
        Path *minReducePath = shortestPath(&virtualSource, false);
//        if(!checkVehicleTimeWindow(minReducePath)){
//            cout << "Min Reduced Cost Path violates time window" << endl;
//            break;
//        }
        if (minReducePath == nullptr)
            break;
//    Path *minReducePath = shortestPath(&virtualSource, true, blockNodes);
        computePathCost(minReducePath);
        cout << "Step " << step << " , " << pathToString(minReducePath) << " , reduced cost "
             << computePathReducedCost(minReducePath) << "  , numOfPath " << pathSet.size() << endl;

//        vector<Path*> newPaths;
        newPaths.push_back(minReducePath);


        unhandledTasks.clear();
        for (auto &node : nodes) {
            if (find(minReducePath->content.begin(), minReducePath->content.end(), &node.second) ==
                minReducePath->content.end()) {
                unhandledTasks.push_back(&node.second);
            }
        }
        while (!unhandledTasks.empty()) {
            buildPath(newPaths, true, true, false, true);
        }

//        vector<Path *> heuristicPaths;
////
//        buildPathByHeuristic(heuristicPaths);
//        neighborSearch(heuristicPaths);
//
//        double totalCost = 0;
//        for (auto &path : heuristicPaths) {
//            computePathCost(path);
//            totalCost += path->cost;
//        }

//        for (auto &path : heuristicPaths) {
//            newPaths.push_back(path);
//        }


        buildMathModel(newPaths);

        maxDetourRemovingSearch(solution);

//        solutionCost = computeSolutionCost(solution);

        for (auto &path : newPaths) {
            pathSet.push_back(path);
        }
//
        newPaths.clear();

        step++;
    }

    cout << "Optimal Solution" << endl;
    for (auto &path : solution) {
        if (!checkVehicleTimeWindow(path)) {
            cout << pathToString(path) << "   violate time window";
            continue;
        }

        if (!checkVehicleCapacity(path)) {
            cout << pathToString(path) << "   violate capacity";
            continue;
        }
        cout << pathToString(path) << endl;
    }


    XPRBdelprob(model);

    for (auto it = pathSet.begin(); it != pathSet.end(); ++it) {
        delete *it;
    }
    pathSet.clear();

    milliseconds end = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
    );

    cout << "Time " << (end.count() - start.count()) << endl;
    return 0;
}



