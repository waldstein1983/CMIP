//
// Created by baohuaw on 8/14/17.
//

#define MAX 10000000;

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <xprb.h>
#include <iterator>
#include <sstream>
#include <fstream>


using namespace std;

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
    double originalCost;

    Arc(Node *source, Node *target, double cost, double time) : source(source), target(target), cost(cost), time(time) {}
};

struct Label {
    Node *node;
    Label *preLabel;
    double cost;
    double arrivalTime;
    int curDemand;

    bool operator<(const Label &other) const {
        return arrivalTime < other.arrivalTime;
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
//    int id;
//    Node *source;
//    Node *target;
    vector<Node *> content;
    double cost;

    Path(vector<Node *> &content, double cost) : content(content), cost(cost) {}

    bool operator<(const Path *other) const {
        return cost < other->cost;
    }
};

//struct Vehicle {
//    int id;
//    Path *path;
//    int capacity;
//};

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

//map<Vehicle *, vector<Path>> vehiclePaths;
vector<Path *> pathSet;

//map<Node, vector<Path>> allPaths;

vector<Label *> NPS;

//vector<Vehicle> vehicles;

const int capacity = 200;

map<Node *, vector<Label *>> nodeLabels;

map<Node *, vector<Arc *>> outArcs;

vector<Arc *> arcs;

Node virtualSource;
const int sourceNodeId = 0;
const int targetNodeId = 26;
Node virtualTarget;
map<int, Node> nodes;
vector<Node *> unhandledTasks;

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

int distance(double x1, double y1, double x2, double y2) {
    return (int) sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
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
    return true;
}

bool isPathExisting(Path* p){
    for(auto &path : pathSet){
        if(path->content.size() == p->content.size()){
            bool samePath = true;
            for(int i = 0;i < path->content.size();i++){
                if(path->content[i]->id != p->content[i]->id){
                    samePath = false;
                    break;
                }
            }
            if(samePath){
                return true;
            }
        }
    }

    return false;
}

Path *shortestPathWithoutCycle(Node *source, bool includeCyclePath, vector<Node *> blockNodes) {
    auto *sourceLabel = new Label(source, nullptr, 0, 0, 0);
    nodeLabels[source].push_back(sourceLabel);
    NPS.push_back(sourceLabel);

    while (!NPS.empty()) {
//        sort(NPS.begin(), NPS.end());
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


//        Label *targetLabel = NPS[0];
//        NPS.erase(NPS.begin());

        for (auto &arc : outArcs[targetLabel->node]) {

//            if (!(find(blockNodes.begin(), blockNodes.end(), arc->target) == blockNodes.end()))
//                continue;
//
//            if (!includeCyclePath && isPathCycle(targetLabel, arc->target))
//                continue;

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
//            cout << "Path from " << source->id << " to " << label->node->id << "   -> " << path << "   cost: "
//                 << label->cost << endl;

//            if(label->cost < minCost){
//                Path p(pathId, content, (int)label->cost);
//                shortest = &p;
////
//            }
            Path *p = new Path(content, (int) label->cost);

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
        if (path->cost < minCost && !isPathExisting(path)) {
            minCost = path->cost;
            shortest = path;
        }
    }

    if(shortest == nullptr)
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
            Path *p = new Path(content,  label->cost);

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
//    sort(paths.begin(), paths.end());
//    string path;
//    for (auto &pathNode : paths[0]->content) {
//        path += to_string(pathNode->id) + "->";
//    }
////    path = path + " -> " + to_string(targetNodeId);
//    cout << "least reduced cost path: " << path << endl;
//    Path* shortest = paths[0];

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
        if (lineId == 1) {
            cout << str << endl;
        } else if (lineId < 10) {
            if (lineId < 7)
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
            Arc* arc = new Arc(&source.second, &target.second,
                               distance(source.second.x, source.second.y, target.second.x,
                                        target.second.y),
                               distance(source.second.x, source.second.y, target.second.x,
                                        target.second.y));
            outArcs[&source.second].push_back(arc);
            arcs.push_back(arc);
        }
    }
}

int detourDistance(Path *path, int pos, Node *insertion) {
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

    int detour = distance(prev_x, prev_y, insertion->x, insertion->y);
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
            int detour = detourDistance(path, pos, task);
            double score;
            score = -detour;

//            score -= task->latestTime - task->earliestTime;
            score += -task->cost;

//            double critical = 1 - ((double) (task->latestTime - task->earliestTime) / task->earliestTime);
//
//            score += critical * 7;

            score += distance(virtualSource.x, virtualSource.y, task->x, task->y);
            if (!checkTimeWindow(path, task, pos)) {
                continue;
            }

//            cout << "Score " << score << endl;
            tuples.push_back(new InsertingTuple(path, pos, task, score));
        }
    }

    if (!tuples.empty()) {
//        sort(tuples.begin(), tuples.end());
        InsertingTuple *maxScoreTuple = nullptr;
        double maxScore = -MAX;
        for (auto &tuple : tuples) {
            if (tuple->score > maxScore) {
                maxScore = tuple->score;
                maxScoreTuple = tuple;
            }
        }

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

//void remove()
void computePathCost(Path *path) {
    int pathCost = 0;
    pathCost += distance(virtualSource.x, virtualSource.y, path->content[0]->x, path->content[0]->y);
    for (int i = 1; i < path->content.size(); i++) {
        pathCost += distance(path->content[i - 1]->x, path->content[i - 1]->y,
                             path->content[i]->x, path->content[i]->y);
    }
    pathCost += distance(virtualTarget.x, virtualTarget.y, path->content[path->content.size() - 1]->x,
                         path->content[path->content.size() - 1]->y);
    path->cost = pathCost;
}

void buildPath() {
    double maxScore = -MAX;
    Node *maxScoreTask = nullptr;
    for (auto &task : unhandledTasks) {
        double detourDuration = distance(virtualSource.x, virtualSource.y, task->x, task->y);
        double score = 0;
        score -= detourDuration;
        score -= task->latestTime - task->earliestTime;
        score += -task->cost;
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
    Path *path = new Path(content, 0);
//    vehicle->path->content.push_back(maxScoreTask);

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

    computePathCost(path);
    pathSet.push_back(path);
}

//void initVehicles() {
//    for (int i = 1; i <= 10000; i++) {
//        vector<Node *> nodes;
//        vehicles.push_back({i, new Path(nodes, 0), 200});
//    }
//}



void buildPathBasedOnUnhandledTasks() {


    while (!unhandledTasks.empty()) {
        buildPath();
    }

//    for (auto &vehicle : vehicles) {
//        if (vehicle.path->content.empty())
//            continue;
//
//        if (find(pathSet.begin(), pathSet.end(), vehicle.path) != pathSet.end())
//            continue;
//
//
//        cout << "Vehicle " << vehicle.id << endl;
//        string path = to_string(virtualSource.id);
//        for (int i = 0; i < vehicle.path->content.size(); i++) {
//            path = path + " -> " + to_string(vehicle.path->content[i]->id);
//        }
//        path = path + " -> " + to_string(virtualTarget.id);
//        cout << path << endl;
//
////        vehicle.path->id = static_cast<int>(pathSet.size());
//        computePathCost(vehicle.path);
//        pathSet.push_back(vehicle.path);
//
////        vehiclePaths[&vehicle].push_back(vehicle.path);
//    }
}

void buildMathModel() {
//    XPRBresetprob(model);
//    x_v_p.clear();
//    ctrs.clear();
    XPRBprob model = XPRBnewprob("vrptw");
    map<Path *, XPRBvar> x;
//    map<Node *, double> duals;
    map<Node *, XPRBctr> ctrs;

    for (int i = 0; i < pathSet.size(); i++) {
        x[pathSet[i]] = XPRBnewvar(model, XPRB_PL, XPRBnewname("x_%d", i), 0, 1);
    }

    XPRBctr obj = XPRBnewctr(model, "Obj", XPRB_N);
    for (auto &path : pathSet) {
        XPRBaddterm(obj, x[path], path->cost);
    }

    XPRBsetobj(model, obj);
//    XPRBprintctr(obj);


    for (auto &node : nodes) {
        XPRBctr fulfill = XPRBnewctr(model, "Demand Fulfillment Once", XPRB_E);
        for (auto &path : pathSet) {

            if (find(path->content.begin(), path->content.end(), &(node.second)) != path->content.end()) {
                XPRBaddterm(fulfill, x[path], 1);
            }
        }
        XPRBaddterm(fulfill, nullptr, 1);
        ctrs[&node.second] = fulfill;
//        XPRBprintctr(fulfill);
    }
    XPRBsetsense(model, XPRB_MINIM);

    cout << "Solve Math Model.... " << endl;
    XPRBsetmsglevel(model, 1);
    XPRBlpoptimise(model, "");
    for (auto &path : pathSet) {
        if (abs(XPRBgetsol(x[path]) - 1) <= 0.00001)
            cout << XPRBgetsol(x[path])
                 << "   cost " << path->cost << endl;
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

    XPRBdelprob(model);
}


void buildMathModel2() {
//    XPRBresetprob(model);
//    x_v_p.clear();
//    ctrs.clear();
    XPRBprob model = XPRBnewprob("vrptw");
    map<Path *, XPRBvar> x;
    XPRBvar xd;
    XPRBvar xc;
//    map<Node *, double> duals;
    map<Node *, XPRBctr> ctrs;
    XPRBctr numOfRoute;
    XPRBctr objValue;

    for (int i = 0; i < pathSet.size(); i++) {
        x[pathSet[i]] = XPRBnewvar(model, XPRB_PL, XPRBnewname("p_%d", i), 0, 1);
    }

    xd = XPRBnewvar(model, XPRB_PL, XPRBnewname("x_d"), 0, 100000);

    xc = XPRBnewvar(model, XPRB_PL, XPRBnewname("x_c"), 0, 100000);

    XPRBctr obj = XPRBnewctr(model, "Obj", XPRB_N);
    for (auto &path : pathSet) {
        XPRBaddterm(obj, x[path], path->cost);
    }

    XPRBsetobj(model, obj);
//    XPRBprintctr(obj);


    for (auto &node : nodes) {
        XPRBctr fulfill = XPRBnewctr(model, "Demand Fulfillment Once", XPRB_G);
        for (auto &path : pathSet) {
            if (find(path->content.begin(), path->content.end(), &(node.second)) != path->content.end()) {
                XPRBaddterm(fulfill, x[path], 1);
            }
        }
        XPRBaddterm(fulfill, nullptr, 1);
        ctrs[&node.second] = fulfill;
    }

    numOfRoute = XPRBnewctr(model, "Number of Route", XPRB_E);
    for (auto &path : pathSet) {
        XPRBaddterm(numOfRoute, x[path], 1);
    }
    XPRBaddterm(numOfRoute, xd, -1);
    XPRBaddterm(numOfRoute, nullptr, 0);

    objValue = XPRBnewctr(model, "Obj Value", XPRB_E);
    for (auto &path : pathSet) {
        XPRBaddterm(objValue, x[path], path->cost);
    }
    XPRBaddterm(objValue, xc, -1);
    XPRBaddterm(objValue, nullptr, 0);


    XPRBsetsense(model, XPRB_MINIM);

    cout << "Solve Math Model.... " << endl;
    XPRBsetmsglevel(model, 4);
    XPRBlpoptimise(model, "");
    for (auto &path : pathSet) {
        if (abs(XPRBgetsol(x[path]) - 1) <= 0.00001)
            cout << XPRBgetsol(x[path]) << "   cost " << path->cost << endl;
    }
    cout << "X_d " << XPRBgetsol(xd)  << endl;
    cout << "X_c " << XPRBgetsol(xc)  << endl;


    cout << "Optimal Cost " << XPRBgetobjval(model) << endl;

    double dual_xc = XPRBgetdual (objValue);
    cout << "dual for xc " << dual_xc << endl;
    for(auto &arc : arcs){

//        arc->originalCost = arc->cost;
        arc->cost = arc->time;

        arc->cost = arc->cost * (1 - XPRBgetdual (objValue));

//        if(arc->target == &virtualTarget)
//            continue;
//        arc->cost = arc->cost - XPRBgetdual(ctrs[(arc->target)]);
    }

    for (auto &node : nodes) {
        node.second.cost = 0;
        node.second.cost = -XPRBgetdual(ctrs[&node.second]);

        cout << "dual for node " << XPRBgetdual(ctrs[&node.second]) << endl;
    }

    XPRBdelprob(model);
}


void findMinReducedCostPath() {
//    Path* minReducePath = shortestPath(&virtualSource);
    vector<Node *> blockNodes;


    Path *minReducePath = shortestPathWithoutCycle(&virtualSource, true, blockNodes);
//    Path *minReducePath = shortestPath(&virtualSource, true, blockNodes);
    computePathCost(minReducePath);
    pathSet.push_back(minReducePath);

//    for(auto &arc : arcs){
//        arc->cost = arc->time;
//    }

    unhandledTasks.clear();
//
    for (auto &node : nodes) {
        if (find(minReducePath->content.begin(), minReducePath->content.end(), &node.second) ==
            minReducePath->content.end()) {
            unhandledTasks.push_back(&node.second);
        }
    }
//
    buildPathBasedOnUnhandledTasks();
}

void buildSimplePath(){
    for(auto &node : nodes){
        vector<Node*> content;
        content.push_back(&node.second);
        auto *p = new Path(content,0);
        computePathCost(p);
        pathSet.push_back(p);
    }

}

int main() {

    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/vrptw/solomon_25/C101.txt";
    readFromFile(fileName);
    initArcs();

//    vector<Node *> blockNodes;
//    Path *minReducePath = shortestPathWithoutCycle(&virtualSource, true, blockNodes);

//    for (auto &node : nodes) {
//        unhandledTasks.push_back(&node.second);
//    }
//    buildPathBasedOnUnhandledTasks();
    buildSimplePath();
    buildMathModel();


    int step =1;
    while(step < 100){
        vector<Node *> blockNodes;
        Path *minReducePath = shortestPathWithoutCycle(&virtualSource, true, blockNodes);
//    Path *minReducePath = shortestPath(&virtualSource, true, blockNodes);
        computePathCost(minReducePath);
        pathSet.push_back(minReducePath);
        buildMathModel();
        step++;
    }


//    int step = 1;
//    while (step < 10) {
//        for (auto &node : nodes) {
//            unhandledTasks.push_back(&node.second);
//        }
//        buildPathBasedOnUnhandledTasks();
//        buildMathModel();
//        step++;
//    }


    findMinReducedCostPath();
//
//
//    for (auto &path : pathSet) {
//        string content;
//        for (auto &node : path->content) {
//            content = content + " -> ";
//            content = content + to_string(node->id);
//        }
//        cout << content << endl;
//    }

    buildMathModel2();

    for (auto it = pathSet.begin(); it != pathSet.end(); ++it) {
        delete *it;
    }
    pathSet.clear();
    return 0;
}



