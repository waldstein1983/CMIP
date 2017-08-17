//
// Created by baohuaw on 8/14/17.
//

#define MAX 100000;

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
    int earliestTime;
    int latestTime;
    int serviceTime;
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
    int cost;
    int time;
};

struct Label {
    Node *node;
    Label *preLabel;
    double cost;
    int arrivalTime;
    int curDemand;

    bool operator<(const Label &other) const {
        return arrivalTime < other.arrivalTime;
    }

    bool operator==(const Label &other) {
        return this == &other;

    }

    Label(Node *node, Label *preLabel, double cost, int arrivalTime, int curDemand) : node(node), preLabel(preLabel),
                                                                                      cost(cost),
                                                                                      arrivalTime(arrivalTime),
                                                                                      curDemand(curDemand) {}
};

struct Path {
//    int id;
//    Node *source;
//    Node *target;
    vector<Node *> content;
    int cost;

    Path(vector<Node *> &content, int cost) : content(content), cost(cost) {}

    bool operator<(const Path *other) const {
        return cost < other->cost;
    }
};

struct Vehicle {
    int id;
    Path *path;
    int capacity;
};

struct InsertingTuple {
    Vehicle *vehicle;
    int pos;
    Node *insertion;
    double score;

    InsertingTuple(Vehicle *vehicle, int pos, Node *insertion, double score) : vehicle(vehicle), pos(pos),
                                                                               insertion(insertion), score(score) {}

    bool operator<(const InsertingTuple &other) const {
        return score < other.score;
    }
};

//map<Vehicle *, vector<Path>> vehiclePaths;
vector<Path *> pathSet;

//map<Node, vector<Path>> allPaths;

vector<Label *> NPS;

vector<Vehicle> vehicles;

const int capacity = 200;

map<Node *, vector<Label *>> nodeLabels;

map<Node *, vector<Arc>> outArcs;

Node virtualSource;
const int sourceNodeId = 0;
const int targetNodeId = 11;
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
    int departTime = 0;
    int arriveTime = 0;

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
        } else if (arriveTime > path->content[0]->latestTime) {
            return false;
        } else {
            departTime = arriveTime + path->content[i]->serviceTime;
        }
    }
    return true;
}

Path *shortestPathWithoutCycle(Node *source, bool includeCyclePath, vector<Node *> blockNodes) {
    auto *sourceLabel = new Label(source, nullptr, 0, 0, 0);
    nodeLabels[source].push_back(sourceLabel);
    NPS.push_back(sourceLabel);

    while (!NPS.empty()) {
//        sort(NPS.begin(), NPS.end());
        Label *targetLabel = nullptr;
        int minArriveTime = MAX;
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

            if (!(find(blockNodes.begin(), blockNodes.end(), arc.target) == blockNodes.end()))
                continue;

            if (!includeCyclePath && isPathCycle(targetLabel, arc.target))
                continue;

            if (arc.target->id != targetLabel->node->id &&
                targetLabel->arrivalTime + targetLabel->node->serviceTime + arc.time <= arc.target->latestTime &&
                targetLabel->curDemand + arc.target->demand <= capacity) {
                if (nodeLabels.count(arc.target) == 0 ||
                    (nodeLabels.count(arc.target) != 0 && nodeLabels[arc.target].empty())) {
                    Label *newLabel = new Label(arc.target, targetLabel,
                                                targetLabel->cost + arc.cost + arc.target->cost,
                                                max(targetLabel->arrivalTime + targetLabel->node->serviceTime +
                                                    arc.time, arc.target->earliestTime),
                                                targetLabel->curDemand + arc.target->demand);
                    nodeLabels[arc.target].push_back(newLabel);
                    NPS.push_back(newLabel);
                } else {
                    bool addNewLabel = true;
                    for (auto it = nodeLabels[arc.target].begin(); it != nodeLabels[arc.target].end();) {
                        if ((*it)->cost >= targetLabel->cost + arc.cost + arc.target->cost &&
                            (*it)->arrivalTime >=
                            targetLabel->arrivalTime + targetLabel->node->serviceTime + arc.time &&
                            (*it)->curDemand >= targetLabel->curDemand + arc.target->demand) {

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
                            it = nodeLabels[arc.target].erase(it);
                        } else {
                            if ((*it)->cost <= targetLabel->cost + arc.cost + arc.target->cost &&
                                (*it)->arrivalTime <=
                                targetLabel->arrivalTime + targetLabel->node->serviceTime + arc.time &&
                                (*it)->curDemand <= targetLabel->curDemand + arc.target->demand) {
                                addNewLabel = false;
                                break;
                            }
                            ++it;
                        }
                    }

                    if (addNewLabel) {
                        Label *newLabel = new Label(arc.target, targetLabel,
                                                    targetLabel->cost + arc.cost + arc.target->cost,
                                                    max(targetLabel->arrivalTime + targetLabel->node->serviceTime +
                                                        arc.time, arc.target->earliestTime),
                                                    targetLabel->curDemand + arc.target->demand);
                        nodeLabels[arc.target].push_back(newLabel);
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
//        if (node.first->id != targetNodeId)
//            continue;
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
    int minCost = MAX;
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
            double x = stof(tokens[1]);
            double y = stof(tokens[2]);
            int demand = stoi(tokens[3]);
            int readyTime = stoi(tokens[4]);
            int dueTime = stoi(tokens[5]);
            int serviceTime = stoi(tokens[6]);
            virtualSource = {sourceNodeId, readyTime, dueTime, serviceTime, demand, 0, x, y};

            virtualTarget = {targetNodeId, readyTime, dueTime, serviceTime, demand, 0, x, y};

        } else {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));
            int customerId = stoi(tokens[0]);
            double x = stof(tokens[1]);
            double y = stof(tokens[2]);
            int demand = stoi(tokens[3]);
            int readyTime = stoi(tokens[4]);
            int dueTime = stoi(tokens[5]);
            int serviceTime = stoi(tokens[6]);
            nodes[customerId] = {customerId, readyTime, dueTime, serviceTime, demand, 0, x, y};
        }
        lineId++;
    }
    in.close();
}



void initArcs() {
    for (auto &node : nodes) {
        outArcs[&virtualSource].push_back(
                {&virtualSource, &node.second, distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y),
                 distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y)});
    }

//    outArcs[&virtualSource].push_back(
//            {&virtualSource, &virtualTarget, distance(virtualSource.x, virtualSource.y, virtualTarget.x, virtualTarget.y),
//             distance(virtualSource.x, virtualSource.y, virtualTarget.x, virtualTarget.y)});

    for (auto &node : nodes) {
        outArcs[&node.second].push_back(
                {&node.second, &virtualTarget, distance(virtualTarget.x, virtualTarget.y, node.second.x, node.second.y),
                 distance(virtualSource.x, virtualSource.y, node.second.x, node.second.y)});
    }

    for (auto &source : nodes) {
        for (auto &target : nodes) {
            outArcs[&source.second].push_back({&source.second, &target.second,
                                               distance(source.second.x, source.second.y, target.second.x,
                                                        target.second.y),
                                               distance(source.second.x, source.second.y, target.second.x,
                                                        target.second.y)});
        }
    }
}

int detourDistance(Vehicle *vehicle, int pos, Node *insertion) {
    double prev_x, prev_y, next_x, next_y;
    if (pos == 0) {
        prev_x = virtualSource.x;
        prev_y = virtualSource.y;
    } else {
        prev_x = vehicle->path->content[pos - 1]->x;
        prev_y = vehicle->path->content[pos - 1]->y;
    }

    if (pos == vehicle->path->content.size()) {
        next_x = virtualTarget.x;
        next_y = virtualTarget.y;
    } else {
        next_x = vehicle->path->content[pos]->x;
        next_y = vehicle->path->content[pos]->y;
    }

    int detour = distance(prev_x, prev_y, insertion->x, insertion->y);
    detour += distance(insertion->x, insertion->y, next_x, next_y);
    detour -= distance(prev_x, prev_y, next_x, next_y);
    return detour;
}

bool checkVehicleCapacity(Vehicle *vehicle, Node *insertion) {
    int total = 0;
    for (auto &node : vehicle->path->content) {
        total += node->demand;
    }
    return total + insertion->demand <= vehicle->capacity;
}

bool checkVehicleTimeWindow(Vehicle *vehicle, Node *insertion, int pos) {
    vector<Node *> temp;
    for (auto &node : vehicle->path->content) {
        temp.push_back(node);
    }

    temp.insert(temp.begin() + pos, insertion);

    int departTime = 0;
    int arriveTime;

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

Vehicle *findEmptyVehicle(Node *insertion) {
    for (auto &vehicle : vehicles) {
        if (!vehicle.path->content.empty()) {
            continue;
        }
        if (!checkVehicleCapacity(&vehicle, insertion)) {
            continue;
        }

        if (!checkVehicleTimeWindow(&vehicle, insertion, 0)) {
            continue;
        }
        return &vehicle;
    }
    return nullptr;
}

InsertingTuple *findNextInsertingTuple(Vehicle *vehicle) {
    vector<InsertingTuple *> tuples;
    for (auto &task : unhandledTasks) {
        if (!checkVehicleCapacity(vehicle, task)) {
            continue;
        }

        for (int pos = 0; pos < vehicle->path->content.size() + 1; pos++) {
            int detour = detourDistance(vehicle, pos, task);
            double score;
            score = -detour;

            double critical = 1 - ((double) (task->latestTime - task->earliestTime) / task->earliestTime);
            score += critical * 7;

            score += distance(virtualSource.x, virtualSource.y, task->x, task->y);
            if (!checkVehicleTimeWindow(vehicle, task, pos)) {
                continue;
            }

            tuples.push_back(new InsertingTuple(vehicle, pos, task, score));
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

        InsertingTuple *target = new InsertingTuple(maxScoreTuple->vehicle, maxScoreTuple->pos,
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

void buildPath() {
    double maxScore = -MAX;
    Node *maxScoreTask = nullptr;
    for (auto &task : unhandledTasks) {
        double detourDuration = distance(virtualSource.x, virtualSource.y, task->x, task->y);
        double score = 0;
        score -= detourDuration;
        score -= task->latestTime - task->earliestTime;
        score += -task->cost;
        double critical = 1 - ((double) (task->latestTime - task->earliestTime) / task->earliestTime);
        score += critical * 7;
        if (maxScore < score) {
            maxScore = score;
            maxScoreTask = task;
        }
    }


    Vehicle *vehicle = findEmptyVehicle(maxScoreTask);
    if (vehicle == nullptr)
        return;
    vehicle->path->content.push_back(maxScoreTask);

    for (auto it = unhandledTasks.begin(); it != unhandledTasks.end();) {
        if ((*it)->id == maxScoreTask->id) {
            it = unhandledTasks.erase(it);
        } else {
            ++it;
        }
    }

    InsertingTuple *insertingTuple = findNextInsertingTuple(vehicle);

    while (insertingTuple != nullptr) {
        insertingTuple->vehicle->path->content.insert(
                insertingTuple->vehicle->path->content.begin() + insertingTuple->pos,
                insertingTuple->insertion);

        for (auto it = unhandledTasks.begin(); it != unhandledTasks.end();) {
            if ((*it)->id == insertingTuple->insertion->id) {
                it = unhandledTasks.erase(it);
            } else {
                ++it;
            }
        }

        insertingTuple = findNextInsertingTuple(vehicle);
    }
}

void initVehicles() {
    for (int i = 1; i <= 100; i++) {
        vector<Node *> nodes;
        vehicles.push_back({i, new Path(nodes, 0), 200});
    }
}

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

void buildPathBasedOnUnhandledTasks() {


    while (!unhandledTasks.empty()) {
        buildPath();
    }

    for (auto &vehicle : vehicles) {
        if (vehicle.path->content.empty())
            continue;

        if (find(pathSet.begin(), pathSet.end(), vehicle.path) != pathSet.end())
            continue;


        cout << "Vehicle " << vehicle.id << endl;
        string path = to_string(virtualSource.id);
        for (int i = 0; i < vehicle.path->content.size(); i++) {
            path = path + " -> " + to_string(vehicle.path->content[i]->id);
        }
        path = path + " -> " + to_string(virtualTarget.id);
        cout << path << endl;

//        vehicle.path->id = static_cast<int>(pathSet.size());
        computePathCost(vehicle.path);
        pathSet.push_back(vehicle.path);

//        vehiclePaths[&vehicle].push_back(vehicle.path);
    }
}

void buildMathModel() {
//    XPRBresetprob(model);
//    x_v_p.clear();
//    ctrs.clear();
    XPRBprob model = XPRBnewprob("vrptw");
    map<Vehicle *, map<Path *, XPRBvar>> x_v_p;
//    map<Node *, double> duals;
    map<Node *, XPRBctr> ctrs;

    for (auto &vehicle : vehicles) {
        for (int i = 0; i < pathSet.size(); i++) {
            x_v_p[&vehicle][pathSet[i]] = XPRBnewvar(model, XPRB_PL, XPRBnewname("x_%d_%d", vehicle.id, i), 0, 1);
        }
    }

    XPRBctr obj = XPRBnewctr(model, "Obj", XPRB_N);
    for (auto &vehicle : vehicles) {
        for (auto &path : pathSet) {
            XPRBaddterm(obj, x_v_p[&vehicle][path], path->cost);
        }
    }
    XPRBsetobj(model, obj);
//    XPRBprintctr(obj);


    for (auto &node : nodes) {
        XPRBctr fulfill = XPRBnewctr(model, "Demand Fulfillment Once", XPRB_E);
        for (auto &vehicle : vehicles) {
            for (auto &path : pathSet) {

                if (find(path->content.begin(), path->content.end(), &(node.second)) != path->content.end()) {
                    XPRBaddterm(fulfill, x_v_p[&vehicle][path], 1);
                }
            }
        }
        XPRBaddterm(fulfill, nullptr, 1);
        ctrs[&node.second] = fulfill;
//        XPRBprintctr(fulfill);
    }
    XPRBsetsense(model, XPRB_MINIM);

    cout << "Solve Math Model.... " << endl;
    XPRBsetmsglevel(model, 4);
    XPRBlpoptimise(model, "");
    for (auto &vehicle : vehicles) {
        for (auto &path : pathSet) {
            if (abs(XPRBgetsol(x_v_p[&vehicle][path]) - 1) <= 0.00001)
                cout << "Vehicle " << vehicle.id << "  ->  "
                     << XPRBgetsol(x_v_p[&vehicle][path])
                     << "   cost " << path->cost << endl;
        }
    }

    cout << "Optimal Cost " << XPRBgetobjval(model) << endl;

    for (auto &node : nodes) {
        node.second.cost = 0;
        node.second.cost = -XPRBgetdual(ctrs[&node.second]);
    }
}


void findMinReducedCostPath() {
//    Path* minReducePath = shortestPath(&virtualSource);
    vector<Node *> blockNodes;


    Path *minReducePath = shortestPathWithoutCycle(&virtualSource, false, blockNodes);
    computePathCost(minReducePath);
    pathSet.push_back(minReducePath);

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


int main() {
    initVehicles();

    string fileName = "/home/local/ANT/baohuaw/CLionProjects/CMIP/data/vrptw/solomon_25/C101.txt";
    readFromFile(fileName);
    initArcs();

//    vector<Node *> blockNodes;
//    Path *minReducePath = shortestPathWithoutCycle(&virtualSource, true, blockNodes);

    for (auto &node : nodes) {
        unhandledTasks.push_back(&node.second);
    }
    buildPathBasedOnUnhandledTasks();
    buildMathModel();

    int step = 1;
    while(step < 10){
        for (auto &node : nodes) {
            unhandledTasks.push_back(&node.second);
        }
        buildPathBasedOnUnhandledTasks();
        buildMathModel();
        step++;
    }


//    findMinReducedCostPath();
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

//    buildMathModel();

    for (auto it = pathSet.begin(); it != pathSet.end(); ++it) {
        delete *it;
    }
    pathSet.clear();
    return 0;
}



