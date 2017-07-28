//
// Created by baohuaw on 7/28/17.
//

#ifndef CMIP_PROBLEM_H
#define CMIP_PROBLEM_H



#include <map>
#include <fstream>
#include <sstream>

#define MAX DBL_MAX

using namespace std;

int numFacility;
int numCustomer;

double ub = MAX;

map<int, double> openingCosts;
map<int, map<int, double>> servingCosts;


void readFromFile(const string &fileName) {
    ifstream in(fileName);
    if (!in) {
        cout << "Cannot open input file.\n";
        return;
    }

    string str;
    int lineId = 0;
    while (getline(in, str)) {
        istringstream iss(str);
        if (lineId == 0) {
            cout << str << endl;
        } else if (lineId == 1) {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));

            numFacility = stoi(tokens[0]);
            numCustomer = stoi(tokens[1]);
        } else {
            vector<string> tokens;
            copy(istream_iterator<string>(iss), istream_iterator<string>(),
                 back_inserter(tokens));

            int facilityId = stoi(tokens[0]);
            double openingCost = stof(tokens[1]);
            openingCosts[facilityId] = openingCost;

            for (int j = 1; j <= tokens.size() - 2; j++) {
                servingCosts[facilityId][j] = stof(tokens[j - 1 + 2]);
            }
        }

        lineId++;
    }
    in.close();
}


#endif //CMIP_PROBLEM_H
