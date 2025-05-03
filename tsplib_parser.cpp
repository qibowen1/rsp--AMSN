#include "tsplib_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

namespace fs = std::filesystem;
using namespace std;

RSPGraph parseTSPLIB(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "�޷����ļ�: " << filename << endl;
        exit(1);
    }

    RSPGraph graph;
    string line;
    int dimension = 0;
    bool node_section = false;

    while (getline(file, line)) {
        // ��������β�ո�
        line.erase(line.find_last_not_of(" \t") + 1);
        line.erase(0, line.find_first_not_of(" \t"));

        // ͳһתΪСд����
        string lowerLine = line;
        transform(lowerLine.begin(), lowerLine.end(), lowerLine.begin(), ::tolower);

        if (lowerLine.find("dimension") != string::npos) {
            size_t colon = line.find(':');
            if (colon != string::npos) {
                string dimStr = line.substr(colon + 1);
                // ���˷������ַ�
                dimStr.erase(remove_if(dimStr.begin(), dimStr.end(),
                    [](char c) { return !isdigit(c); }), dimStr.end());
                if (!dimStr.empty()) {
                    dimension = stoi(dimStr);
                    graph.initialize(dimension);
                    graph.nodes.resize(dimension);
                }
            }
        }
        else if (lowerLine.find("node_coord_section") != string::npos) {
            node_section = true;
            continue; // ��������
        }
        else if (line.find("EOF") != string::npos) {
            break;
        }

        if (node_section) {
            istringstream iss(line);
            int id;
            double x, y;
            iss >> id >> x >> y;
            id--; // ת��Ϊ0-based����
            graph.nodes[id] = { x, y };
        }
    }

    // ����������
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            double dx = graph.nodes[i].first - graph.nodes[j].first;
            double dy = graph.nodes[i].second - graph.nodes[j].second;
            double dist = hypot(dx, dy);
            // �޸�Ϊ������������������
            int dist_int = static_cast<int>(round(dist));  // �����������ת��[1,3](@ref)

            // ���¾������
            graph.routing_cost[i][j] = dist_int;
            graph.assign_cost[i][j] = dist_int;
        }
    }

    return graph;
}