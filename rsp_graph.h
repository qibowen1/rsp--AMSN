#pragma once
#include <vector>
#include <unordered_map>
#include <utility>
#include <limits>

struct RSPGraph {
    std::vector<std::pair<double, double>> nodes;
    std::vector<std::vector<double>> routing_cost;
    std::vector<std::vector<double>> assign_cost;
    int depot = 0;

    void initialize(int size) {
        routing_cost.resize(size, std::vector<double>(size));
        assign_cost.resize(size, std::vector<double>(size));
    }
};

struct Solution {
    std::vector<int> ring;
    std::unordered_map<int, int> assignments;
    std::vector<bool> in_ring;
    double routing_cost = 0;
    double assign_cost = 0;

    double total_cost() const { return routing_cost + assign_cost; }
    void initialize(int size) {
        in_ring.resize(size, false);
        ring.push_back(0);
        in_ring[0] = true;
    }
};

