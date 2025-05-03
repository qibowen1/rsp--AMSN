#pragma once
#include <string>
struct Result {
    std::string filename;
    int alpha;
    double total_cost;
    double routing_cost;
    double assign_cost;
    double duration_s;
	double offset;//与基准的偏差
    std::string ring_path;
};
