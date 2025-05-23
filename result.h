#pragma once
#include <string>
struct Result {
    std::string filename;
    int alpha;
    double total_cost;
    double routing_cost;
    double assign_cost;
    double duration_s;
	double offset;//���׼��ƫ��

    // ��ӹ��캯��������������
    Result(std::string f, int a, double t, double r, double as, double d, double o)
        : filename(f), alpha(a),
        total_cost(round(t * 100) / 100),  // ����2λС��
        routing_cost(round(r * 100) / 100),
        assign_cost(round(as * 100) / 100),
        duration_s(round(d * 100) / 100),
        offset(round(o * 100) / 100) {
    }
};
