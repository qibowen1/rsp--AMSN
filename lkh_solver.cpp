#include "lkh_wrapper.h"

#include "lkh_wrapper.h"
#include "tsplib_parser.h"
using namespace std;
Solution solve_with_lkh(const RSPGraph& graph, const std::string& input_dir, const std::string& output_dir, const std::string tspfilename) {
    Solution solution;

    // 1. ׼�� TSP �ļ�
    string tsp_file = input_dir +"/" + tspfilename;
    string tour_file = output_dir + "/temp_lkh_output.txt";

    // 2. ���� LKH
    char* args[] = {
        (char*)"LKH",
        (char*)"PROBLEM_FILE", (char*)tsp_file.c_str(),
        (char*)"OUTPUT_TOUR_FILE", (char*)tour_file.c_str(),
        (char*)"RUNS", (char*)"1",
        nullptr  // ��������
    };
    LKH_Solver(6, args);

    // 3. ���� LKH �������ʵ�ִ˺�����
    solution = parse_lkh_tour(tour_file, graph);

    return solution;
}