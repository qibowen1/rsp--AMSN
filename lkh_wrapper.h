#pragma once
#include "rsp_graph.h"
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
	// ���� LKH ������������ LKH.c ���޸�Ϊ�ɵ��õĽӿڣ�
	void LKH_Solver(int argc, char* argv[]);

#ifdef __cplusplus
}
#endif

Solution solve_with_lkh(const RSPGraph& graph, const std::string& input_dir, const std::string& output_dir, const std::string tspfilename);