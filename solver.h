#pragma once
#include "rsp_graph.h"
#include <memory>
using namespace std;

class AMNSSolver {
private:
	RSPGraph graph;
	int alpha;
	int max_search_iter;
	double RCL_ratio;// RCL比例
	int neighborhood_types;
	Solution best;
	static const int VNS_MAX_LEVEL = 1;


public:
	AMNSSolver(const RSPGraph& g, int a, int max_search_iter, double RCL_ratio, int neighborhood_types);
	void solve(int benchmark_opt, int MAX_ITER);
	Solution getBestSolution() const;

private:
	//初始解
	Solution greedy_construction(int solutionPoolSize) const;

	//localsearch


	Solution variable_neighborhood_VNS(Solution s) const;

	Solution shaking(Solution s, int k) const;

	Solution local_search_VND(Solution s, bool fastmodel) const;//VDN

	Solution local_search_FULL_VND(Solution s, bool fast_model) const;//完全遍历

	//遍历领域所有解
	Solution exhaustive_add_node(Solution s) const;

	Solution exhaustive_drop_node(Solution s) const;

	Solution exhaustive_add_drop(Solution s) const;

	Solution exhaustive_two_opt(Solution s) const;

	//领域随机解

	Solution generate_random_add_drop_neighbors(Solution s) const;

	Solution generate_random_2opt_neighbors(Solution s) const;

	void Random_drop_one(Solution& s) const;

	Solution Random_drop_one2(Solution s) const;

	void Random_add_one(Solution& s) const;

	Solution Random_add_one2(Solution s) const;

	//tools

	int find_closest_ring_node(const Solution& s, int u) const;

	void enhanced_full_two_opt(Solution& s) const;

	void randomized_two_opt(Solution& s) const;

	int random_unused_node(const Solution& s) const;

	void insert_node(Solution& s, int v) const;

	void evaluate(Solution& s) const;

	void reallocateNoCircle(Solution& s) const;

};