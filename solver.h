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
    void solve();
    Solution getBestSolution() const;

private:
	//初始解
    Solution greedy_construction(int solutionPoolSize) const;

	Solution Random_construction() const;

	Solution GVNS(Solution s,int maxIter) const;
	//localsearch
    Solution variable_neighborhood_search(Solution s) const;
	Solution variable_neighborhood_VNS(Solution s) const;
	Solution local_search_v2(Solution s) const;//一次生成四个领域


    Solution shaking(Solution s, int k) const;

    void drop_redistribute(Solution& s, int deleted_node) const;

    Solution local_search_VND(Solution s) const;//VDN

	void generate_promising_neighbors(Solution s, vector<Solution>& H) const;
	
	//领域操作
	vector<Solution> generate_add_neighbors(Solution s) const;

	vector<Solution> generate_add_neighbors_version2(Solution s) const;

	vector<Solution> generate_Dinf_drop_neighbors(Solution s) const;

	vector<Solution> generate_drop_neighbors2(Solution s) const;

	vector<Solution> generate_add_drop_neighbors(Solution s) const;

	vector<Solution> generate_random_add_drop_neighbors(Solution s) const;

	vector<Solution> generate_opt_neighbors(Solution s) const;

	Solution generate_random_2opt_neighbors(Solution s) const;

	//vector<Solution> generate_3opt_neighbors(Solution s) const;

	void Random_drop_one(Solution& s) const;

	Solution Random_drop_one2(Solution s) const;

	Solution Dinf_drop_op(Solution s) const;

	void Random_add_one(Solution& s) const;

	Solution Random_add_one2(Solution s) const;
	int find_closest_ring_node(const Solution& s, int u) const;

	void enhanced_two_opt(Solution& s) const;

	void enhanced_full_two_opt(Solution& s) const;

	void randomized_two_opt(Solution& s) const;

	//void randomized_three_opt(Solution& s) const;

	//double calculate_segment_cost(const vector<int>& path, int start, int end) const;

	int random_unused_node(const Solution& s) const;

	void insert_node(Solution& s, int v) const;

	void evaluate(Solution& s) const;

};