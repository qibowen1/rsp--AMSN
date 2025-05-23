#include"solver.h"
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <queue>
#include <numeric>
#include <limits>
#include <unordered_set>  // 必须包含该头文件
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <iomanip>
#include <ctime>
#include <functional>
using namespace std;
// 在需要使用shuffle的函数内部定义
std::random_device rd;  // 随机种子
std::mt19937 rng(rd()); // 使用Mersenne Twister引擎
const double INF = numeric_limits<double>::max(); // 定义无穷大常量

// 实现AMNSSolver所有方法
/*
针对Ring Star Problem（RSP），​​自适应多阶段邻域搜索算法（Adaptive Multi-phase Neighborhood Search, AMNS）​​，
RSP问题是在图中寻找一个环形路径，使得每个非节点都被分配到环上的某个节点，是的环边的花费和分配边的花费之和最小化。
结合动态邻域选择、以平衡全局探索与局部 开发能力。
*/
/*
a: 3 5 7 9
*/
AMNSSolver::AMNSSolver(const RSPGraph& g, int a, int max_search_iter, double RCL_ratio, int neighborhood_types, string filename)
    : graph(g), alpha(a), max_search_iter(max_search_iter), RCL_ratio(RCL_ratio), neighborhood_types(neighborhood_types), filename(filename){
}

void AMNSSolver::solve(int benchmark_opt, int MAX_ITER) {
	double min_cost = INF;
    glo_benchmark_opt = benchmark_opt;
	build_static_cache(graph.nodes.size()*0.1); // 预计算静态缓存，k=5

    for (int run = 0; run < MAX_ITER; ++run) {
        Solution temp_best;
        Solution current;
        int solutionPoolSize = 1; // 初始解池大小
        current = sample_greedy_construction(solutionPoolSize); // greedy_construction
        temp_best = current;
        cout << "初始解成本: " << current.total_cost() << endl;
        Solution neighbor = variable_neighborhood_VNS(current); // variable_neighborhood_VNS
        evaluate(neighbor);              // 重新计算新解的成本
        double delta = neighbor.total_cost() - current.total_cost();
        if (delta < 0) {
            current = neighbor;
            temp_best = current;
        }
        enhanced_full_two_opt(temp_best); // 执行全局2-opt
        if (temp_best.total_cost() < min_cost) {
			min_cost = temp_best.total_cost();
            best = temp_best;
        }
        cout << "in迭代 #"<<run<<"次，" << " 当前最优解成本: " << temp_best.total_cost()
            << " (路由: " << temp_best.routing_cost
            << ", 分配: " << temp_best.assign_cost << ")" 
            <<"  全局最优解："<<best.total_cost() << endl;
        // 检查是否达到基准
        if (best.total_cost() == benchmark_opt) {
            cout << "★ 达到基准最优 " << benchmark_opt << "，提前终止 ★" << endl;
            break;
        }
    }
    
}

// 预计算静态缓存 附近k个节点 会考虑 邻域
void AMNSSolver::build_static_cache(int k){
    for (int u = 0; u < graph.nodes.size(); ++u) {
        vector<pair<double, int>> all_distances;

        // 计算到所有其他节点的距离
        for (int v = 0; v < graph.nodes.size(); ++v) {
            if (u != v) {
                all_distances.emplace_back(graph.routing_cost[u][v], v);
            }
        }

        // 按距离排序并保留前k个
        sort(all_distances.begin(), all_distances.end());
        for (int i = 0; i < min(k, (int)all_distances.size()); ++i) {
            static_near_cache[u].push_back(all_distances[i].second);
        }
    }
}


Solution AMNSSolver::getBestSolution() const {
    return best; // 返回当前最优解 
}


Solution AMNSSolver::greedy_construction(int solutionPoolSize) const {
    Solution bestInit;
    double bestCost = INF;

    for (int iter = 0; iter < solutionPoolSize; iter++) {
        Solution s;
        s.initialize(graph.nodes.size());
        s.ring.push_back(graph.depot);  // 初始时环只有depot
        s.in_ring[graph.depot] = true;

        // 初始化非环节点集合
        vector<int> non_ring_nodes;
        for (int i = 0; i < graph.nodes.size(); ++i) {
            if (i != graph.depot) {
                non_ring_nodes.push_back(i);
                s.assignments[i] = graph.depot;  // 初始都分配到depot
            }
        }

        bool improved = true;
        double prev_cost = INF;

        // 主循环：随机选择节点并插入最佳位置
        while (improved && !non_ring_nodes.empty()) {
            // 1. 随机选择一个不在环中的节点
            int random_index = rand() % non_ring_nodes.size();
            int v = non_ring_nodes[random_index];

            // 3. 执行插入并更新解
            Solution temp = s;
            insert_node(temp, v);         // 执行节点插入

            // 计算新成本
            evaluate(temp);
            double new_cost = temp.total_cost();

            // 4. 判断是否接受插入
            if (s.ring.size() >= graph.nodes.size() * ((10 - alpha) / 10.0) && new_cost >= prev_cost) {//graph.nodes.size()*((10-alpha)/10.0) 
                improved = false;  // 停止条件
            }
            else {
                s = temp;
                prev_cost = new_cost;
                // 从候选集中移除已插入节点
                non_ring_nodes.erase(non_ring_nodes.begin() + random_index);
            }
        }

        // 保留最佳解
        if (s.total_cost() < bestCost) {
            bestCost = s.total_cost();
            bestInit = s;
        }
    }

    return bestInit;
}

Solution AMNSSolver::sample_greedy_construction(int solutionPoolSize) const {
    Solution s;
    s.initialize(graph.nodes.size());     // 初始化解决方案结构

    vector<int> candidates;               // 候选节点集合（非根节点）
    for (int i = 0; i < graph.nodes.size(); ++i)
        if (i != graph.depot) candidates.push_back(i);

    double prev_cost = INF;  // 记录前一次的总成本
    int stagnation = 0;      // 成本未改善的连续次数            
    const int max_stagnation = 2; // 允许的最大停滞次数

    while (!candidates.empty()) {         // 逐步构建环路径
        vector<pair<int, double>> candidate_scores; // 候选节点评分集合
        // 计算每个候选节点的插入评分
        for (int v : candidates) {
            double min_ring = INF, min_assign = INF;
            // 1. 计算插入环的最小路由成本变化, 最小环增量
            for (int i = 0; i < s.ring.size(); ++i) {
                int next = (i + 1) % s.ring.size();
                double delta = graph.routing_cost[s.ring[i]][v]
                    + graph.routing_cost[v][s.ring[next]]
                        - graph.routing_cost[s.ring[i]][s.ring[next]];
                    min_ring = min(min_ring, delta);
            }

            // 2. 计算该节点能带来的分配成本收益
            double assign_gain = 0;
            for (int u = 0; u < graph.nodes.size(); ++u) {
                if (!s.in_ring[u]) {      // 只考虑当前未分配的节点
                    double new_cost = graph.assign_cost[u][v];
                    double old_cost = graph.assign_cost[u][s.assignments[u]];
					if (old_cost > new_cost) {
						assign_gain += old_cost - new_cost;
					}
                }
            }
            //assign_gain = assign_gain / (graph.nodes.size() - s.ring.size()); // 平均分配收益
            // 综合评分（路由成本变化+分配收益）
            candidate_scores.emplace_back(v, (10 - alpha) * assign_gain-alpha * min_ring);
        }

        // 构建限制候选列表（RCL）
        sort(candidate_scores.begin(), candidate_scores.end(),
            [](auto& a, auto& b) { return a.second > b.second; });

        int rcl_size = (int)(candidate_scores.size() * RCL_ratio); // RCL大小取前80%
        if (rcl_size == 0) break;
        int selected = candidate_scores[rand() % rcl_size].first; // 随机选择RCL中的节点
        insert_node(s, selected);         // 执行节点插入
        candidates.erase(remove(candidates.begin(), candidates.end(), selected), candidates.end());
        evaluate(s);  // 重新计算总成本
        /* 新增终止条件判断 */
        if (s.ring.size() >= 3) { // 环长度达标后开始检查成本变化
            if (s.total_cost() >= prev_cost&&stagnation>1) {
                break; // 成本未改善达到阈值则终止
            }
            else if(s.total_cost() < prev_cost){
                stagnation = 0;       // 重置未改善计数器
                prev_cost = s.total_cost(); // 更新基准成本
            }
            else {
                stagnation++;
            }
        }
        else {
            prev_cost = s.total_cost(); // 环长度不足时继续构建
        }
    }
    return s;
}


Solution AMNSSolver::variable_neighborhood_VNS(Solution s) const {//2
     int k_max; // 动态上限graph.nodes.size()*0.05
    if (graph.nodes.size() < 100) {
        k_max = graph.nodes.size() * 0.18;
    }
    else if(graph.nodes.size() >= 100&&graph.nodes.size()<=200){
		k_max = graph.nodes.size() * 0.05;
    }
    int k = 1;
    int total_iter = 1000;
    int it = 0;
    Solution best_sol = s;
	int no_improve = 0;
	int max_no_improve = 5; // 最大不改善次数
    Solution current = s;
    current = local_search_FULL_VND(current, false);//快速
    cout << "VNS后解成本: " << current.total_cost() << endl;
    do {
        // Shaking阶段：根据当前k值生成扰动解
        Solution shaken_sol = shaking(current, k); //current  shaking(current, k)
        // 使用VND进行局部搜索
        Solution local_optima = local_search_FULL_VND(shaken_sol,false);//时间有点长
        it++;
        // 接受准则
        if (local_optima.total_cost() < best_sol.total_cost()) {
            best_sol = local_optima;
            current = local_optima;
            //k = 1;
			no_improve = 0; // 重置不改善次数
        }
        else {
			no_improve++;
            k ++;  // 循环扰动强度
        }
        if (it > total_iter) {
            cout << "超过最大次数,退出" << endl;
        }
    } while (k <= k_max);
    cout << "shake+vns后的成本: " << current.total_cost() << endl;
    return best_sol;
}


Solution AMNSSolver::local_search_FULL_VND(Solution s, bool fast_model) const {
    const int NEIGHBORHOOD_TYPES = 6;
    int k = 0;
    int it_add = 0, it_drop = 0, it_swap = 0, add_drop = 0, two_opt = 0, tree_opt=0;//add 8ms drop 13ms, it_swap 2ms, add_drop 17 ms, two_opt 4ms
    // 耗时统计（毫秒）
    double time_add = 0, time_drop = 0, time_swap = 0, time_add_drop = 0, time_two_opt = 0, time_tree_opt = 0;
    while (k < NEIGHBORHOOD_TYPES) {
        Solution current_best = s;
        auto start_time = chrono::high_resolution_clock::now(); // 开始计时
        bool improved = false;
            Solution neighbor = s;
            switch (k) {
            case 5: neighbor = exhaustive_add_node(s); it_add++; break;// exhaustive_cache_add_node exhaustive_add_node
            case 3: neighbor = exhaustive_drop_node(s); it_drop++;break;
            case 4:neighbor = fast_swap_ring_nonring(s); it_swap++;break;
            case 2: neighbor = exhaustive_add_drop(s); add_drop++;break; //exhaustive_cache_add_drop  exhaustive_add_drop
            case 0: neighbor = exhaustive_two_opt(s); two_opt++;break;
            case 1: neighbor = randomized_three_opt(s); tree_opt++;break;
            default:break;
            }
            auto end_time = chrono::high_resolution_clock::now(); // 结束计时
            auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
            double elapsed_ms = duration.count() / 1000.0; // 转换为毫秒

            // 记录各操作的耗时
            switch (k) {
            case 5: time_add += elapsed_ms; break;
            case 3: time_drop += elapsed_ms; break;
            case 4: time_swap += elapsed_ms; break;
            case 2: time_add_drop += elapsed_ms; break;
            case 0: time_two_opt += elapsed_ms; break;
            case 1: time_tree_opt += elapsed_ms; break;
            default: break;
            }

            if (neighbor.total_cost() < current_best.total_cost()) {
                current_best = neighbor;
                improved = true;
            }

        if (improved) {
            s = current_best;
            k = 0; // 回到第一个邻域
        }
        else {
            k++;    // 尝试下一个邻域
        }
    }
    // 输出统计结果
    cout << fixed << setprecision(2);
    cout << "=== 文件名:"<<filename<<" alpha:"<<alpha << " 邻域搜索统计 == = " << endl;
    cout << "操作类型\t调用次数\t总耗时(ms)\t平均耗时(ms)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "add_node\t" << it_add << "\t\t" << time_add << "\t\t" << (it_add ? time_add / it_add : 0) << endl;
    cout << "drop_node\t" << it_drop << "\t\t" << time_drop << "\t\t" << (it_drop ? time_drop / it_drop : 0) << endl;
    cout << "swap_node\t" << it_swap << "\t\t" << time_swap << "\t\t" << (it_swap ? time_swap / it_swap : 0) << endl;
    cout << "add_drop\t" << add_drop << "\t\t" << time_add_drop << "\t\t" << (add_drop ? time_add_drop / add_drop : 0) << endl;
    cout << "two_opt\t\t" << two_opt << "\t\t" << time_two_opt << "\t\t" << (two_opt ? time_two_opt / two_opt : 0) << endl;
    cout << "three_opt\t" << tree_opt << "\t\t" << time_tree_opt << "\t\t" << (tree_opt ? time_tree_opt / tree_opt : 0) << endl;
    cout << "================================================" << endl;
    return s;
}

//all
Solution AMNSSolver::exhaustive_add_node(Solution s) const {
    Solution best_solution = s;
	double min_cost = INF;//不一定比s小，后面以概率接受 如果比s小 换成s.total_cost()

    // 遍历所有非环节点
    for (int v = 0; v < graph.nodes.size(); ++v) {

        if (!s.in_ring[v] && v != graph.depot) {
            Solution temp = s;
            insert_node(temp, v);
            evaluate(temp);

            if (temp.total_cost() < min_cost) {
                min_cost = temp.total_cost();
                best_solution = temp;
            }
        }
    }

    return best_solution;
}

Solution AMNSSolver::exhaustive_cache_add_node(Solution s) const {
    Solution best_solution = s;
    double min_cost = INF;

    // 遍历所有非环节点
    for (int v = 0; v < graph.nodes.size(); ++v) {

        if (!s.in_ring[v] && v != graph.depot) {
            Solution temp = s;
            insert_cache_node(temp, v);
            evaluate(temp);

            if (temp.total_cost() < min_cost) {
                min_cost = temp.total_cost();
                best_solution = temp;
            }
        }
    }

    return best_solution;
}

Solution AMNSSolver::exhaustive_drop_node(Solution s) const {
    if (s.ring.size() <= 3) return s; // 保持最小环大小

    Solution best_solution = s;
    double min_cost = INF;

    // 遍历所有非depot的环节点
    for (size_t pos = 1; pos < s.ring.size(); ++pos) {
        int v = s.ring[pos];
		if (v == graph.depot) continue; // 跳过depot
        Solution temp = s;

        // 删除节点
        temp.ring.erase(temp.ring.begin() + pos);
        temp.in_ring[v] = false;

        // 阶段2：构建受影响节点集合
        unordered_set<int> affected_nodes;

        // 收集原分配至删除节点的非环节点
        for (auto& entry : temp.assignments) {//非环点 u 环点 v
            int u = entry.first;
            int v_2 = entry.second;
            if (v_2 == v && !s.in_ring[u]) {
                affected_nodes.insert(u);
            }
        }

        // 添加被删节点自身（需重新分配）
        affected_nodes.insert(v);

        // 阶段3：批量重分配
        for (int u : affected_nodes) {
            int new_assign = find_closest_ring_node(temp, u);
            // 有效性验证（防止无效分配）
            if (new_assign == -1 || !s.in_ring[new_assign]) {
                cerr << "Invalid assignment for node " << u << endl;
                continue;
            }
            // 更新分配关系
            temp.assignments[u] = new_assign;
        }

        evaluate(temp);
        if (temp.total_cost() < min_cost) {
            min_cost = temp.total_cost();
            best_solution = temp;
        }
    }

    return best_solution;
}

Solution AMNSSolver::exhaustive_add_drop(Solution s) const {
    Solution best_solution = s;
	double min_cost = INF;
    // 先尝试所有可能的删除
    Solution after_drop = exhaustive_drop_node(s);
    if (after_drop.total_cost() < min_cost) {//best_solution.total_cost()
        best_solution = after_drop;
    }

    // 再尝试所有可能的添加
    Solution after_add = exhaustive_add_node(after_drop);
    if (after_add.total_cost() < best_solution.total_cost()) {
        best_solution = after_add;
    }

    return best_solution;
}

Solution AMNSSolver::exhaustive_cache_add_drop(Solution s) const {
    Solution best_solution = s;
    double min_cost = INF;
    // 先尝试所有可能的删除
    Solution after_drop = exhaustive_drop_node(s);
    if (after_drop.total_cost() < min_cost) {//best_solution.total_cost()
        best_solution = after_drop;
    }

    // 再尝试所有可能的添加
    Solution after_add = exhaustive_cache_add_node(after_drop);
    if (after_add.total_cost() < best_solution.total_cost()) {
        best_solution = after_add;
    }

    return best_solution;
}

Solution AMNSSolver::exhaustive_two_opt(Solution s) const {
    bool improved;
    int i = 1;
    do {
        improved = false;
        double best_delta = 0;
        int best_i = -1, best_j = -1;

        // 遍历所有可能的2-opt交换
        for (int i = 0; i < s.ring.size() - 1; ++i) {
            for (int j = i + 2; j < s.ring.size(); ++j) {
                // 计算交换前后的成本变化
                double delta = graph.routing_cost[s.ring[i]][s.ring[j]]
                    + graph.routing_cost[s.ring[i + 1]][s.ring[(j + 1) % s.ring.size()]]
                    - graph.routing_cost[s.ring[i]][s.ring[i + 1]]
                    - graph.routing_cost[s.ring[j]][s.ring[(j + 1) % s.ring.size()]];

                if (delta < best_delta) {
                    best_delta = delta;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        // 执行最佳交换
        if (best_delta < -1e-6) { // 考虑浮点误差
            reverse(s.ring.begin() + best_i + 1, s.ring.begin() + best_j + 1);
            s.routing_cost += best_delta;
            improved = true;
        }
    } while (improved);

    return s;
}


Solution AMNSSolver::randomized_three_opt(Solution s) const {
    const int MAX_ATTEMPTS = 50; // 适当增加尝试次数
    int attempts = 0;
    bool improved;
    const int n = s.ring.size();
    Solution temp_best = s;
    if (n < 6) return s; // 至少需要6个节点才能进行3-opt

    // 预计算原始总成本
    double original_total = s.total_cost();

    do {
        improved = false;

        // 1. 随机选择三个不同的边（确保足够间距）
        int i = rand() % (n - 5);
        int j = (i + 2 + rand() % (n - i - 4)) % n;
        int k = (j + 2 + rand() % (n - j - 2)) % n;

        // 确保i < j < k
        if (i > j) swap(i, j);
        if (j > k) swap(j, k);
        if (i > j) swap(i, j);

        // 2. 获取相关节点
        int a = s.ring[i];
        int b = s.ring[(i + 1) % n];
        int c = s.ring[j];
        int d = s.ring[(j + 1) % n];
        int e = s.ring[k];
        int f = s.ring[(k + 1) % n];

        // 3. 计算原始成本（只计算被修改的边）
        double original = alpha * (graph.routing_cost[a][b] +
            graph.routing_cost[c][d] +
            graph.routing_cost[e][f]);

        // 4. 定义7种可能的3-opt情况
        enum ThreeOptCase {
            CASE_0 = 0, // 原始情况
            CASE_1,     // 2-opt (i-j)
            CASE_2,     // 2-opt (j-k)
            CASE_3,     // 2-opt (i-k)
            CASE_4,     // 3-opt (保持方向)
            CASE_5,     // 3-opt (反转中间段)
            CASE_6      // 3-opt (反转最后段)
        };

        // 存储各种情况的成本和路径
        vector<double> costs(7, 0);
        vector<vector<int>> segments(7, s.ring);

        // 情况0: 原始情况（不做任何改变）
        costs[CASE_0] = original;

        // 情况1: 2-opt (i-j)
        reverse(segments[CASE_1].begin() + i + 1, segments[CASE_1].begin() + j + 1);
        costs[CASE_1] = alpha * (graph.routing_cost[a][c] +
            graph.routing_cost[b][e] +
            graph.routing_cost[d][f]);

        // 情况2: 2-opt (j-k)
        reverse(segments[CASE_2].begin() + j + 1, segments[CASE_2].begin() + k + 1);
        costs[CASE_2] = alpha * (graph.routing_cost[a][b] +
            graph.routing_cost[c][e] +
            graph.routing_cost[d][f]);

        // 情况3: 2-opt (i-k)
        reverse(segments[CASE_3].begin() + i + 1, segments[CASE_3].begin() + k + 1);
        costs[CASE_3] = alpha * (graph.routing_cost[a][e] +
            graph.routing_cost[d][b] +
            graph.routing_cost[c][f]);

        // 情况4: 3-opt (保持方向)
        rotate(segments[CASE_4].begin() + i + 1,
            segments[CASE_4].begin() + j + 1,
            segments[CASE_4].begin() + k + 1);
        costs[CASE_4] = alpha * (graph.routing_cost[a][d] +
            graph.routing_cost[e][c] +
            graph.routing_cost[b][f]);

        // 情况5: 3-opt (反转中间段)
        reverse(segments[CASE_5].begin() + i + 1, segments[CASE_5].begin() + j + 1);
        rotate(segments[CASE_5].begin() + i + 1,
            segments[CASE_5].begin() + j + 1,
            segments[CASE_5].begin() + k + 1);
        costs[CASE_5] = alpha * (graph.routing_cost[a][d] +
            graph.routing_cost[e][b] +
            graph.routing_cost[c][f]);

        // 情况6: 3-opt (反转最后段)
        reverse(segments[CASE_6].begin() + j + 1, segments[CASE_6].begin() + k + 1);
        rotate(segments[CASE_6].begin() + i + 1,
            segments[CASE_6].begin() + j + 1,
            segments[CASE_6].begin() + k + 1);
        costs[CASE_6] = alpha * (graph.routing_cost[a][e] +
            graph.routing_cost[d][c] +
            graph.routing_cost[b][f]);

        // 5. 找出最佳改进
        int best_case = CASE_0;
        double best_gain = 0;
        for (int c = 1; c < 7; ++c) {
            double gain = original - costs[c];
            if (gain > best_gain) {
                best_gain = gain;
                best_case = c;
            }
        }

        // 6. 应用改进
        if (best_gain > 1e-6) {
            s.ring = segments[best_case];

            // 更新路由成本（差分更新）
            s.routing_cost += (costs[best_case] - original);

            // 重新计算分配成本（因为环结构改变了）
            reallocateNoCircle(s);
            evaluate(s);

            improved = true;
            attempts = 0;
        }
        else {
            attempts++;
        }

    } while (improved || attempts < MAX_ATTEMPTS);
    return s;
}

// 交换环节点和非环节点（随机抽样版本） 当前邻域最小的 不一定比s小
Solution AMNSSolver::fast_swap_ring_nonring(Solution s) const {
    if (s.ring.size() <= 3) return s;
	Solution best_solution = s;
    double min = INF;
    const int depot = graph.depot;
	const int trials = graph.nodes.size() / 5; // 随机尝试次数 图的大小决定
    double original_cost = s.total_cost();

    for (int t = 0; t < trials; ++t) {
        Solution temp = s;

        // 随机选择环节点（非depot）
        int v_pos = 1 + rand() % (temp.ring.size() - 1);
        int v = temp.ring[v_pos];
		if (v == depot) continue; // 跳过depot
        // 随机选择一个分配到v的非环节点
        vector<int> candidates;
        for (const auto& [u, assigned] : temp.assignments) {
            if (assigned == v && !temp.in_ring[u] && u != depot) {
                candidates.push_back(u);
            }
        }
        if (candidates.empty()) continue;
        int u = candidates[rand() % candidates.size()];

        // 执行交换
        temp.ring[v_pos] = u;
        temp.in_ring[v] = false;
        temp.in_ring[u] = true;

        // 重新分配v和原属于v的节点
        unordered_set<int> affected;
        affected.insert(v);
        for (const auto& [node, assigned] : temp.assignments) {
            if (assigned == v && !temp.in_ring[node]) {
                affected.insert(node);
            }
        }

        // 批量重分配
        for (int node : affected) {
            int new_assign = find_closest_ring_node(temp, node);
            if (new_assign != -1) {
                temp.assignments[node] = new_assign;
            }
        }

        evaluate(temp);

        // 有效性检查
        if (temp.total_cost() < min &&
            find(temp.ring.begin(), temp.ring.end(), depot) != temp.ring.end()) {
			best_solution = temp;
			min = temp.total_cost();
        }
    }
    return best_solution;
}

Solution AMNSSolver::shaking(Solution s, int k) const {
    Solution shaken = s;
    const int min_ring_size = 3;

    // 收集环上和非环节点
    vector<int> ring_nodes;
    vector<int> non_ring_nodes;
    for (int i = 0; i < graph.nodes.size(); ++i) {
        if (shaken.in_ring[i] && i != graph.depot) {
            ring_nodes.push_back(i);
        }
        else if (!shaken.in_ring[i]) {
            non_ring_nodes.push_back(i);
        }
    }

    // 根据扰动强度k决定操作次数
    for (int step = 0; step < k + 1; ++step) {
        if (ring_nodes.empty() || non_ring_nodes.empty()) {
            break; // 无法继续扰动
        }

        // 随机选择两个不同节点
        int node1 = ring_nodes[rand() % ring_nodes.size()];
        int node2 = non_ring_nodes[rand() % non_ring_nodes.size()];

        // 情况1: 两个节点都在环上 - 删除第一个
        if (shaken.in_ring[node1] && shaken.in_ring[node2]) {
            if (shaken.ring.size() > min_ring_size) {
                auto it = find(shaken.ring.begin(), shaken.ring.end(), node1);
                if (it != shaken.ring.end()) {
                    shaken.ring.erase(it);
                    shaken.in_ring[node1] = false;
                    reallocateNoCircle(shaken);
                }
            }
        }
        // 情况2: 两个节点都不在环上 - 添加第二个
        else if (!shaken.in_ring[node1] && !shaken.in_ring[node2]) {
            insert_node(shaken, node2);
        }
        // 情况3: 一个在环一个不在环 - 删除环上的，添加非环的
        else {
            // 确保node1在环上，node2不在环上
            if (!shaken.in_ring[node1]) {
                swap(node1, node2);
            }

            // 删除环上的node1
            if (shaken.ring.size() > min_ring_size) {
                auto it = find(shaken.ring.begin(), shaken.ring.end(), node1);
                if (it != shaken.ring.end()) {
                    shaken.ring.erase(it);
                    shaken.in_ring[node1] = false;
                    reallocateNoCircle(shaken);
                }
            }

            // 添加非环的node2
            insert_node(shaken, node2);
        }
    }

    // 确保解的有效性
    if (shaken.ring.size() < min_ring_size) {
        // 如果扰动导致环太小，恢复原始解
        return s;
    }

    evaluate(shaken);
    return shaken;
}


void AMNSSolver::insert_node(Solution& s, int v) const {
    double min_cost = INF;
    int best_pos = -1;

    // 寻找插入后路由成本最小的位置
    for (int i = 0; i < s.ring.size(); ++i) {
        int j = (i + 1) % s.ring.size();
        double delta = graph.routing_cost[s.ring[i]][v]
            + graph.routing_cost[v][s.ring[j]]
                - graph.routing_cost[s.ring[i]][s.ring[j]];
            if (delta < min_cost) {
                min_cost = delta;
                best_pos = i + 1;                // 插入到i和j之间
            }
    }

    s.ring.insert(s.ring.begin() + best_pos, v); // 插入节点
    s.in_ring[v] = true;                   // 更新标记

    // 更新分配关系（如果新节点更近）
    for (int u = 0; u < graph.nodes.size(); ++u) {
        if (!s.in_ring[u] && graph.assign_cost[u][v] < graph.assign_cost[u][s.assignments[u]]) {
            s.assignments[u] = v;          // 更新分配关系
        }
    }
}

void AMNSSolver::insert_cache_node(Solution& s, int v) const {
    double min_cost = INF;
    int best_pos = -1;
    bool has_valid_candidate = false;

    // 只检查静态缓存中的节点（但必须是当前环中的节点）
    for (int candidate : static_near_cache.at(v)) {
        // 检查候选节点是否在环中
        auto pos = find(s.ring.begin(), s.ring.end(), candidate);
        if (pos == s.ring.end()) continue;

        has_valid_candidate = true;

        // 检查该环节点前后的插入位置
        for (int offset = 0; offset <= 1; ++offset) {
            int insert_pos = (pos - s.ring.begin() + offset) % s.ring.size();
            int prev = s.ring[(insert_pos - 1 + s.ring.size()) % s.ring.size()];
            int next = s.ring[insert_pos % s.ring.size()];

            double delta = graph.routing_cost[prev][v] + graph.routing_cost[v][next]
                - graph.routing_cost[prev][next];

                if (delta < min_cost) {
                    min_cost = delta;
                    best_pos = insert_pos;
                }
        }
    }

    // 如果没有有效候选（缓存中的节点都不在环中）
    if (!has_valid_candidate) {
        return; // 直接不考虑该节点插入
    }

    if (best_pos != -1) {
        s.ring.insert(s.ring.begin() + best_pos, v);
        s.in_ring[v] = true;

        // 更新分配关系
        for (int u = 0; u < graph.nodes.size(); ++u) {
            if (!s.in_ring[u] && graph.assign_cost[u][v] < graph.assign_cost[u][s.assignments[u]]) {
                s.assignments[u] = v;
            }
        }
    }
}

int AMNSSolver::find_closest_ring_node(const Solution& s, int u) const {
    int best = -1;
    double min_cost = INF;
    for (int v : s.ring) {                 // 遍历环中所有节点
        if (graph.assign_cost[u][v] < min_cost) {
            min_cost = graph.assign_cost[u][v];
            best = v;
        }
    }
    return best;
}


void AMNSSolver::enhanced_full_two_opt(Solution& s) const {
    bool improved;
    do {
        improved = false;
        double best_gain = 0;
        int best_i = -1, best_j = -1;
        double pregain = 0;
        // 遍历所有可能的边交换
        for (int i = 0; i < s.ring.size() - 1; ++i) {
            for (int j = i + 2; j < s.ring.size(); ++j) {
                // 计算原始路径的成本
                double original = graph.routing_cost[s.ring[i]][s.ring[i + 1]]
                    + graph.routing_cost[s.ring[j]][s.ring[(j + 1) % s.ring.size()]];
                // 计算交换后的新路径成本
                double modified = graph.routing_cost[s.ring[i]][s.ring[j]]
                    + graph.routing_cost[s.ring[i + 1]][s.ring[(j + 1) % s.ring.size()]];
                // 计算增益
                double gain = original - modified;

                // 找到当前最优的交换
                if (gain > 0 && gain > pregain) {
                    best_gain = gain;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        // 执行最佳交换（只要增益为正）
        if (best_gain > 0) {
            reverse(s.ring.begin() + best_i + 1, s.ring.begin() + best_j + 1);
            improved = true; // 标记本轮有改进，继续下一轮搜索
        }
    } while (improved); // 持续优化直到无法改进
}


void AMNSSolver::reallocateNoCircle(Solution& s) const {
	// 遍历所有非环节点
	for (int u = 0; u < graph.nodes.size(); ++u) {
		if (!s.in_ring[u]) { // 如果是非环节点
			int new_assign = find_closest_ring_node(s, u); // 找到最近的环节点
			if (new_assign != -1 && s.in_ring[new_assign]) { // 确保新分配有效
				s.assignments[u] = new_assign; // 更新分配关系
			}
		}
	}
}

void AMNSSolver::evaluate(Solution& s) const {
    // 计算路由成本（环路径总长度）
    s.routing_cost = 0;
    for (int i = 0; i < s.ring.size(); ++i) {
        int j = (i + 1) % s.ring.size();
        s.routing_cost += alpha * graph.routing_cost[s.ring[i]][s.ring[j]];
    }

    // 计算分配成本（所有非环节点到最近环节点的成本）
    s.assign_cost = 0;
    for (int u = 0; u < graph.nodes.size(); ++u) {
        if (!s.in_ring[u]) {              // 只考虑非环节点
            int v = s.assignments[u];     // 获取当前分配的环节点
            s.assign_cost += (10 - alpha) * graph.assign_cost[u][s.assignments[u]];
        }
    }
}