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
max_search_iter:对每个初始解的迭代次数 3 ls
RCL_ratio:rcl比例 0.8
neighborhood_types:邻域类型数量 5 4
*/
AMNSSolver::AMNSSolver(const RSPGraph& g, int a, int max_search_iter, double RCL_ratio, int neighborhood_types)
    : graph(g), alpha(a), max_search_iter(max_search_iter), RCL_ratio(RCL_ratio), neighborhood_types(neighborhood_types) {
}

void AMNSSolver::solve(int benchmark_opt, int MAX_ITER) {
	double min_cost = INF;
	build_static_cache(5); // 预计算静态缓存，k=5
    for (int run = 0; run < MAX_ITER; ++run) {
        Solution temp_best;
        Solution current;
        int solutionPoolSize = 5; // 初始解池大小
        current = greedy_construction(solutionPoolSize); // greedy_construction
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
            if (s.ring.size() >= 3 && new_cost >= prev_cost) {
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

Solution AMNSSolver::variable_neighborhood_VNS(Solution s) const {//2
    const int k_max = 9; // 动态上限graph.nodes.size()*0.05
    int k = 1;
    Solution best_sol = s;
    Solution current = local_search_FULL_VND(s, false);//快速
    cout << "VNS后解成本: " << current.total_cost() << endl;
    do {
        // Shaking阶段：根据当前k值生成扰动解
        Solution shaken_sol = shaking(current, k); //current  shaking(current, k)
        // 使用VND进行局部搜索
        Solution local_optima = local_search_FULL_VND(shaken_sol,false);//时间有点长

        // 接受准则
        if (local_optima.total_cost() < best_sol.total_cost()) {
            best_sol = local_optima;
            current = local_optima;
            k = 1;
        }
        else {
            k ++;  // 循环扰动强度
        }
    } while (k <= k_max);
    cout << "shake+vns后的成本: " << current.total_cost() << endl;
    return best_sol;
}


//方法3 随机变量邻域搜索，使用VND进行局部搜索，直到没有改进为止

//Solution AMNSSolver::local_search_VND(Solution s, bool fast_model) const {
//    int k = 0;
//    const int max_iter = fast_model ? 3 : 5; // 根据模式调整迭代次数
//    const int localsearchMax = fast_model ? graph.nodes.size() / 5 : graph.nodes.size() / 2;
//    do {
//        Solution local_best = s;
//        int choose = rand() % max_iter;
//        for (int i = 0;i < localsearchMax;i++) {
//            vector<Solution> H;
//            Solution temp = s;
//            // 生成当前类型的邻域
//            switch (choose) {
//            case 0: temp = Random_drop_one2(s);break;
//            case 1: temp = Random_add_one2(s);break;
//            case 2: temp = generate_random_add_drop_neighbors(s); break;
//            case 3: temp = generate_random_2opt_neighbors(s); break;
//            default:break;
//            }
//            evaluate(temp);
//            if (temp.total_cost() < local_best.total_cost()) {
//                // 计算当前邻域的解
//                local_best = temp;
//            }
//        }
//        // 如果找到改进解
//        if (local_best.total_cost() < s.total_cost()) {
//            k = 0;
//            s = local_best;
//        }
//        else k++;
//    } while (k < max_iter);  // 只要任一邻域有改进就继续
//
//    return s;
//}

Solution AMNSSolver::local_search_VND(Solution s, bool fast_model) const {
    const int NEIGHBORHOOD_TYPES = 4; // 邻域类型数量
    int k = 0; // 当前邻域索引

    while (k < NEIGHBORHOOD_TYPES) {
        Solution current_best = s;
        bool improved = false;
        const int max_trials = fast_model ? graph.nodes.size() / 4 : graph.nodes.size();

        for (int i = 0; i < max_trials; ++i) {
            Solution neighbor = s;

            // 显式调用不同邻域操作
            switch (k) {
            case 3: neighbor = generate_random_2opt_neighbors(s); break;
            case 0: neighbor = Random_add_one2(s); break;
            case 1: neighbor = Random_drop_one2(s); break;
            case 2: neighbor = generate_random_add_drop_neighbors(s); break;
            }

            evaluate(neighbor);
            if (neighbor.total_cost() < current_best.total_cost()) {
                current_best = neighbor;
                improved = true;
            }
        }

        if (improved) {
            s = current_best;
            k = 0; // 重置到第一个邻域
        }
        else {
            k++;    // 尝试下一个邻域
        }
    }

    return s;
}

Solution AMNSSolver::local_search_FULL_VND(Solution s, bool fast_model) const {
    const int NEIGHBORHOOD_TYPES = 6;
    int k = 0;

    while (k < NEIGHBORHOOD_TYPES) {
        Solution current_best = s;
        bool improved = false;
            Solution neighbor = s;
            switch (k) {
            case 0: neighbor = exhaustive_add_node(s); break;
            case 1: neighbor = exhaustive_drop_node(s); break;
            case 2: neighbor = exhaustive_add_drop(s); break;
			case 3: neighbor = batch_drop_nodes(s,3); break;
            case 4: neighbor = batch_add_nodes(s,3); break;
            case 5: neighbor = exhaustive_two_opt(s); break;
            default:break;
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

    return s;
}

//all
Solution AMNSSolver::exhaustive_add_node(Solution s) const {
    Solution best_solution = s;
    double min_cost = s.total_cost();

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
    double min_cost = s.total_cost();

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
    double min_cost = s.total_cost();

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

    // 先尝试所有可能的删除
    Solution after_drop = exhaustive_drop_node(s);
    if (after_drop.total_cost() < best_solution.total_cost()) {
        best_solution = after_drop;
    }

    // 再尝试所有可能的添加
    Solution after_add = exhaustive_add_node(after_drop);
    if (after_add.total_cost() < best_solution.total_cost()) {
        best_solution = after_add;
    }

    return best_solution;
}

Solution AMNSSolver::exhaustive_two_opt(Solution s) const {
    bool improved;
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

//批量LS
Solution AMNSSolver::batch_add_nodes(Solution s, int batch_size) const {
    vector<int> candidates;
    // 筛选最有潜力的非环节点（按分配成本降序） 分配成本高的可能加入环能更好
    for (int u = 0; u < graph.nodes.size(); ++u) {
        if (!s.in_ring[u] && u != graph.depot) {
            candidates.push_back(u);
        }
    }
    sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        return graph.assign_cost[a][s.assignments[a]] > graph.assign_cost[b][s.assignments[b]];
        });

    // 尝试添加前K个节点
    for (int i = 0; i < min(batch_size, (int)candidates.size()); ++i) {
        int v = candidates[i];
        double min_insert_cost = INF;
        int best_pos = -1;

        // 寻找全局最优插入位置
        for (int j = 0; j < s.ring.size(); ++j) {
            int next = (j + 1) % s.ring.size();
            double delta = alpha * (graph.routing_cost[s.ring[j]][v]
                + graph.routing_cost[v][s.ring[next]]
                    - graph.routing_cost[s.ring[j]][s.ring[next]]);
            if (delta < min_insert_cost) {
                min_insert_cost = delta;
                best_pos = j + 1;
            }
        }

        // 执行插入并更新分配
        if (best_pos != -1) {
            s.ring.insert(s.ring.begin() + best_pos, v);
            s.in_ring[v] = true;
        }

    }
	reallocateNoCircle(s); // 重新分配
	evaluate(s); // 重新计算成本
    return s;
}

Solution AMNSSolver::batch_drop_nodes(Solution s, int batch_size) const {
    if (s.ring.size() <= 3 + batch_size) return s; // 保持最小环大小

    vector<pair<double, int>> ring_scores;
    // 计算环内节点效用分数（越低越应删除）
    for (int i = 0; i < s.ring.size(); ++i) {
        int v = s.ring[i];
        if (v == graph.depot) continue;

        double routing_contribution = alpha * (
            graph.routing_cost[s.ring[(i - 1 + s.ring.size()) % s.ring.size()]][v] +
            graph.routing_cost[v][s.ring[(i + 1) % s.ring.size()]]
            );
        double assign_contribution = (10 - alpha) * accumulate_assigned_cost(s, v);
        ring_scores.emplace_back(assign_contribution-routing_contribution, v);//效用分数 = (10-α) × 分配成本贡献-α × 路由成本贡献,效用分数越小说明该节点路由成本高 但是分配成本低,删除可能更好
    }

    // 选择效用最低的K个节点
    sort(ring_scores.begin(), ring_scores.end());
    vector<int> to_remove;
    for (int i = 0; i < min(batch_size, (int)ring_scores.size()); ++i) {
        to_remove.push_back(ring_scores[i].second);
    }

    // 批量删除并重分配
    for (int v : to_remove) {
        auto it = find(s.ring.begin(), s.ring.end(), v);
        if (it != s.ring.end()) {
            s.ring.erase(it);
            s.in_ring[v] = false;
        }
    } 
	reallocateNoCircle(s); // 重新分配
	evaluate(s); // 重新计算成本
    return s;
}


// 交换环节点和非环节点（随机抽样版本） 当前邻域最小的 不一定比s小
Solution AMNSSolver::fast_swap_ring_nonring(Solution s) const {
    if (s.ring.size() <= 3) return s;
	Solution best_solution = s;
    double min = INF;
    const int depot = graph.depot;
	const int trials = graph.nodes.size() / 10; // 随机尝试次数 图的大小决定
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

// 计算分配到该环节点的所有非环节点的总分配成本
double AMNSSolver::accumulate_assigned_cost(const Solution& s, int ring_node) const {
    double total = 0.0;
    for (const auto& [u, assigned_v] : s.assignments) {
        if (assigned_v == ring_node && !s.in_ring[u]) { // 只统计非环节点的分配
            total += graph.assign_cost[u][ring_node];
        }
    }
    return total;
}

//random
Solution AMNSSolver::generate_random_add_drop_neighbors(Solution s) const {
    Solution neighbor = s;
    Solution res = s;
    int temp = s.total_cost();
    Random_drop_one(neighbor);
    evaluate(neighbor);
    if (neighbor.total_cost() < temp) {
		temp = neighbor.total_cost();
        res = neighbor;
    }
    Random_add_one(neighbor);
    evaluate(neighbor);
    if (neighbor.total_cost() < temp) {
        res = neighbor;
    }

    return res;
}



Solution AMNSSolver::generate_random_2opt_neighbors(Solution s) const {
    randomized_two_opt(s);
    evaluate(s);
    return s;
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

void AMNSSolver::Random_add_one(Solution& s) const {
    int v = random_unused_node(s);         // 随机选择不在环中的节点
    if (v == -1) return;                  // 无可用节点时退出
    insert_node(s, v);                     // 执行插入操作
}

Solution AMNSSolver::Random_add_one2(Solution s) const {
    int v = random_unused_node(s);         // 随机选择不在环中的节点
    if (v == -1) return s;                  // 无可用节点时退出
    insert_node(s, v);                     // 执行插入操作
    return s;
}

void AMNSSolver::Random_drop_one(Solution& s) const {
    if (s.ring.size() <= 3) return;

    // 阶段1：删除节点
    int pos = 1 + rand() % (s.ring.size() - 1);
    int deleted_node = s.ring[pos];
    if (deleted_node == 0) return;
    s.ring.erase(s.ring.begin() + pos);
    s.in_ring[deleted_node] = false;

    // 阶段2：构建受影响节点集合
    unordered_set<int> affected_nodes;

    // 收集原分配至删除节点的非环节点
    for (auto& entry : s.assignments) {//非环点 u 环点 v
        int u = entry.first;
        int v = entry.second;
        if (v == deleted_node) {
            affected_nodes.insert(u);
        }
    }

    // 添加被删节点自身（需重新分配）
    affected_nodes.insert(deleted_node);

    // 阶段3：批量重分配
    for (int u : affected_nodes) {
        int new_assign = find_closest_ring_node(s, u);

        // 有效性验证（防止无效分配）
        if (new_assign == -1 || !s.in_ring[new_assign]) {
            cerr << "Invalid assignment for node " << u << endl;
            continue;
        }
        // 更新分配成本（差分计算优化性能）
        if (s.assignments.count(u)) {
            s.assign_cost -= (10 - alpha) * graph.assign_cost[u][s.assignments[u]];
        }
        s.assign_cost += (10 - alpha) * graph.assign_cost[u][new_assign];
        // 更新分配关系
        s.assignments[u] = new_assign;
    }
}

Solution AMNSSolver::Random_drop_one2(Solution s) const {
    if (s.ring.size() <= 3) return s;

    // 阶段1：删除节点
    int pos = 1 + rand() % (s.ring.size() - 1);
    int deleted_node = s.ring[pos];
    if (deleted_node == 0) return s;
    s.ring.erase(s.ring.begin() + pos);
    s.in_ring[deleted_node] = false;

    // 阶段2：构建受影响节点集合
    unordered_set<int> affected_nodes;

    // 收集原分配至删除节点的非环节点
    for (auto& entry : s.assignments) {//非环点 u 环点 v
        int u = entry.first;
        int v = entry.second;
        if (v == deleted_node && !s.in_ring[u]) {
            affected_nodes.insert(u);
        }
    }

    // 添加被删节点自身（需重新分配）
    affected_nodes.insert(deleted_node);

    // 阶段3：批量重分配
    for (int u : affected_nodes) {
        int new_assign = find_closest_ring_node(s, u);
        // 有效性验证（防止无效分配）
        if (new_assign == -1 || !s.in_ring[new_assign]) {
            cerr << "Invalid assignment for node " << u << endl;
            continue;
        }
        // 更新分配成本（差分计算优化性能）
        if (s.assignments.count(u)) {
            s.assign_cost -= (10 - alpha) * graph.assign_cost[u][s.assignments[u]];
        }
        s.assign_cost += (10 - alpha) * graph.assign_cost[u][new_assign];
        // 更新分配关系
        s.assignments[u] = new_assign;
    }
    return s;
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

void AMNSSolver::randomized_two_opt(Solution& s) const {

    if (s.ring.size() < 4) return;
    int i = rand() % (s.ring.size() - 1);
    int j = rand() % (s.ring.size() - 1);
    // 随机选择两个不同的边
    while (abs(i - j) < 2) {
        i = rand() % (s.ring.size() - 1);
        j = rand() % (s.ring.size() - 1);
    }
    // 确保i < j
    if (i > j) swap(i, j);
    //j++; // 调整为j+1的索引

    // 计算当前成本
    double original = graph.routing_cost[s.ring[i]][s.ring[i + 1]]
        + graph.routing_cost[s.ring[j]][s.ring[(j + 1) % s.ring.size()]];

    // 计算交换后成本
    double modified = graph.routing_cost[s.ring[i]][s.ring[j]]
        + graph.routing_cost[s.ring[i + 1]][s.ring[(j + 1) % s.ring.size()]];

    // 如果成本降低
    if (modified < original) {
        reverse(s.ring.begin() + i + 1, s.ring.begin() + j + 1);
    }
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

int AMNSSolver::random_unused_node(const Solution& s) const {
    vector<int> candidates;
    for (int i = 0; i < graph.nodes.size(); ++i) {
        if (!s.in_ring[i] && i != graph.depot) {
            candidates.push_back(i);
        }
    }
    if (candidates.empty()) return -1;
    return candidates[rand() % candidates.size()];
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