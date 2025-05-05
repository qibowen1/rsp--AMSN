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
    for (int run = 0; run < MAX_ITER; ++run) {
        Solution temp_best;
        Solution current;
        int solutionPoolSize = 5; // 初始解池大小
        int GVNS_SIZE = 1;
        current = greedy_construction(solutionPoolSize); // 生成贪婪初始解 Random_construction greedy_construction
        temp_best = current;
        int search_iter = 0;
        cout << "初始解成本: " << current.total_cost() << endl;
        Solution neighbor = variable_neighborhood_VNS(current); // 生成邻域解 variable_neighborhood_VNS
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

Solution AMNSSolver::GVNS(Solution s, int maxIter) const {//基本无提升？
    Solution temp_best = s;
    for (int i = 0;i < maxIter;i++) {
        int k = rand() % 3; // 随机选择shake 长度
        Solution after_shake = shaking(s, k); // 执行扰动
        Solution after_ls = local_search_VND(after_shake, true); // 执行局部搜索
        if (after_ls.total_cost() < temp_best.total_cost()) {
            temp_best = after_ls; // 更新当前最优解
            s = after_ls;
        }
    }
    return temp_best;
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



Solution AMNSSolver::Random_construction() const {
    Solution s;
    s.initialize(graph.nodes.size());
    vector<int> candidates;
    for (int i = 0; i < graph.nodes.size(); ++i)
        if (i != graph.depot) candidates.push_back(i);
    s.in_ring[graph.depot] = true;
    // 随机选择节点构建环（至少3个节点）
    shuffle(candidates.begin(), candidates.end(), rng);
    int ring_size = max(3, (int)(candidates.size() * 0.8)); // 随机环大小（50%节点）
    std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
    // 构建初始环
    for (int i = 0; i < ring_size - 1; ++i) { // depot已包含
        if (candidates.empty()) break;

        // 直接取打乱后的前N个节点（避免重复）
        int v = candidates.back();
        s.ring.push_back(v);
        s.in_ring[v] = true;
        candidates.pop_back(); // 确保节点只选一次
    }
    //s.ring.push_back(graph.depot); // 闭合环

    // 分配非环节点
    for (int u = 0; u < graph.nodes.size(); ++u) {
        if (!s.in_ring[u]) {
            int closest = find_closest_ring_node(s, u);
            s.assignments[u] = closest;
        }
    }
    evaluate(s);
    return s;
}

Solution AMNSSolver::variable_neighborhood_search(Solution s) const {//1
    Solution best_sol = s;
    Solution improved = s;
    Solution temp = s;
    Solution current = s;
    int stage = 0;//停滞次数
    int max_stage = 3;
    int tempbest = s.total_cost();
    do {
        tempbest = current.total_cost();
        temp = local_search_VND(current,true);//local_search_VND local_search_v2
        if (temp.total_cost() <= tempbest) {
            if (temp.total_cost() == tempbest) {
                stage++;
                if (stage >= max_stage) {
                    break; // 达到最大停滞次数，退出
                }
            }
            else {
                stage = 0; // 重置停滞计数器
            }
            tempbest = temp.total_cost();
            improved = temp;
            current = temp;
        }
        else {
            stage++;
            break;
        }
    } while (true);
    if (best_sol.total_cost() > tempbest) {
        best_sol = improved;
    }
    return best_sol;
}

Solution AMNSSolver::variable_neighborhood_VNS(Solution s) const {//2
    const int k_max = 5;  // 最大邻域层级
    int k = 0;            // 当前扰动强度
    Solution best_sol = s;
    Solution current = s;
    cout << "VNS初始解成本: " << current.total_cost() << endl;
    do {
        // Shaking阶段：根据当前k值生成扰动解
        Solution shaken_sol = shaking(current, k);
        evaluate(shaken_sol); // 评估扰动解
        // 使用VND进行局部搜索（使用您现有的local_search_v2）
        Solution local_optima = local_search_VND(shaken_sol,true);

        // 接受准则
        if (local_optima.total_cost() < best_sol.total_cost()) {
            best_sol = local_optima;
            current = local_optima;
            k = 1;  // 重置扰动强度
        }
        else {
            k++;    // 增强扰动强度
        }
    } while (k <= k_max);
    cout << "VNS后的成本: " << current.total_cost() << endl;
    return best_sol;
}


//方法3 随机变量邻域搜索，使用VND进行局部搜索，直到没有改进为止
Solution AMNSSolver::local_search_VND(Solution s, bool fastmodel) const {
    int k = 0;
    int max_iter = 4;
    int itao = 2;
    int localsearchMax = graph.nodes.size() / itao;
    do {
        Solution local_best = s;
        int choose = rand() % max_iter;
        for (int i = 0;i < localsearchMax;i++) {
            vector<Solution> H;
            Solution temp = s;
            // 生成当前类型的邻域
            switch (choose) {
            case 0: temp = Random_drop_one2(s);
            case 1: temp = Random_add_one2(s);
            case 2: H = generate_random_add_drop_neighbors(s); break;
            case 3:temp = generate_random_2opt_neighbors(s); break;
            default:break;
            }
            evaluate(temp);
            if (temp.total_cost() < local_best.total_cost()) {
                // 计算当前邻域的解
                local_best = temp;
            }
            if (!H.empty()) {
                // 找到当前邻域中的最优解
                auto best_it = min_element(H.begin(), H.end(),
                    [](const auto& a, const auto& b) {
                        return a.total_cost() < b.total_cost();
                    });
                int temp_res = best_it->total_cost();
                if (temp_res < local_best.total_cost()) {
                    local_best = *best_it;
                }
            }
        }
        // 如果找到改进解
        if (local_best.total_cost() < s.total_cost()) {
            k = 0;
            s = local_best;
        }
        else k++;
    } while (k < max_iter);  // 只要任一邻域有改进就继续

    return s;
}



//random
vector<Solution> AMNSSolver::generate_random_add_drop_neighbors(Solution s) const {
    vector<Solution> neighbors;
    Solution neighbor = s;
    int temp = neighbor.total_cost();
    Random_drop_one(neighbor);
    evaluate(neighbor);
    if (neighbor.total_cost() < temp) {
        neighbors.push_back(neighbor);
    }
    Random_add_one(neighbor);
    evaluate(neighbor);
    if (neighbor.total_cost() < temp) {
        neighbors.push_back(neighbor);
    }

    return neighbors;
}


Solution AMNSSolver::generate_random_2opt_neighbors(Solution s) const {
    randomized_two_opt(s);
    evaluate(s);
    return s;
}

Solution AMNSSolver::shaking(Solution s, int k) const {
    Solution shaken = s;
    const int depot = graph.depot;
    const int min_ring_size = 3;
    // 根据扰动强度k决定操作次数
    int perturbation_steps = k + 1;
    for (int step = 0; step < perturbation_steps; ++step) {
        // 随机选择扰动类型
        int op_type = rand() % 4;
        switch (op_type) {
        case 0: // 随机交换环内两个节点
            if (shaken.ring.size() > min_ring_size) {
                int i = 1 + rand() % (shaken.ring.size() - 1);
                int j = 1 + rand() % (shaken.ring.size() - 1);
                std::swap(shaken.ring[i], shaken.ring[j]);
            }
            break;
        case 2: // 随机反转环的子段
            if (shaken.ring.size() > min_ring_size) {
                int i = 1 + rand() % (shaken.ring.size() - 2);
                int j = i + 1 + rand() % (shaken.ring.size() - i - 1);
                std::reverse(shaken.ring.begin() + i, shaken.ring.begin() + j);
            }
            break;
        default: break; // 无效操作
        }
    }
    // 重新计算成本
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