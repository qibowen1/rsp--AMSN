#include "solver.h"
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
结合动态邻域选择、以平衡全局探索与局部 开发能力。
*/
/*
a: 3 5 7 9
max_search_iter:对每个初始解的迭代次数 3 ls
RCL_ratio:rcl比例 0.8
neighborhood_types:邻域类型数量 5 4
*/
AMNSSolver::AMNSSolver(const RSPGraph& g, int a, int max_search_iter, double RCL_ratio, int neighborhood_types) 
    : graph(g), alpha(a),max_search_iter(max_search_iter),RCL_ratio(RCL_ratio),neighborhood_types(neighborhood_types) {}

void AMNSSolver::solve() { 
	Solution current;
	int solutionPoolSize = 10; // 初始解池大小
	current = greedy_construction(solutionPoolSize); // 生成贪婪初始解 Random_construction greedy_construction
    best = current;
    int search_iter = 0;
    cout << "初始解成本: " << current.total_cost() << endl;
    Solution neighbor = variable_neighborhood_search(current); // 生成邻域解
    evaluate(neighbor);              // 重新计算新解的成本
    double delta = neighbor.total_cost() - current.total_cost();
    if (delta < 0) {
        current = neighbor;
        best = current;
    }
    Solution after_GVNS=GVNS(neighbor, 2); // 执行GVNS
	if (after_GVNS.total_cost() < best.total_cost()) {
		best = after_GVNS;
	}
    cout << "in迭代 #" << " 当前最优解成本: " << best.total_cost()
        << " (路由: " << best.routing_cost
        << ", 分配: " << best.assign_cost << ")" << endl;
}

Solution AMNSSolver::GVNS(Solution s, int maxIter) const {//基本无提升？
    Solution temp_best = s;
    for (int i = 0;i < maxIter;i++) {
        int k = rand() % 3; // 随机选择shake 长度
        Solution after_shake = shaking(s, k); // 执行扰动
        evaluate(after_shake); // 评估扰动后的解
        Solution after_ls = local_search_VND(after_shake); // 执行局部搜索
        evaluate(after_ls); // 评估局部搜索后的解
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
        for (int i = 0;i < solutionPoolSize;i++) {
            Solution s;
            s.initialize(graph.nodes.size());     // 初始化解决方案结构
            vector<int> candidates;               // 候选节点集合（非根节点）
            for (int i = 0; i < graph.nodes.size(); ++i)
                if (i != graph.depot) candidates.push_back(i);
            double prev_cost = INF;  // 记录前一次的总成本
            int stagnation = 0;      // 成本未改善的连续次数            

            while (!candidates.empty()) {         // 逐步构建环路径
                vector<pair<int, double>> candidate_scores; // 候选节点评分集合
                // 计算每个候选节点的插入评分
                for (int v : candidates) {
                    double min_ring = INF, min_assign = INF;
                    // 1. 计算插入环的最小路由成本变化
                    for (int i = 0; i < s.ring.size(); ++i) {
                        int next = (i + 1) % s.ring.size();
                        double delta = graph.routing_cost[s.ring[i]][v]
                            + graph.routing_cost[v][s.ring[next]]
                                - graph.routing_cost[s.ring[i]][s.ring[next]];
                            min_ring = min(min_ring, delta);
                    }

                    // 2. 计算该节点能带来的分配成本收益
                    double assign_gain = 0;
					int avg_assign = 0;
                    for (int u = 0; u < graph.nodes.size(); ++u) {
                        if (!s.in_ring[u]) {      // 只考虑当前未分配的节点
                            /*double new_cost = graph.assign_cost[u][v];
                            double old_cost = graph.assign_cost[u][s.assignments[u]];
                            assign_gain += new_cost - old_cost;*/
							avg_assign += graph.assign_cost[u][v];
                        }
                    }
					avg_assign = avg_assign / (graph.nodes.size() - s.ring.size()); // 平均分配收益
                    // 综合评分（路由成本变化+分配收益）
                    candidate_scores.emplace_back(v, alpha * min_ring + (10 - alpha) * avg_assign);//
                }
                // 构建限制候选列表（RCL）
                sort(candidate_scores.begin(), candidate_scores.end(),
                    [](auto& a, auto& b) { return a.second < b.second; });

                int rcl_size = (int)(candidate_scores.size() * RCL_ratio); // RCL大小取前80%
                if (rcl_size == 0) break;
                int selected = candidate_scores[rand() % rcl_size].first; // 随机选择RCL中的节点
                insert_node(s, selected);         // 执行节点插入
                candidates.erase(remove(candidates.begin(), candidates.end(), selected), candidates.end());
                
                evaluate(s);  // 重新计算总成本
                /* 新增终止条件判断 */
                if (s.ring.size() >= 3) { // 环长度达标后开始检查成本变化
                    if (s.total_cost() >= prev_cost) {
                        break; // 成本未改善达到阈值则终止
                    }
                    else {
                        prev_cost = s.total_cost(); // 更新基准成本
                    }
                }
                else {
                    prev_cost = s.total_cost(); // 环长度不足时继续构建
                }
                
        }
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
    int max_stage = 2;
	int tempbest = s.total_cost();
        do {
            tempbest = current.total_cost();
            temp = local_search_VND(current);//local_search_VND local_search_v2
            if (temp.total_cost() <= tempbest) {
                if (temp.total_cost() == tempbest) {
					stage++;
					if (stage >= max_stage) {
						break; // 达到最大停滞次数，退出
					}
				}else {
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
		cout << "扰动解成本: " << shaken_sol.total_cost() << endl;
        // 使用VND进行局部搜索（使用您现有的local_search_v2）
        Solution local_optima = local_search_VND(shaken_sol);

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

Solution AMNSSolver::local_search_v2(Solution s) const {//3
    Solution current = s;
    do {
        // 增量生成邻域
        vector<Solution> H;
        generate_promising_neighbors(current, H);

        auto best_it = min_element(H.begin(), H.end(),
            [](const auto& a, const auto& b) { return a.total_cost() < b.total_cost(); });

        if (H.size() != 0 && best_it->total_cost() < current.total_cost()) {
            current = *best_it;
        }
        else {
            break; // 无改进则终止
        }
    } while (true);
    return current;
}

Solution AMNSSolver::local_search_VND(Solution s) const {
    int k = 0;
    int max_iter = 4;
    int localsearchMax = graph.nodes.size()/5;
    do {
        Solution local_best = s;
		int choose = rand() % max_iter;
        for(int i=0;i<localsearchMax;i++){
            vector<Solution> H;
			Solution temp = s;
            // 生成当前类型的邻域
            switch (choose) {
            case 0: temp = Random_drop_one2(s);
            case 1: temp = Random_add_one2(s);
            case 2: H = generate_random_add_drop_neighbors(s); break;
			case 3:temp=generate_random_2opt_neighbors(s); break;
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
        if (local_best.total_cost()<s.total_cost()) {
             k = 0;
             s = local_best;
         }else k++;
    } while (k<max_iter);  // 只要任一邻域有改进就继续

    return s;
}

// 其他generate_xxx_neighbors函数类似修改

void AMNSSolver::generate_promising_neighbors(Solution s, vector<Solution>& H) const {
    // 1. 添加邻域
    auto add_neighbors = generate_add_neighbors(s);
    H.insert(H.end(), add_neighbors.begin(), add_neighbors.end());
    // 2. 删除邻域
    auto drop_neighbors = generate_Dinf_drop_neighbors(s);
    H.insert(H.end(), drop_neighbors.begin(), drop_neighbors.end());
    // 3. 添加+删除组合邻域
    //auto add_drop_neighbors = generate_add_drop_neighbors(s);
    //H.insert(H.end(), add_drop_neighbors.begin(), add_drop_neighbors.end());
    // 4. 优化邻域
    auto opt_neighbors = generate_opt_neighbors(s);
    H.insert(H.end(), opt_neighbors.begin(), opt_neighbors.end());

}

vector<Solution> AMNSSolver::generate_add_neighbors(Solution s) const {
    vector<Solution> neighbors;
    vector<int> candidates;
    
    // 收集所有未使用的节点
    for (int i = 0; i < graph.nodes.size(); ++i) {
        if (!s.in_ring[i] && i != graph.depot) {
            candidates.push_back(i);
        }
    }
    if (candidates.empty()) return neighbors;
    for (int v:candidates) {
        Solution neighbor = s;
        insert_node(neighbor, v);//每次都插入吗？
        evaluate(neighbor);
        if (neighbor.total_cost() < s.total_cost()) {
            neighbors.push_back(neighbor);
        }   
    }
    return neighbors;
}


vector<Solution> AMNSSolver::generate_Dinf_drop_neighbors(Solution s) const {
    vector<Solution> neighbors;

    // 排除根节点和最小环长约束
    if (s.ring.size() <= 3) return neighbors;
    Solution neighbor = s;

    neighbor= Dinf_drop_op(neighbor);
    if (neighbor.total_cost() < s.total_cost())
        neighbors.push_back(neighbor);
    return neighbors;
}

vector<Solution> AMNSSolver::generate_drop_neighbors2(Solution s) const {
    vector<Solution> neighbors;

    // 排除根节点和最小环长约束
    if (s.ring.size() <= 3) return neighbors;

    // 尝试删除每个非根节点
    for (int i = 1; i < s.ring.size(); ++i) {
        Solution neighbor = s;
        int deleted_node = neighbor.ring[i];
        drop_redistribute(neighbor,deleted_node);
        evaluate(neighbor);
        if (neighbor.total_cost() < s.total_cost())
            neighbors.push_back(neighbor);
    }
    return neighbors;
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

vector<Solution> AMNSSolver::generate_opt_neighbors(Solution s) const {
    vector<Solution> neighbors;

    // 1. 增强版优化（1个高质量解）
    Solution enhanced_neighbor = s;
    enhanced_full_two_opt(enhanced_neighbor); //enhanced_full_two_opt  randomized_two_opt
    evaluate(enhanced_neighbor);
    if (enhanced_neighbor.total_cost() < s.total_cost())
        neighbors.push_back(enhanced_neighbor);
    return neighbors;
}

Solution AMNSSolver::generate_random_2opt_neighbors(Solution s) const {
	randomized_two_opt(s);
	evaluate(s);
	return s;
}

//Solution AMNSSolver::shaking(Solution s, int k) const {
//    Solution temp = s;
//    const int depot = graph.depot; // 根节点索引
//    if (k == 0) return temp;   //不扰动
//
//    for (int i = 0; i < k; ++i) {
//        // 随机选择两个不同节点
//        int vi = rand() % (graph.nodes.size() - 1) + 1; // 生成1~n-1的随机数
//        int vj = rand() % (graph.nodes.size() - 1) + 1;
//        while (vj == vi) vj = rand() % (graph.nodes.size() - 1) + 1;
//        bool vi_in = temp.in_ring[vi];
//        bool vj_in = temp.in_ring[vj];
//
//        // 执行对应操作
//        if (vi_in && !vj_in) {      // Case 1: vi在环，vj不在
//            Random_drop_one(temp);
//            Random_add_one(temp);
//        }
//        else if (vi_in && vj_in) {  // Case 2: 两者都在环
//            Random_drop_one(temp);
//        }
//        else if (!vi_in && !vj_in) { // Case 3: 两者都不在环
//            Random_add_one(temp);
//        }
//    }
//    return temp;
//}

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

        //case 1: // 随机删除-插入环外节点
        //    if (shaken.ring.size() > min_ring_size) {
        //        int pos = 1 + rand() % (shaken.ring.size() - 1);
        //        int deleted = shaken.ring[pos];
        //        shaken.ring.erase(shaken.ring.begin() + pos);
        //        shaken.in_ring[deleted] = false;

        //        int v = random_unused_node(shaken);
        //        if (v != -1) {
        //            insert_node(shaken, v);
        //        }
        //    }
        //    break;

        case 2: // 随机反转环的子段
            if (shaken.ring.size() > min_ring_size) {
                int i = 1 + rand() % (shaken.ring.size() - 2);
                int j = i + 1 + rand() % (shaken.ring.size() - i - 1);
                std::reverse(shaken.ring.begin() + i, shaken.ring.begin() + j);
            }
            break;

        //case 3: // 破坏性扰动：删除连续k+1个节点
        //    if (shaken.ring.size() > min_ring_size + k) {
        //        int start = 1 + rand() % (shaken.ring.size() - k - 1);
        //        shaken.ring.erase(shaken.ring.begin() + start, shaken.ring.begin() + start + k + 1);
        //        // 随机添加新节点补足环
        //        for (int i = 0; i <= k; ++i) {
        //            int v = random_unused_node(shaken);
        //            if (v != -1) insert_node(shaken, v);
        //        }
        //    }
        //    break;
		default: break; // 无效操作
        }
	}
	// 重新计算成本
	evaluate(shaken);
	return shaken;
}

void AMNSSolver::drop_redistribute(Solution& s, int pos) const {
    // 阶段2：构建受影响节点集合
    unordered_set<int> affected_nodes;
	int deleted_node = s.ring[pos];
    s.ring.erase(s.ring.begin() + pos);
    s.in_ring[deleted_node] = false;
    // 收集原分配至删除节点的非环节点
    for (auto& entry : s.assignments) {//非环点 u 环点 v
        int u = entry.first;
        int v = entry.second;
        if (v == deleted_node && !s.in_ring[u]) {
            affected_nodes.insert(u);
        }
    }
    affected_nodes.insert(deleted_node); // 被删节点自身

    // 重新分配
    for (int u : affected_nodes) {
        if (s.in_ring[u]) continue;

        int new_assign = find_closest_ring_node(s, u);
        if (new_assign == -1) continue;

        // 更新分配关系
        s.assign_cost += (10 - alpha) * (graph.assign_cost[u][new_assign] -
            graph.assign_cost[u][s.assignments[u]]);
        s.assignments[u] = new_assign;
    }
}



//vector<Solution> AMNSSolver::generate_3opt_neighbors(Solution s) const {
//    vector<Solution> neighbors;
//
//    // 1. 增强版优化（1个高质量解）
//    Solution enhanced_neighbor = s;
//    randomized_three_opt(enhanced_neighbor);
//    evaluate(enhanced_neighbor);
//    if (enhanced_neighbor.total_cost() < s.total_cost())
//        neighbors.push_back(enhanced_neighbor);
//
//    return neighbors;
//}

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
Solution AMNSSolver::Dinf_drop_op(Solution s) const {
	Solution best_res = s;
    int best_save = 0;
	for (int i = 1; i < s.ring.size(); ++i) {//0节点除外
		Solution temp = s;
        int cur_best_save = 0;
        int pre = i - 1;
		int post = (i + 1) % temp.ring.size();
		int deleted_node = temp.ring[i];
		int pre_node = temp.ring[pre];
		int post_node = temp.ring[post];
		int saving = alpha*(graph.routing_cost[pre_node][deleted_node] + graph.routing_cost[deleted_node][post_node] -
            graph.routing_cost[pre_node][post_node]);//saving值

        temp.ring.erase(temp.ring.begin() + i);
        temp.in_ring[deleted_node] = false;
        // 构建受影响节点集合
        unordered_set<int> affected_nodes;
        // 收集原分配至删除节点的非环节点
        for (auto& entry : temp.assignments) {//非环点 u 环点 v
            int u = entry.first;
            int v = entry.second;
            if (v == deleted_node) {
                affected_nodes.insert(u);
            }
        }
        affected_nodes.insert(deleted_node);

        // 添加被删节点自身（需重新分配）
        affected_nodes.insert(deleted_node);
        int penalty = 0;
        // 阶段3：批量重分配
        for (int u : affected_nodes) {
            int new_assign = find_closest_ring_node(temp, u);

            // 有效性验证（防止无效分配）
            if (new_assign == -1 || !temp.in_ring[new_assign]) {
                cerr << "Invalid assignment for node " << u << endl;
                continue;
            }
			penalty += (10 - alpha) * (graph.assign_cost[u][new_assign] - graph.assign_cost[u][temp.assignments[u]]);
            // 更新分配关系
            temp.assignments[u] = new_assign;
        }
		cur_best_save = saving - penalty;
		evaluate(temp); // 重新计算总成本
        if (cur_best_save > 0 &&cur_best_save>best_save) {
			best_res = temp; // 更新最优解
			best_save = cur_best_save; // 更新最优节省值
        }
	}
    return best_res;

}

void AMNSSolver::Random_drop_one(Solution& s) const {
    if (s.ring.size() <= 3) return;

    // 阶段1：删除节点
    int pos = 1 + rand() % (s.ring.size() - 1);
    int deleted_node = s.ring[pos];
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
    s.ring.erase(s.ring.begin() + pos);
    s.in_ring[deleted_node] = false;

    // 阶段2：构建受影响节点集合
    unordered_set<int> affected_nodes;

    // 收集原分配至删除节点的非环节点
    for (auto& entry : s.assignments) {//非环点 u 环点 v
        int u = entry.first;
        int v = entry.second;
        if (v == deleted_node&& !s.in_ring[u]) {
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

void AMNSSolver::enhanced_two_opt(Solution& s) const {
    double best_gain = 0;
    int best_i = -1, best_j = -1;

    // 遍历所有可能的边交换
    for (int i = 0; i < s.ring.size() - 1; ++i) {
        for (int j = i + 2; j < s.ring.size(); ++j) {
            double original = graph.routing_cost[s.ring[i]][s.ring[i + 1]]
                + graph.routing_cost[s.ring[j]][s.ring[(j + 1) % s.ring.size()]];
            double modified = graph.routing_cost[s.ring[i]][s.ring[j]]
                + graph.routing_cost[s.ring[i + 1]][s.ring[(j + 1) % s.ring.size()]];
            double gain = original - modified;

            if (gain > best_gain) {
                best_gain = gain;
                best_i = i;
                best_j = j;
            }
        }
    }

    // 执行最佳交换
    if (best_gain > 0) {
        reverse(s.ring.begin() + best_i + 1, s.ring.begin() + best_j + 1);
    }
}

void AMNSSolver::enhanced_full_two_opt(Solution& s) const {
    bool improved;
    do {
        improved = false;
        double best_gain = 0;
        int best_i = -1, best_j = -1;

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
                if (gain > 0) {
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



//void AMNSSolver::randomized_three_opt(Solution& s) const {
//    const int MAX_ATTEMPTS = 30; // 最大无改进尝试次数
//    int attempts = 0;
//    bool improved;
//    double original_cost = s.routing_cost;
//
//    do {
//        improved = false;
//
//        // 随机选择三个不同的边
//        int i = rand() % (s.ring.size() - 5);
//        int j = i + 1 + rand() % (s.ring.size() - i - 4);
//        int k = j + 1 + rand() % (s.ring.size() - j - 3);
//
//        // 获取当前路径段的实际节点索引
//        int a = s.ring[i];
//        int b = s.ring[(i + 1) % s.ring.size()];
//        int c = s.ring[j];
//        int d = s.ring[(j + 1) % s.ring.size()];
//        int e = s.ring[k];
//        int f = s.ring[(k + 1) % s.ring.size()];
//
//        // 计算原始段的成本
//        double original_segment = alpha * (
//            graph.routing_cost[a][b] +
//            graph.routing_cost[c][d] +
//            graph.routing_cost[e][f]
//            );
//
//        // 七种可能的连接方式（实际有效四种）
//        vector<pair<double, vector<int>>> candidates;
//
//        // 方式1：反转中间两段
//        vector<int> path1 = s.ring;
//        reverse(path1.begin() + i + 1, path1.begin() + j + 1);
//        reverse(path1.begin() + j + 1, path1.begin() + k + 1);
//        double cost1 = calculate_segment_cost(path1, i, k);
//
//        // 方式2：反转第三段
//        vector<int> path2 = s.ring;
//        reverse(path2.begin() + j + 1, path2.begin() + k + 1);
//        double cost2 = calculate_segment_cost(path2, i, k);
//
//        // 方式3：反转前两段
//        vector<int> path3 = s.ring;
//        reverse(path3.begin() + i + 1, path3.begin() + j + 1);
//        double cost3 = calculate_segment_cost(path3, i, k);
//
//        // 方式4：反转全部三段
//        vector<int> path4 = s.ring;
//        reverse(path4.begin() + i + 1, path4.begin() + k + 1);
//        double cost4 = calculate_segment_cost(path4, i, k);
//
//        // 找到最佳改进
//        double best_cost = min({ cost1, cost2, cost3, cost4 });
//        double gain = original_segment - best_cost;
//
//        if (gain > 0) {
//            // 应用最佳改进
//            if (best_cost == cost1) s.ring = path1;
//            else if (best_cost == cost2) s.ring = path2;
//            else if (best_cost == cost3) s.ring = path3;
//            else s.ring = path4;
//
//            // 更新路由成本（差分更新）
//            s.routing_cost += (best_cost - original_segment) / alpha;
//            improved = true;
//            attempts = 0; // 重置计数器
//        }
//        else {
//            attempts++;
//        }
//
//    } while (improved || attempts < MAX_ATTEMPTS);
//}

//double AMNSSolver::calculate_segment_cost(const vector<int>& path, int start, int end) const {
//    double cost = 0;
//    for (int i = start; i <= end; ++i) {
//        int j = (i + 1) % path.size();
//        cost += graph.routing_cost[path[i]][path[j]];
//    }
//    return alpha * cost;
//}


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