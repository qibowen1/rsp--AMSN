#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <numeric>
#include <limits>
#include <unordered_set>  // 必须包含该头文件
#include <fstream>
#include <sstream>
#include <string>
#include <random>  // 随机数库
using namespace std;

const double INF = numeric_limits<double>::max(); // 定义无穷大常量

// 在需要使用shuffle的函数内部定义
std::random_device rd;  // 随机种子
std::mt19937 rng(rd()); // 使用Mersenne Twister引擎
/*
针对Ring Star Problem（RSP），​​自适应多阶段邻域搜索算法（Adaptive Multi-phase Neighborhood Search, AMNS）​​，
结合动态邻域选择、模拟退火机制和分配优化策略，以平衡全局探索与局部 开发能力。
*/
// 图结构定义：包含节点坐标、路由成本和分配成本矩阵
struct RSPGraph {
    vector<pair<double, double>> nodes;      // 节点坐标集合
    vector<vector<double>> routing_cost;     // 路由成本矩阵（环边权重）
    vector<vector<double>> assign_cost;       // 分配成本矩阵（星型边权重）
    int depot = 0;                           // 根节点索引（默认为0）

    void initialize(int size) {              // 初始化矩阵大小
        routing_cost.resize(size, vector<double>(size));
        assign_cost.resize(size, vector<double>(size));
    }
};

// 解决方案结构定义
struct Solution {
    vector<int> ring;                        // 环路径节点索引集合
    unordered_map<int, int> assignments;     // 分配关系：非环节点->最近环节点
    vector<bool> in_ring;                     // 标记节点是否在环中
    double routing_cost = 0;                 // 环路径总成本
    double assign_cost = 0;                  // 分配总成本

    double total_cost() const {               // 计算解决方案总成本
        return routing_cost + assign_cost;
    }

    void initialize(int size) {               // 初始化解决方案
        in_ring.resize(size, false);
        ring.push_back(0);                    // 初始时环只包含根节点
        in_ring[0] = true;
    }
};

// 自适应多阶段邻域搜索算法类
class AMNSSolver {
private:
    RSPGraph graph;                          // 问题图实例
	int alpha ;                        // 分配权重系数（可调节）
    Solution best;                            // 历史最优解

    // 算法参数配置
    const int NUM_OPS = 6;                   // 邻域操作类型总数
    // 新增VNS相关参数
    const int VNS_MAX_LEVEL = 3;        // 最大邻域层级
    const int VNS_LOCAL_SEARCH_STEPS = 50; // 本地搜索步数

public:
    AMNSSolver(const RSPGraph& g, int a) : graph(g), alpha(a) { // 构造函数初始化

    }

    // 主求解函数（max_iter: 最大迭代次数）
    void solve(int max_iter) {
        Solution current;
        best = current;   
        current = greedy_construction(); // 生成贪婪初始解
        double curr_min = current.total_cost();
        double temp = 1000.0; // 初始温度
        const double cooling_rate = 0.995;
        for (int iter = 0; iter < max_iter; ++iter) {
			Solution neighbor = variable_neighborhood_search(current); // 生成邻域解
            evaluate(neighbor);              // 重新计算新解的成本
			//cout << "迭代 #" << iter << " 当前成本: " << neighbor.total_cost() << endl;
			//cout << "neighbor.total_cost():" << neighbor.total_cost() << "current.total_cost():" << current.total_cost() <<"当前最优："<<curr_min<<endl;
            // 模拟退火接受准则
            double delta = neighbor.total_cost() - current.total_cost();
            if (delta < 0 || exp(-delta / temp) >(rand() % 1000) / 1000.0) {
                current = neighbor;
            }

            // 更新温度
            temp *= cooling_rate;
            // 自适应退火调度
            

            // 定期重置为历史最优解
            if (iter!=0&&iter % 100 == 0) {
                current = best;
            }
            if (current.total_cost() < curr_min) { // 更新全局最优解
                curr_min = current.total_cost();
                best = current;
            }

        }
		cout << "in迭代 #" << max_iter << " 当前最优解成本: " << curr_min
			<< " (路由: " << best.routing_cost
			<< ", 分配: " << best.assign_cost << ")" << endl;
    }

    Solution getBestSolution() {
        return best; // 返回当前最优解 
    }

private:
    // 贪婪随机化构造初始解（GRASP策略）
    Solution greedy_construction() {
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
                for (int u = 0; u < graph.nodes.size(); ++u) {
                    if (!s.in_ring[u]) {      // 只考虑当前未分配的节点
                        double new_cost = graph.assign_cost[u][v];
                        double old_cost = graph.assign_cost[u][s.assignments[u]];
                        assign_gain += new_cost-old_cost;
                    }
                }
				assign_gain = assign_gain / (graph.nodes.size() - s.ring.size()); // 平均分配收益
                // 综合评分（路由成本变化+分配收益）
                candidate_scores.emplace_back(v, alpha * min_ring + (10-alpha) * assign_gain);
            }

            // 构建限制候选列表（RCL）
            sort(candidate_scores.begin(), candidate_scores.end(),
                [](auto& a, auto& b) { return a.second < b.second; });

            int rcl_size = (int)(candidate_scores.size() * 0.2); // RCL大小取前80%
            if (rcl_size == 0) break;
            int selected = candidate_scores[rand() % rcl_size].first; // 随机选择RCL中的节点
			//cout << "选择节点: " << selected << endl;
            insert_node(s, selected);         // 执行节点插入
            candidates.erase(remove(candidates.begin(), candidates.end(), selected), candidates.end());
            evaluate(s);  // 重新计算总成本
            /* 新增终止条件判断 */
            if (s.ring.size() >= 3) { // 环长度达标后开始检查成本变化
                if (s.total_cost() >= prev_cost) {
                    break; // 成本未改善达到阈值则终止
                }
                else {
                    stagnation = 0;       // 重置未改善计数器
                    prev_cost = s.total_cost(); // 更新基准成本
                }
            }
            else {
                prev_cost = s.total_cost(); // 环长度不足时继续构建
            }
        }

        /* 确保环的最小长度约束 */
        while (s.ring.size() < 3 && !candidates.empty()) {
            int v = candidates.back();
            insert_node(s, v);
            candidates.pop_back();
            evaluate(s);
        }

        return s;
    }
    Solution variable_neighborhood_search(Solution s) {
        Solution best_sol=s;
        double best_cost = s.total_cost();
        Solution current = s;         // 维护当前解
        for (int k = 1; k <= VNS_MAX_LEVEL; ++k) {
            Solution shaken = shaking(current, k); // 扰动阶段
            Solution improved = local_search_v2(shaken); // 本地搜索
            if (improved.total_cost() < current.total_cost()) {
                current = improved;  // 更新当前解为改进解
                best_sol = improved.total_cost() < best_cost ? improved : best_sol;
                k = 1; // 重置邻域层级
            }
        }
        return best_sol;
    }

    /*
    Solution shaking(Solution s, int k) {
    // Case 1: vi在环且vj不在 → 替换操作
    // Case 2: 两者都在环 → 删除vi
    // Case 3: 两者都不在 → 添加vj

    */
    Solution shaking(Solution s, int k) {
        Solution temp = s;
        const int depot = graph.depot; // 根节点索引

        for (int i = 0; i < k; ++i) {
            // 随机选择两个不同节点
            int vi = rand() % (graph.nodes.size() - 1) + 1; // 生成1~n-1的随机数
            int vj = rand() % (graph.nodes.size() - 1) + 1;
            while (vj == vi) vj = rand() % (graph.nodes.size() - 1) + 1;

            bool vi_in = temp.in_ring[vi];
            bool vj_in = temp.in_ring[vj];

            // 执行对应操作
            if (vi_in && !vj_in) {      // Case 1: vi在环，vj不在
                // 删除vi（确保不是根节点）
                if (vi != depot && temp.ring.size() > 3) {
                    auto it = find(temp.ring.begin(), temp.ring.end(), vi);
                    if (it != temp.ring.end()) {
                        temp.ring.erase(it);
                        temp.in_ring[vi] = false;
                        drop_redistribute(temp, vi); // 自定义删除后处理
                    }
                }
                // 添加vj
                if (!temp.in_ring[vj]) {
                    insert_node(temp, vj);
                }
            }
            else if (vi_in && vj_in) {  // Case 2: 两者都在环
                // 删除vi（确保不是根节点）
                if (vi != depot && temp.ring.size() > 3) {
                    auto it = find(temp.ring.begin(), temp.ring.end(), vi);
                    if (it != temp.ring.end()) {
                        temp.ring.erase(it);
                        temp.in_ring[vi] = false;
                        drop_redistribute(temp, vi);
                    }
                }
            }
            else if (!vi_in && !vj_in) { // Case 3: 两者都不在环
                // 添加vj
                if (!temp.in_ring[vj]) {
                    insert_node(temp, vj);
                }
            }
        }
        return temp;
    }

    // 新增辅助函数：删除节点后的处理
    void drop_redistribute(Solution& s, int deleted_node) {

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
        affected_nodes.insert(deleted_node); // 被删节点自身

        // 重新分配
        for (int u : affected_nodes) {
            if (s.in_ring[u]) continue;

            int new_assign = find_closest_ring_node(s, u);
            if (new_assign == -1) continue;

            // 更新分配关系
            s.assign_cost += graph.assign_cost[u][new_assign] -
                graph.assign_cost[u][s.assignments[u]];
            s.assignments[u] = new_assign;
        }
    }
    // 新的局部搜索过程（所有领域）
    Solution local_search_v2(Solution s) {
        Solution current = s;
        evaluate(current);
        vector<Solution> H;
        do {
            // 生成四个邻域的解集
            vector<Solution> Hadd = generate_add_neighbors(current);
            vector<Solution> Hdrop = generate_drop_neighbors(current);
            vector<Solution> Hadd_drop = generate_add_drop_neighbors(current);
            vector<Solution> Hopt = generate_opt_neighbors(current);

            // 合并邻域解集
            H.clear();
            H.insert(H.end(), Hadd.begin(), Hadd.end());
            H.insert(H.end(), Hdrop.begin(), Hdrop.end());
            H.insert(H.end(), Hadd_drop.begin(), Hadd_drop.end());
            H.insert(H.end(), Hopt.begin(), Hopt.end());

            if (!H.empty()) {
                // 选择最优候选解
                auto best_it = min_element(H.begin(), H.end(),
                    [](const Solution& a, const Solution& b) {
                        return a.total_cost() < b.total_cost();
                    });
                current = *best_it;
                evaluate(current);
            }
        } while (!H.empty());

        return current;
    }

    // 生成添加邻域（Nadd）
    vector<Solution> generate_add_neighbors(Solution s) {
        vector<Solution> neighbors;
        vector<int> candidates;

        // 收集所有未使用的节点
        for (int i = 0; i < graph.nodes.size(); ++i) {
			if (!s.in_ring[i] && i != graph.depot) {
				candidates.push_back(i);
			}
        }
            
            // 为每个候选节点生成可能的插入解
            for (int v : candidates) {
                Solution neighbor = s;
                insert_node(neighbor, v);
                evaluate(neighbor);
                if (neighbor.total_cost() < s.total_cost())
                    neighbors.push_back(neighbor);
            }
            return neighbors;
    }

    // 生成删除邻域（Ndrop）
    vector<Solution> generate_drop_neighbors(Solution s) {
        vector<Solution> neighbors;

        // 排除根节点和最小环长约束
        if (s.ring.size() <= 3) return neighbors;

        // 尝试删除每个非根节点
        for (int i = 1; i < s.ring.size(); ++i) {
            Solution neighbor = s;
            int deleted_node = neighbor.ring[i];
            drop_and_redistribute(neighbor);
            evaluate(neighbor);
            if (neighbor.total_cost() < s.total_cost())
                neighbors.push_back(neighbor);
        }
        return neighbors;
    }

    // 生成添加/删除组合邻域
    vector<Solution> generate_add_drop_neighbors(Solution s) {
        vector<Solution> neighbors;
        // 实现组合操作，例如先添加后删除或相反
        // 这里展示一个简单实现：随机添加+删除组合
        for (int i = 0; i < 5; ++i) { // 生成5个候选
            Solution neighbor = s;
            if (neighbor.ring.size() > 3) {
                drop_and_redistribute(neighbor);
                insert_random_node(neighbor);
                evaluate(neighbor);
                if (neighbor.total_cost() < s.total_cost())
                    neighbors.push_back(neighbor);
            }
        }
        return neighbors;
    }

    // 生成优化邻域（2-opt等）
    vector<Solution> generate_opt_neighbors(Solution s) {
        vector<Solution> neighbors;

        // 1. 增强版优化（1个高质量解）
        Solution enhanced_neighbor = s;
        randomized_two_opt(enhanced_neighbor);
        evaluate(enhanced_neighbor);
        if (enhanced_neighbor.total_cost() < s.total_cost())
            neighbors.push_back(enhanced_neighbor);

        // 2. 保留原随机优化（3个随机解）
        for (int i = 0; i < 3; ++i) {
            Solution rand_neighbor = s;
            two_opt(rand_neighbor);
            evaluate(rand_neighbor);
            if (rand_neighbor.total_cost() < s.total_cost())
                neighbors.push_back(rand_neighbor);
        }

        return neighbors;
    }

    // 邻域操作1：随机插入未使用的节点
    void insert_random_node(Solution& s) {
        int v = random_unused_node(s);         // 随机选择不在环中的节点
        if (v == -1) return;                  // 无可用节点时退出
        insert_node(s, v);                     // 执行插入操作
    }

    // 邻域操作2：随机删除环中节点（非根节点）
    // 修改后的完整代码段
    void drop_and_redistribute(Solution& s) {
        if (s.ring.size() <= 2) return;

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
                s.assign_cost -= graph.assign_cost[u][s.assignments[u]];
            }
            s.assign_cost += graph.assign_cost[u][new_assign];

            // 更新分配关系
            s.assignments[u] = new_assign;

        }

    }

    // 邻域操作3：交换环内和环外的两个节点
    void swap_nodes(Solution& s) {
        if (s.ring.size() < 2) return;
        int i = 1 + rand() % (s.ring.size() - 1); // 随机选择环内位置
        int j = 1 + rand() % (s.ring.size() - 1);
        swap(s.ring[i], s.ring[j]);           // 交换节点位置
    }

    // 邻域操作4：2-opt边交换优化
    void two_opt(Solution& s) {
        if (s.ring.size() < 4) return;        // 至少需要4个节点进行边交换
        int i = 1 + rand() % (s.ring.size() - 3); // 随机选择起始位置
        int j = i + 1 + rand() % (s.ring.size() - i - 1); // 选择结束位置
        reverse(s.ring.begin() + i, s.ring.begin() + j + 1); // 反转子路径
    }

    void enhanced_two_opt(Solution& s) {
        double best_gain = 0;
        int best_i = -1, best_j = -1;

        // 遍历所有可能的边交换
        for (int i = 0; i< s.ring.size() - 1; ++i) {
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

    void randomized_two_opt(Solution& s) {
        const int max_attempts = 50; // 最大尝试次数
        int attempts = 0;
        bool improved = true;

        while (improved || attempts < max_attempts) {
            improved = false;

            // 随机选择两个不同的边
            int i = rand() % (s.ring.size() - 1);
            int j = rand() % (s.ring.size() - 1);
            if (abs(i - j) < 2) continue; // 避免相邻边

            // 确保i < j
            if (i > j) swap(i, j);
            j++; // 调整为j+1的索引

            // 计算当前成本
            double original = graph.routing_cost[s.ring[i]][s.ring[i + 1]]
                + graph.routing_cost[s.ring[j]][s.ring[(j + 1) % s.ring.size()]];

            // 计算交换后成本
            double modified = graph.routing_cost[s.ring[i]][s.ring[j]]
                + graph.routing_cost[s.ring[i + 1]][s.ring[(j + 1) % s.ring.size()]];

            // 如果成本降低
            if (modified < original) {
                reverse(s.ring.begin() + i + 1, s.ring.begin() + j + 1);
                improved = true;
                attempts = 0; // 重置尝试计数器
            }
            else {
                attempts++;
            }
        }
    }

    // 随机选择未使用的节点（返回-1表示无可用节点）
    int random_unused_node(const Solution& s) {
        vector<int> candidates;
        for (int i = 0; i < graph.nodes.size(); ++i) {
            if (!s.in_ring[i] && i != graph.depot) {
                candidates.push_back(i);
            }
        }
        if (candidates.empty()) return -1;
        return candidates[rand() % candidates.size()];
    }

    // 评估解决方案的总成本
    void evaluate(Solution& s) {
        // 计算路由成本（环路径总长度）
        s.routing_cost = 0;
        for (int i = 0; i < s.ring.size(); ++i) {
            int j = (i + 1) % s.ring.size();
            s.routing_cost += alpha*graph.routing_cost[s.ring[i]][s.ring[j]];
        }

        // 计算分配成本（所有非环节点到最近环节点的成本）
        s.assign_cost = 0;
        for (int u = 0; u < graph.nodes.size(); ++u) {
            if (!s.in_ring[u]) {              // 只考虑非环节点
                s.assign_cost += (10-alpha)*graph.assign_cost[u][s.assignments[u]];
            }
        }
    }

    // 辅助函数：查找离节点u最近的环节点
    int find_closest_ring_node(Solution& s, int u) {
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

    // 插入节点到环的最佳位置
    void insert_node(Solution& s, int v) {
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
                s.assign_cost += graph.assign_cost[u][v] - graph.assign_cost[u][s.assignments[u]];
                s.assignments[u] = v;          // 更新分配关系
            }
        }
    }
    
};



RSPGraph parseTSPLIB(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "无法打开文件: " << filename << endl;
        exit(1);
    }

    RSPGraph graph;
    string line;
    int dimension = 0;
    bool node_section = false;

    while (getline(file, line)) {
        // 清理行首尾空格
        line.erase(line.find_last_not_of(" \t") + 1);
        line.erase(0, line.find_first_not_of(" \t"));

        // 统一转为小写处理
        string lowerLine = line;
        transform(lowerLine.begin(), lowerLine.end(), lowerLine.begin(), ::tolower);

        if (lowerLine.find("dimension") != string::npos) {
            size_t colon = line.find(':');
            if (colon != string::npos) {
                string dimStr = line.substr(colon + 1);
                // 过滤非数字字符
                dimStr.erase(remove_if(dimStr.begin(), dimStr.end(),
                    [](char c) { return !isdigit(c); }), dimStr.end());
                if (!dimStr.empty()) {
                    dimension = stoi(dimStr);
                    graph.initialize(dimension);
                    graph.nodes.resize(dimension);
                }
            }
        }
        else if (lowerLine.find("node_coord_section") != string::npos) {
            node_section = true;
			continue; // 跳过该行
        }
        else if (line.find("EOF") != string::npos) {
            break;
        }

        if (node_section) {
            istringstream iss(line);
            int id;
            double x, y;
            iss >> id >> x >> y;
            id--; // 转换为0-based索引
            graph.nodes[id] = { x, y };
        }
    }

    // 计算距离矩阵
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            double dx = graph.nodes[i].first - graph.nodes[j].first;
            double dy = graph.nodes[i].second - graph.nodes[j].second;
            double dist = hypot(dx, dy);
            // 修改为四舍五入后的整数距离
            int dist_int = static_cast<int>(round(dist));  // 添加四舍五入转换[1,3](@ref)

            // 更新距离矩阵
            graph.routing_cost[i][j] = dist_int;
            graph.assign_cost[i][j] = dist_int;
        }
    }

    return graph;
}
int main() {
    const int ITER_PER_RUN = 50;  // 每次运行迭代次数
    Solution global_best;          // 全局最优解
    double min_cost = INF;         // 最小成本跟踪
	int max_iter = 100;           // 最大迭代次数
    RSPGraph graph = parseTSPLIB("E:/code/c/tsplib/eil51.tsp"); //berlin52 eil51
    for (int run = 0; run < max_iter; ++run) {
		cout << "\n======= 运行 #" << run + 1 << " =======\n";
        // 创建新求解器实例
        AMNSSolver solver(graph, 3);

        // 执行迭代优化
        solver.solve(ITER_PER_RUN);

        // 获取当前最优解
        Solution current = solver.getBestSolution();

        // 更新全局最优
        if (current.total_cost() < min_cost) {
            global_best = current;
            min_cost = current.total_cost();
        }
		cout << "迭代 #" << run << " 当前最优解成本: " << min_cost
			<< " (路由: " << global_best.routing_cost
			<< ", 分配: " << global_best.assign_cost << ")" << endl;
    }
    // 输出最终结果
    cout << "\n======= 全局最优解 =======";
    cout << "\n总成本: " << global_best.total_cost()
        << " (路由: " << global_best.routing_cost
        << ", 分配: " << global_best.assign_cost << ")\n";

    // 输出环路径
    cout << "\n环路径: 0";
    for (int node : global_best.ring) {
        if (node != 0) cout << " -> " << node; // 跳过重复的0
    }
    cout << " -> 0" << endl;

    // 输出分配关系
    cout << "\n分配关系:" << endl;
    for (int u = 0; u < graph.nodes.size(); ++u) {
        if (!global_best.in_ring[u]) {
            cout << "节点" << u << "\t→ 环节点" << global_best.assignments[u] << endl;
        }
    }

    return 0;
}