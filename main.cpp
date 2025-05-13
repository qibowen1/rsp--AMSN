#include "main.h"
#include <unordered_map>
namespace fs = std::filesystem;
using namespace std;

const double INF = numeric_limits<double>::max(); // 定义无穷大常量


/*
针对Ring Star Problem（RSP），​​自适应多阶段邻域搜索算法（Adaptive Multi-phase Neighborhood Search, AMNS）​​，
结合动态邻域选择和分配优化策略，以平衡全局探索与局部 开发能力。
*/

int main() {
    std::srand(std::time(nullptr));  // 用当前时间初始化种子
	const int MAX_ITER = 25; // 每个文件的最大迭代次数
	const int MAX_SEARCH_ITER = 1; // 每个初始解的迭代次数
    const vector<int> ALPHAS = { 5,7,9 };  // 需要测试的alpha值列表
	double RCL_ratio = 0.8; // RCL比例  
	int neighborhood_types = 5; // 邻域类型数量    
    std::unordered_map<BenchmarkKey, double> benchmark_map;
    const string TSP_DIR = "E:/code/c/c++/input_data/tsp_part"; // 修改为你的TSP目录路径
	const string TSP_TOUR_DIR = "E:/code/c/c++/out_tour"; // alpha=3时的LKH输出目录
    vector<Result> all_results;  // 存储所有结果
    load_benchmark_data("E:/code/c/c++/opt.txt", benchmark_map);

    // 生成带时间的文件名
    auto now = chrono::system_clock::now();
    time_t now_time = chrono::system_clock::to_time_t(now);
    struct tm local_time;
    localtime_s(&local_time, &now_time);  // Windows安全版本
    // 创建日期格式的文件夹名 (年月日)
    ostringstream dirname_ss;
    dirname_ss << "E:/code/c/c++/rsp-result/results_"
        << put_time(&local_time, "%Y%m%d");  // 年月日

    // 创建文件夹
    filesystem::path dir_path(dirname_ss.str());
    if (!filesystem::exists(dir_path)) {
        filesystem::create_directories(dir_path);
    }
    // 创建带时间的文件名
    ostringstream filename_ss;
    filename_ss << dirname_ss.str() << "/results_"
        << put_time(&local_time, "%m%d_%H%M")  // 月日_时分
        << ".csv";

    // 创建并打开输出文件
    ofstream out_file(filename_ss.str());
    if (!out_file) {
        cerr << "无法创建输出文件！" << endl;
        return 1;
    }
	out_file << "文件名, Alpha, 总成本, 路由成本, 分配成本, 运行时间(s), 差额(%)\n";

    // 遍历目录中的所有TSP文件
    for (const auto& entry : fs::directory_iterator(TSP_DIR)) {
        // 跳过非普通文件或非.tsp后缀的文件
        if (!fs::is_regular_file(entry.path()) || entry.path().extension() != ".tsp")
            continue;
        // 提取纯文件名（包含后缀）
        string filename_output = entry.path().filename().string(); // 例如"st70.tsp"
        string filename = entry.path().string();// 完整路径
        cout << "\n======= 处理文件: " << filename_output << " =======\n";
        // 解析当前TSP文件
        RSPGraph graph = parseTSPLIB(filename);
        Solution global_best;
        double min_cost = INF;
        // 3. 获取文件名stem（例如从"st70.tsp"得到"st70"）
        string filename_stem = entry.path().stem().string();
		for (int alpha : ALPHAS) {
			cout << "\n======= 测试 alpha = " << alpha << " =======\n";
            Solution global_best;
            double min_cost = INF;
            // 5. 修改运行循环
            bool early_stop = false;
            // 4. 查找基准最优值
            BenchmarkKey key{ filename_stem, alpha };
            auto it = benchmark_map.find(key);
            if (it == benchmark_map.end()) {
                cerr << "警告：未找到基准数据 " << filename_stem << " alpha=" << alpha << endl;
                continue;
            }
            double benchmark_opt = it->second;
            // 记录开始时间
            auto start_time = chrono::high_resolution_clock::now();
            // 对当前文件执行多次独立求解
                AMNSSolver solver(graph, alpha, MAX_SEARCH_ITER,RCL_ratio,neighborhood_types); // alpha
                solver.solve(benchmark_opt, MAX_ITER);

                // 更新全局最优解
                Solution current = solver.getBestSolution();
                if (current.total_cost() < min_cost) {
                    global_best = current;
                    min_cost = current.total_cost();
                }
            // 计算运行时间
            auto end_time = chrono::high_resolution_clock::now();
             auto duration = chrono::duration_cast<chrono::duration<double>>(end_time - start_time);
            double duration_sec = duration.count();
			double offset = 100*(global_best.total_cost() - benchmark_opt)/ benchmark_opt;
            // 输出最终结果
            
            cout << "\n​**​* 文件 " << filename_output <<"alpha="<<alpha << " 的最优解 ​**​*\n";
            cout << "总成本: " << global_best.total_cost()
                << "\n路由成本: " << global_best.routing_cost
                << "\n分配成本: " << global_best.assign_cost
                << "\n运行时间: " << duration_sec << " 秒"
                << "\n差额(%):" << offset << endl;

            // 存储结果
            all_results.push_back({
                filename_output,
                alpha,
                global_best.total_cost(),
                global_best.routing_cost,
                global_best.assign_cost,
                duration_sec,  // 修改字段名
				offset
                });

            // 立即写入当前结果到文件

            out_file << filename_output << ","
                << alpha << ","
                << global_best.total_cost() << ","
                << global_best.routing_cost << ","
                << global_best.assign_cost << ","
                << duration_sec << ",\""
                << offset << "\",\"";
            out_file << "\"\n";  // 结束引号和换行
		}
    }
    // 统一输出结果
    cout << "\n\n======= 最终结果汇总 =======";
    for (const auto& res : all_results) {
        cout << "\n\n文件: " << res.filename
            << "\nAlpha: " << res.alpha
            << "\n总成本: " << res.total_cost
            << "\n路由成本: " << res.routing_cost
            << "\n分配成本: " << res.assign_cost
            << "\n运行时间: " << res.duration_s << "s";
    }
    // 关闭文件
    out_file.close();

    // 可选：在控制台显示保存路径
    cout << "\n结果已保存至：results.csv" << endl;
    cout << "\n\n======= 汇总结束 =======\n";
    return 0;
}