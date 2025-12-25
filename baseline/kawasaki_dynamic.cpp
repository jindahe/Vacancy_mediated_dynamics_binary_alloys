#include "functions_kawasaki.h"
#include <chrono>
#include <iostream>
#include <filesystem>
#include <math.h>

using namespace std;
namespace fs = std::filesystem;

int main() {

    const int L = 128; // 可以跑 128，因为这里也用了局部能量优化
    const double Tc = 2.269;
    const double T = Tc / 2.0;
    const double J = 1.0;
    const int num_mc = pow(2,20); 

    mt19937 gen(chrono::steady_clock::now().time_since_epoch().count());
    
    // 初始化无空位晶格
    vector<int> lattice;
    initialize_lattice_kawasaki(lattice, L, gen);

    double current_total_energy = get_total_energy(lattice, L, J);

    ofstream r_file("../output/t_vs_R_kawasaki.txt");
    ofstream time_log("../output/time_log_kawasaki.txt");

    cout << "Starting Kawasaki Dynamics Simulation (L=" << L << ")..." << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // 预计算指数表 (Optimization)
    // Kawasaki dE 可能范围较大，简单起见直接计算或稍微查表
    // 这里相邻交换涉及的邻居数较少，dE 取值有限，仍可用 exp 直接算

    for (int mcs = 0; mcs <= num_mc; ++mcs) {
        auto step_start = chrono::high_resolution_clock::now();

        // Monte Carlo Step: 尝试 N 次交换
        for (int step = 0; step < L * L; ++step) {
            // 1. 随机选点 A
            int idx1 = uniform_int_distribution<int>(0, L * L - 1)(gen);
            int x1 = idx1 / L;
            int y1 = idx1 % L;

            // 2. 随机选邻居 B
            int dir = uniform_int_distribution<int>(0, 3)(gen);
            int x2 = x1, y2 = y1;
            if (dir == 0) x2++; else if (dir == 1) x2--; 
            else if (dir == 2) y2++; else y2--;
            int idx2 = get_idx(x2, y2, L);

            // 3. 仅当自旋不同时才考虑交换
            if (lattice[idx1] != lattice[idx2]) {
                // 计算能量差 dE (Local Energy Update)
                // 技巧：交换前的能量贡献只需考虑 idx1 和 idx2 周围的邻居
                // 注意：idx1 和 idx2 互为邻居，这条键的能量在交换前后不变 (-1 * 1 = -1, 1 * -1 = -1)
                // 所以只需计算它们与其他 3 个邻居的相互作用
                
                int s1 = lattice[idx1];
                int s2 = lattice[idx2];
                
                // 获取周围邻居和（包含对方，后面抵消或直接算）
                int sum1 = get_neighbor_sum(lattice, idx1, L); // 含 idx2
                int sum2 = get_neighbor_sum(lattice, idx2, L); // 含 idx1

                // E_local = -J * s * sum_neighbors
                // E_old_pair = -J * s1 * sum1 - J * s2 * sum2
                // 注意：直接相加会把 idx1-idx2 键算两次，但我们只关心差值
                // 更精确的 dE 计算：
                // dE = E_new - E_old
                // 交换后，位置 idx1 变成 s2，位置 idx2 变成 s1
                // 邻居们的状态没变
                
                // 为了避免重复计算键，我们只看与“非互为邻居”的那些邻居的能量变化
                // 简便方法：先算出改变前的局部能量，再算出改变后的
                
                double e_old = -J * (s1 * sum1 + s2 * sum2);
                
                // 模拟交换
                lattice[idx1] = s2;
                lattice[idx2] = s1;
                
                int sum1_new = sum1; // 邻居没变
                int sum2_new = sum2;
                
                // 但注意：sum1 是 idx1 的邻居和。其中有一个邻居是 idx2。
                // 在计算 sum1 时读取的是 lattice[idx2]。
                // 刚才我们把 lattice[idx2] 改成了 s1 (原先是 s2)。
                // 所以 sum1 的值其实变了！这里不能简单复用 sum1。
                // 最稳妥的方法：交换后重新 scan 邻居，或者使用修正公式。
                
                // 让我们用最稳妥的重算法 (只算 4+4=8 次访问，依然 O(1))
                sum1_new = get_neighbor_sum(lattice, idx1, L);
                sum2_new = get_neighbor_sum(lattice, idx2, L);
                
                double e_new = -J * (s2 * sum1_new + s1 * sum2_new);
                double dE = e_new - e_old;

                // Metropolis
                if (dE <= 0 || uniform_real_distribution<double>(0, 1)(gen) < exp(-dE / T)) {
                    // 接受：保留当前的交换状态
                    current_total_energy += dE; 
                    // 这里总能量更新可能因为双重计算有微小误差累积，建议定期重算全局能量
                } else {
                    // 拒绝：换回去
                    lattice[idx1] = s1;
                    lattice[idx2] = s2;
                }
            }
        }

        auto step_end = chrono::high_resolution_clock::now();
        double step_time = chrono::duration<double>(step_end - step_start).count();

        // 采样记录
        if ((mcs > 0 && (mcs & (mcs - 1)) == 0) || mcs == 0 || mcs == num_mc) {
            // 为防止累积误差，定期校准能量 (Kawasaki 容易有 dE 累积误差)
            current_total_energy = get_total_energy(lattice, L, J); 
            
            double R = 2.0 / ((current_total_energy / (L * L)) / J + 2.0);
            
            r_file << mcs << "\t" << R << endl;
            vector<double> Cr = calculate_C_r(lattice, L);
            save_C_r(Cr, mcs);
            save_lattice(lattice, L, mcs);

            auto now = chrono::high_resolution_clock::now();
            double elapsed = chrono::duration<double>(now - start_time).count();
            cout << "MCS: " << mcs << " | Time: " << elapsed << "s | R: " << R << endl;
        }

        // 记录日志
        if (mcs % 1000 == 0) time_log << mcs << "\t" << step_time << endl;
    }

    time_log.close();
    r_file.close();
    return 0;
}