#include "functions.h"
#include <chrono>
#include <math.h>

using namespace std;

int main() {
    const int L = 128;
    const double Tc = 2.269, T = Tc / 2.0, J = 1.0;
    const int num_mc = pow(2,20); //33000;
    
    mt19937 gen(chrono::steady_clock::now().time_since_epoch().count());
    vector<int> lattice;
    int vacancy_pos;
    initialize_lattice(lattice, L, vacancy_pos, gen);
    
    double current_energy = get_total_energy(lattice, L, J);
    ofstream r_file("../output/t_vs_R.txt");
    ofstream time_log("../output/time_log_1.txt");
    
    auto total_start = chrono::high_resolution_clock::now();
    
    for (int mcs = 0; mcs <= num_mc; ++mcs) {
        auto step_start = chrono::high_resolution_clock::now();
        // Monte Carlo 步 (空位交换逻辑)
        for (int step = 0; step < L * L; ++step) {
            int dir = uniform_int_distribution<int>(0, 3)(gen);
            int vx = vacancy_pos / L, vy = vacancy_pos % L;
            int nx = vx, ny = vy;
            if (dir == 0) nx++; else if (dir == 1) nx--; else if (dir == 2) ny++; else ny--;
            int neighbor_idx = get_idx(nx, ny, L);

            double e_old = get_local_energy(lattice, neighbor_idx, L, J);
            swap(lattice[vacancy_pos], lattice[neighbor_idx]);
            double e_new = get_local_energy(lattice, vacancy_pos, L, J);
            double dE = e_new - e_old;

            if (dE <= 0 || uniform_real_distribution<double>(0, 1)(gen) < exp(-dE / T)) {
                current_energy += dE;
                vacancy_pos = neighbor_idx;
            } else {
                swap(lattice[vacancy_pos], lattice[neighbor_idx]);
            }
        }

        auto step_end = chrono::high_resolution_clock::now();
        double step_time = chrono::duration<double>(step_end - step_start).count();

        bool is_sample_step = ( (mcs > 0 && (mcs & (mcs - 1)) == 0) || mcs == 0 || mcs == num_mc );

        // 定期记录 (t = 2^n) 符合 Requirement b/c
        if (is_sample_step) {
            // 1. 计算 R 并写入 output/t_vs_R.txt (符合 Requirement c)
            double R_energy = 2.0 / ((current_energy / (L * L)) / J + 2.0);
            r_file << mcs << "\t" << R_energy << endl;

            // 2. 计算并存储 C(r) 到 output/Cr_t_X.txt (符合 Requirement c)
            vector<double> Cr = calculate_C_r(lattice, L);
            save_C_r(Cr, mcs);

            // 3. 存储晶格配置到 output/lattice_t_X.txt
            save_lattice(lattice, L, mcs);

            cout << "MCS: " << mcs << " | Time: " << step_time << " s | R: " << R_energy << endl;
        } else {
            if (mcs % 1000 == 0) {
                cout << "MCS: " << mcs << " | Time: " << step_time << " s" << endl;
            }
        }

        if (mcs % 1000 == 0) {
            time_log << mcs << "\t" << step_time << "\n";
        }
    }
    
    auto total_end = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(total_end - total_start).count();
    cout << "\nTotal simulation time: " << total_time << " s (" << total_time/60.0 << " min)" << endl;
    time_log.close();
    r_file.close();
    return 0;
}