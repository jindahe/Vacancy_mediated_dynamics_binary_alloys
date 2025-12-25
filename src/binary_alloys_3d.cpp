#include "functions_3d.h"
#include <chrono>
#include <filesystem>
#include <math.h>

using namespace std;
namespace fs = std::filesystem;

int main() {
    
    // 1. 参数设置
    const int L = 64; // 3D 计算量大，L建议先设小一点
    const int N = L * L * L;
    const double Tc_3d = 4.51; // 3D 临界温度
    const double T = Tc_3d / 2.0;
    const double J = 1.0;
    const int num_mc = pow(2,20);

    // 2. 初始化
    mt19937 gen(chrono::steady_clock::now().time_since_epoch().count());
    vector<int> lattice;
    int v_pos;
    initialize_lattice_3d(lattice, L, v_pos, gen);
    
    double current_energy = get_total_energy_3d(lattice, L, J);
    ofstream r_file("../output_3d/t_vs_R.txt");
    ofstream time_log("../output_3d/time_log.txt");

    auto total_start = chrono::high_resolution_clock::now();

    // 3. 模拟循环
    for (int mcs = 0; mcs <= num_mc; ++mcs) {
        auto step_start = chrono::high_resolution_clock::now();
        for (int step = 0; step < N; ++step) {
            int dir = uniform_int_distribution<int>(0, 5)(gen); // 6个方向
            int vx = v_pos / (L * L);
            int vy = (v_pos / L) % L;
            int vz = v_pos % L;

            int nx = vx, ny = vy, nz = vz;
            if (dir == 0) nx++; else if (dir == 1) nx--;
            else if (dir == 2) ny++; else if (dir == 3) ny--;
            else if (dir == 4) nz++; else nz--;

            int n_idx = get_idx_3d(nx, ny, nz, L);

            double e_old = get_local_energy_3d(lattice, n_idx, L, J);
            swap(lattice[v_pos], lattice[n_idx]);
            double e_new = get_local_energy_3d(lattice, v_pos, L, J);
            double dE = e_new - e_old;

            if (dE <= 0 || uniform_real_distribution<double>(0,1)(gen) < exp(-dE/T)) {
                current_energy += dE;
                v_pos = n_idx;
            } else {
                swap(lattice[v_pos], lattice[n_idx]);
            }
        }

        auto step_end = chrono::high_resolution_clock::now();
        double step_time = chrono::duration<double>(step_end - step_start).count();

        bool is_sample_step = ( (mcs > 0 && (mcs & (mcs - 1)) == 0) || mcs == 0 || mcs == num_mc );

        // 4. 定期采样 (t = 2^n)
        if (is_sample_step) {
            double R = 3.0 / ((current_energy / N) / J + 3.0);
            r_file << mcs << "\t" << R << endl;
            
            vector<double> Cr = calculate_C_r_3d(lattice, L);
            string c_name = "../output_3d/Cr_t_" + to_string(mcs) + ".txt";
            ofstream c_out(c_name);
            for(size_t r=0; r<Cr.size(); ++r) c_out << r << "\t" << Cr[r] << "\n";

            save_lattice_3d(lattice, L, mcs);
            
            cout << "MCS: " << mcs << " | Time: " << step_time << " s | R: " << R << endl;
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