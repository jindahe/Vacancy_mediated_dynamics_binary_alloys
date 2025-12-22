#ifndef FUNCTIONS_3D_H
#define FUNCTIONS_3D_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

// 三维索引计算：x*L^2 + y*L + z
inline int get_idx_3d(int x, int y, int z, int L) {
    return ((x + L) % L) * L * L + ((y + L) % L) * L + ((z + L) % L);
}

// 三维初始化 (Requirement a)
void initialize_lattice_3d(vector<int>& lattice, int L, int& v_pos, mt19937& g) {
    int N = L * L * L;
    vector<int> atoms;
    atoms.reserve(N);
    for (int i = 0; i < N / 2; ++i) atoms.push_back(1);
    for (int i = 0; i < N - (N/2); ++i) atoms.push_back(-1);
    
    shuffle(atoms.begin(), atoms.end(), g);
    lattice = atoms;
    
    uniform_int_distribution<int> dist(0, N - 1);
    v_pos = dist(g);
    lattice[v_pos] = 0; // 注入空位
}

// 三维局部能量计算：检查 6 个邻居
double get_local_energy_3d(const vector<int>& lattice, int idx, int L, double J) {
    if (lattice[idx] == 0) return 0.0;
    int x = idx / (L * L);
    int y = (idx / L) % L;
    int z = idx % L;
    
    int neighbors[] = {
        get_idx_3d(x+1, y, z, L), get_idx_3d(x-1, y, z, L),
        get_idx_3d(x, y+1, z, L), get_idx_3d(x, y-1, z, L),
        get_idx_3d(x, y, z+1, L), get_idx_3d(x, y, z-1, L)
    };
    
    double energy = 0;
    for (int n : neighbors) {
        if (lattice[n] != 0) {
            energy += (lattice[idx] == lattice[n]) ? -J : J;
        }
    }
    return energy;
}

double get_total_energy_3d(const vector<int>& lattice, int L, double J) {
    double total = 0;
    for (int i = 0; i < L * L * L; ++i) total += get_local_energy_3d(lattice, i, L, J);
    return total / 2.0;
}

// 三维对关联函数 (Requirement b)
vector<double> calculate_C_r_3d(const vector<int>& lattice, int L) {
    int max_r = L / 2;
    vector<double> C(max_r + 1, 0.0);
    vector<long long> count(max_r + 1, 0);

    for (int x = 0; x < L; ++x) {
        for (int y = 0; y < L; ++y) {
            for (int z = 0; z < L; ++z) {
                int idx = get_idx_3d(x, y, z, L);
                if (lattice[idx] == 0) continue;

                for (int r = 1; r <= max_r; ++r) {
                    int neighbors[] = {
                        get_idx_3d(x+r, y, z, L),
                        get_idx_3d(x, y+r, z, L),
                        get_idx_3d(x, y, z+r, L)
                    };
                    for (int n : neighbors) {
                        if (lattice[n] != 0) {
                            C[r] += lattice[idx] * lattice[n];
                            count[r]++;
                        }
                    }
                }
            }
        }
    }
    for (int r = 1; r <= max_r; ++r) if (count[r] > 0) C[r] /= count[r];
    C[0] = 1.0;
    return C;
}

void save_lattice_3d(const vector<int>& lattice, int L, int mcs) {
    // 路径指向 output_3d 目录
    string filename = "../output_3d/lattice_3d_t_" + to_string(mcs) + ".txt";
    ofstream out(filename);
    if (out.is_open()) {
        // 为了节省空间，我们按顺序写入，Python 端再 reshape
        for (int val : lattice) {
            out << val << " ";
        }
        out.close();
        cout << "3D Lattice saved for t = " << mcs << endl;
    }
}

#endif