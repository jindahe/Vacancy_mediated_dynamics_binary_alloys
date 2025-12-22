#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

// 周期性边界索引
inline int get_idx(int x, int y, int L) {
    return ((x + L) % L) * L + ((y + L) % L);
}

// 初始化晶格 (Requirement a)
void initialize_lattice(vector<int>& lattice, int L, int& vacancy_pos, mt19937& g) {
    int N = L * L;
    vector<int> atoms;
    for (int i = 0; i < N / 2; ++i) { atoms.push_back(1); atoms.push_back(-1); }
    shuffle(atoms.begin(), atoms.end(), g);
    lattice = atoms;
    uniform_int_distribution<int> dist(0, N - 1);
    vacancy_pos = dist(g);
    lattice[vacancy_pos] = 0; // 注入空位
}

// 局部能量计算
double get_local_energy(const vector<int>& lattice, int idx, int L, double J) {
    if (lattice[idx] == 0) return 0.0;
    int x = idx / L, y = idx % L;
    int neighbors[] = {get_idx(x+1, y, L), get_idx(x-1, y, L), get_idx(x, y+1, L), get_idx(x, y-1, L)};
    double energy = 0;
    for (int n : neighbors) {
        if (lattice[n] != 0) energy += (lattice[idx] == lattice[n]) ? -J : J;
    }
    return energy;
}

// 全局能量计算
double get_total_energy(const vector<int>& lattice, int L, double J) {
    double total = 0;
    for (int i = 0; i < L * L; ++i) total += get_local_energy(lattice, i, L, J);
    return total / 2.0;
}

// 计算对关联函数 C(r) (符合 Project 15.45 Requirement b)
vector<double> calculate_C_r(const vector<int>& lattice, int L) {
    int max_r = L / 2;
    vector<double> C(max_r + 1, 0.0);
    vector<long long> count(max_r + 1, 0);

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int idx = i * L + j;
            if (lattice[idx] == 0) continue;
            for (int r = 1; r <= max_r; ++r) {
                int idx_h = i * L + (j + r) % L;
                if (lattice[idx_h] != 0) { C[r] += lattice[idx] * lattice[idx_h]; count[r]++; }
                int idx_v = ((i + r) % L) * L + j;
                if (lattice[idx_v] != 0) { C[r] += lattice[idx] * lattice[idx_v]; count[r]++; }
            }
        }
    }
    for (int r = 1; r <= max_r; ++r) if (count[r] > 0) C[r] /= count[r];
    C[0] = 1.0; 
    return C;
}

// 存储 C(r) 到 output 文件夹
void save_C_r(const vector<double>& Cr, int mcs) {
    // 将路径指向 output 目录
    string filename = "../output/Cr_t_" + to_string(mcs) + ".txt"; 
    ofstream out(filename);
    if (out.is_open()) {
        for (size_t r = 0; r < Cr.size(); ++r) {
            out << r << "\t" << fixed << setprecision(6) << Cr[r] << "\n";
        }
        out.close();
    } else {
        cerr << "无法打开文件: " << filename << " (请确保 output 文件夹已创建)" << endl;
    }
}

// 存储晶格快照到 output 文件夹 (符合 Project 15.45 Requirement b)
void save_lattice(const vector<int>& lattice, int L, int mcs) {
    // 将路径指向 output 目录
    string filename = "../output/lattice_t_" + to_string(mcs) + ".txt";
    ofstream out(filename);
    if (out.is_open()) {
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                out << lattice[i * L + j] << " ";
            }
            out << "\n";
        }
        out.close();
    }
}

#endif