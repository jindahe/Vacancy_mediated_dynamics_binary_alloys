#ifndef FUNCTIONS_KAWASAKI_H
#define FUNCTIONS_KAWASAKI_H

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

// 【修改点】：初始化满晶格，无空位 (50% A, 50% B)
void initialize_lattice_kawasaki(vector<int>& lattice, int L, mt19937& g) {
    int N = L * L;
    lattice.resize(N);
    // 前一半填 1，后一半填 -1
    for (int i = 0; i < N; ++i) {
        lattice[i] = (i < N / 2) ? 1 : -1;
    }
    shuffle(lattice.begin(), lattice.end(), g);
}

// 计算单个格点周围的邻居自旋和 (不包括自身相互作用)
// 用于计算 Kawasaki 交换的能量差
int get_neighbor_sum(const vector<int>& lattice, int idx, int L) {
    int x = idx / L;
    int y = idx % L;
    
    int n1 = lattice[get_idx(x + 1, y, L)];
    int n2 = lattice[get_idx(x - 1, y, L)];
    int n3 = lattice[get_idx(x, y + 1, L)];
    int n4 = lattice[get_idx(x, y - 1, L)];
    
    return n1 + n2 + n3 + n4;
}

// 全局能量计算 (用于 R 的计算)
double get_total_energy(const vector<int>& lattice, int L, double J) {
    double energy = 0;
    for (int i = 0; i < L * L; ++i) {
        int neighbors_sum = get_neighbor_sum(lattice, i, L);
        energy += -J * lattice[i] * neighbors_sum;
    }
    return energy / 2.0; // 每条键计算了两次
}

// 对关联函数 (保持不变)
vector<double> calculate_C_r(const vector<int>& lattice, int L) {
    int max_r = L / 2;
    vector<double> C(max_r + 1, 0.0);
    vector<long long> count(max_r + 1, 0);

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int idx = i * L + j;
            for (int r = 1; r <= max_r; ++r) {
                int idx_h = i * L + (j + r) % L;
                C[r] += lattice[idx] * lattice[idx_h]; count[r]++;
                int idx_v = ((i + r) % L) * L + j;
                C[r] += lattice[idx] * lattice[idx_v]; count[r]++;
            }
        }
    }
    for (int r = 1; r <= max_r; ++r) if (count[r] > 0) C[r] /= count[r];
    C[0] = 1.0; 
    return C;
}

// 保存数据辅助函数
void save_lattice(const vector<int>& lattice, int L, int mcs) {
    string filename = "../output/lattice_kawasaki_t_" + to_string(mcs) + ".txt";
    ofstream out(filename);
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) out << lattice[i * L + j] << " ";
        out << "\n";
    }
}
void save_C_r(const vector<double>& Cr, int mcs) {
    string filename = "../output/Cr_kawasaki_t_" + to_string(mcs) + ".txt";
    ofstream out(filename);
    for (size_t r = 0; r < Cr.size(); ++r) out << r << "\t" << Cr[r] << "\n";
}

#endif