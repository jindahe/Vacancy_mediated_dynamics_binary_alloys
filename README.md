# Vacancy-Mediated Dynamics in Binary Alloys

A Monte Carlo simulation project studying spinodal decomposition in binary alloys through vacancy-mediated dynamics.

## 项目简介 (Project Overview)

本项目通过计算机模拟研究二元合金在快速淬火过程中的**旋节分解**（Spinodal Decomposition）现象。与传统的自旋交换动力学不同，本项目采用更贴近真实物理过程的**空位介导动力学**（Vacancy-Mediated Dynamics）来模拟金属原子的扩散和畴生长过程。

This project simulates **spinodal decomposition** in binary alloys after rapid quenching using Monte Carlo methods. Unlike conventional spin-exchange dynamics, this simulation employs **vacancy-mediated dynamics** to more realistically model atomic diffusion and domain growth in metallic alloys.

## 物理背景 (Physical Background)

当二元合金从高温快速淬火到低温不稳定态时，会发生旋节分解过程，两种金属在合金中分离形成畴结构。Lifshitz 和 Slyozov 预测在长时间极限下，线性畴尺寸随时间按 $R \sim t^{1/3}$ 增长。本项目验证这一标度律，并研究空位介导动力学对畴生长的影响。

When a binary alloy is rapidly quenched from high temperature to a low-temperature unstable state, spinodal decomposition occurs as the two metals separate. Lifshitz and Slyozov predicted that at long times, the linear domain size grows as $R \sim t^{1/3}$. This project validates this scaling law and investigates the effects of vacancy-mediated dynamics on domain growth.

## 主要特性 (Key Features)

- **空位介导动力学**: 单个空位与邻近原子交换，更真实地模拟实际合金中的原子扩散
- **2D 和 3D 模拟**: 支持二维和三维晶格系统
- **对数时间采样**: 在 $t = 2^n$ 时刻采样，高效捕捉跨越多个数量级的演化过程
- **多重畴尺寸测量**: 
  - 基于能量的测量: $R = 2/(\langle E \rangle/N + 2)$
  - 基于对关联函数的测量: $C(r) = \langle s_i s_j \rangle$
- **高性能优化**: O(1) 局部能量更新、查表法、内存优化

## 项目结构 (Project Structure)

```
.
├── README.md                  # 项目说明文档
├── final/                     # 最后提交作业
│   ├── presentation.pdf       
│   ├── presentation.pptx      
│   └── report.pdf
├── discription.md            # 详细项目描述
├── baseline/                 # 基线实现（Kawasaki 动力学）
│   ├── kawasaki_dynamic.cpp
│   └── functions_kawasaki.h
├── src/                      # 主要源代码
│   ├── binary_alloys_2d.cpp # 2D 模拟主程序
│   ├── binary_alloys_3d.cpp # 3D 模拟主程序
│   ├── functions.h          # 2D 工具函数
│   ├── functions_3d.h       # 3D 工具函数
│   ├── plot_*.py            # 数据可视化脚本
│   ├── compare_*.py         # 结果对比脚本
│   └── analyze_results.py   # 结果分析脚本
├── output/                   # 2D 模拟输出数据
│   ├── lattice_t_*.txt      # 晶格快照
│   └── Cr_t_*.txt           # 对关联函数数据
├── output_3d/               # 3D 模拟输出数据
└── thesis/                   # 相关论文和文档
```

## 最终交付物 (Final Deliverables)

- 最终报告: [final/report.pdf](final/report.pdf)
- 汇报pdf格式: [final/presentation.pdf](final/presentation.pdf)
- 汇报ppt格式: [final/presentation.pptx](final/presentation.pptx)

## 编译与运行 (Compilation and Execution)

### 环境要求 (Requirements)

- C++ 编译器 (支持 C++17 标准)
- Python 3.x (用于数据可视化)
- Python 库: numpy, matplotlib

### 编译 (Compilation)

```bash
# 编译 2D 模拟程序
cd src
g++ -std=c++17 -O3 -o binary_alloys_2D binary_alloys_2d.cpp

# 编译 3D 模拟程序
g++ -std=c++17 -O3 -o binary_alloys_3d binary_alloys_3d.cpp

# 编译基线程序（Kawasaki 动力学）
cd ../baseline
g++ -std=c++17 -O3 -o kawasaki_dynamic kawasaki_dynamic.cpp
```

### 运行模拟 (Running Simulations)

```bash
# 运行 2D 模拟
cd src
./binary_alloys_2D

# 运行 3D 模拟
./binary_alloys_3d

# 运行基线模拟
cd ../baseline
./kawasaki_dynamic
```

## 数据分析与可视化 (Data Analysis and Visualization)

### 绘制晶格演化 (Plot Lattice Evolution)

```bash
cd src
python plot_lattice.py        # 2D 晶格
python plot_lattice_3d.py     # 3D 晶格
```

### 分析畴尺寸增长 (Analyze Domain Growth)

```bash
python plot_R_from_Energy.py          # 基于能量的 R(t)
python plot_R_from_pair_cor.py        # 基于对关联函数的 R(t) (2D)
python plot_R_from_pair_cor_3d.py     # 基于对关联函数的 R(t) (3D)
```

### 对关联函数分析 (Correlation Function Analysis)

```bash
python plot_correlation_func.py
```

### 结果对比 (Compare Results)

```bash
python compare_growth.py         # 对比不同方法的增长曲线
python compare_energy.py         # 对比能量演化
python compare_correlation.py    # 对比对关联函数
python compare_snapshot.py       # 对比晶格快照
```

## 核心算法 (Core Algorithms)

### 1. 空位介导动力学 (Vacancy-Mediated Dynamics)

```
初始化:
  - 创建 L×L (或 L×L×L) 晶格
  - 随机放置 50% A 原子和 50% B 原子
  - 随机放置一个空位

Monte Carlo 步骤:
  1. 随机选择空位的一个邻居
  2. 计算交换能量差 ΔE
  3. 使用 Metropolis 准则决定是否接受交换:
     - 如果 ΔE ≤ 0: 接受
     - 如果 ΔE > 0: 以概率 exp(-ΔE/kT) 接受
  4. 更新空位位置和晶格配置
```

### 2. 畴尺寸测量 (Domain Size Measurement)

**方法 1: 能量法**
$$R = \frac{2}{\langle E \rangle / N + 2}$$

**方法 2: 对关联函数法**
$$C(r) = \langle s_i s_j \rangle$$
其中 $s_i = +1$ (A 原子), $s_i = -1$ (B 原子)

畴尺寸 R 对应 $C(r)$ 的第一个零点。

## 性能优化 (Performance Optimizations)

1. **O(1) 局部能量更新**: 只计算受影响的邻居位点能量变化
2. **指数运算查表法**: 预计算 Metropolis 判据中的指数值
3. **连续内存存储**: 使用 1D vector 存储 2D/3D 数据
4. **对数时间采样**: 减少数据存储量同时保留关键信息

## 物理参数 (Physical Parameters)

- **格点尺寸 (Lattice Size)**: L = 50-128 (2D), L = 30-50 (3D)
- **相互作用能 (Interaction Energy)**: J = 1.0
- **临界温度 (Critical Temperature)**: $T_c \approx 2.27$ (2D), $T_c \approx 4.5$ (3D)
- **淬火温度 (Quench Temperatures)**: 
  - $T = T_c/2 \approx 1.13$ (主要)
  - $T = 0.2 T_c, 0.7 T_c$ (对比实验)

## 预期结果 (Expected Results)

- 在 log-log 图中，$R(t)$ 应呈现线性关系
- 斜率应接近 1/3，验证 Lifshitz-Slyozov 标度律: $R \sim t^{1/3}$
- 空位介导动力学比传统自旋交换动力学更早进入标度区
- 2D 和 3D 系统应得到相同的增长指数

## 参考文献 (References)

- Lifshitz, I. M., & Slyozov, V. V. (1961). The kinetics of precipitation from supersaturated solid solutions. Journal of Physics and Chemistry of Solids, 19(1-2), 35-50.
- Ising Model and Monte Carlo Simulation
- Statistical Mechanics and Phase Transitions

## 许可证 (License)

This project is for academic and educational purposes.

## 贡献 (Contributing)

欢迎提出问题和改进建议！

Welcome to open issues and submit improvement suggestions!

## 联系方式 (Contact)

如有问题，请通过 GitHub Issues 联系。

For questions, please contact via GitHub Issues.
