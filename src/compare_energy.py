import numpy as np
import matplotlib.pyplot as plt

def get_energy_from_R(R_values, J=1.0):
    """根据公式 R = 2 / (E_avg + 2) 反推平均能量 E_avg"""
    # 逆变换: E_avg = (2 / R) - 2
    # 注意：这里的 E_avg 是无量纲化的 E/J
    # 避免 R=0 的除零错误
    safe_R = np.where(R_values == 0, 1e-6, R_values)
    return (2.0 / safe_R) - 2.0

def plot_energy_comparison():
    try:
        data_vac = np.loadtxt('../output/t_vs_R.txt')
        data_kaw = np.loadtxt('../output/t_vs_R_kawasaki.txt')
    except IOError:
        return

    t_vac, E_vac = data_vac[:, 0], get_energy_from_R(data_vac[:, 1])
    t_kaw, E_kaw = data_kaw[:, 0], get_energy_from_R(data_kaw[:, 1])

    plt.figure(figsize=(10, 6))
    
    plt.plot(t_vac, E_vac, 'b-', label='Vacancy Mediated')
    plt.plot(t_kaw, E_kaw, 'r--', label='Kawasaki Dynamics')

    plt.xscale('log') # 使用半对数坐标，更能看清早期的能量骤降
    plt.xlabel(r'Time $t$ (MCS) [Log Scale]', fontsize=14)
    plt.ylabel(r'Average Energy per Site $\langle E \rangle / N$', fontsize=14)
    plt.title(r'Energy Relaxation Comparison', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, which="both", alpha=0.3)
    
    # 理论基态能量 (对于完全分离的 Ising 模型，内部全为 -J，只有界面有能量)
    # 在 T < Tc 时，能量应趋向于 -2 (每个自旋有4个邻居，E = -4J/2 = -2J)
    plt.axhline(-2.0, color='k', linestyle=':', label='Ground State Limit (-2J)')

    plt.savefig('../output/compare_3_energy.png')
    print("已生成: ../output/compare_3_energy.png")
    plt.show()

if __name__ == "__main__":
    plot_energy_comparison()