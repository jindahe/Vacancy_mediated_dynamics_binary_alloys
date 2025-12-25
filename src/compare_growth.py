import numpy as np
import matplotlib.pyplot as plt

def plot_growth_comparison():
    # 1. 读取数据
    # 假设空位动力学文件为 t_vs_R.txt
    # 假设 Kawasaki 动力学文件为 t_vs_R_kawasaki.txt
    try:
        data_vac = np.loadtxt('../output/t_vs_R.txt')
        data_kaw = np.loadtxt('../output/t_vs_R_kawasaki.txt')
    except IOError:
        print("错误：无法找到数据文件，请检查路径。")
        return

    t_vac, R_vac = data_vac[:, 0], data_vac[:, 1]
    t_kaw, R_kaw = data_kaw[:, 0], data_kaw[:, 1]

    plt.figure(figsize=(10, 7))
    
    # 绘制空位动力学 (蓝色实线)
    plt.loglog(t_vac[1:], R_vac[1:], 'b-o', markersize=4, label='Vacancy Mediated (Optimized)')
    
    # 绘制 Kawasaki 动力学 (红色虚线)
    plt.loglog(t_kaw[1:], R_kaw[1:], 'r--s', markersize=4, label='Spin Exchange (Kawasaki)')

    # 添加理论标度线 t^(1/3)
    # 选取一段中间时间区域绘制参考线
    ref_t = np.logspace(2, 4, 100) # 100 到 10000 MCS
    ref_R = 0.6 * (ref_t ** (1/3)) # 系数 0.6 用于调整参考线高度
    plt.loglog(ref_t, ref_R, 'k:', linewidth=2, label=r'Theory Slope $1/3$')

    plt.xlabel(r'Time $t$ (MCS)', fontsize=14)
    plt.ylabel(r'Domain Size $R(t)$', fontsize=14)
    plt.title(r'Comparison of Domain Growth: Vacancy vs Kawasaki', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, which="both", ls="-", alpha=0.2)
    
    plt.savefig('../output/compare_1_growth.png')
    print("已生成: ../output/compare_1_growth.png")
    plt.show()

if __name__ == "__main__":
    plot_growth_comparison()