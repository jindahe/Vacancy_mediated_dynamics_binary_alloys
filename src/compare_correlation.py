import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def plot_correlation_comparison(target_t=None):
    # 自动寻找最大的共有时间点
    if target_t is None:
        files_vac = glob.glob('../output/Cr_t_*.txt')
        times_vac = set(int(re.search(r't_(\d+)', f).group(1)) for f in files_vac)
        
        files_kaw = glob.glob('../output/Cr_kawasaki_t_*.txt')
        times_kaw = set(int(re.search(r't_(\d+)', f).group(1)) for f in files_kaw)
        
        common_times = times_vac.intersection(times_kaw)
        if not common_times:
            print("错误：未找到两个模拟共有的时间点文件。")
            return
        target_t = max(common_times)

    file_vac = f'../output/Cr_t_{target_t}.txt'
    file_kaw = f'../output/Cr_kawasaki_t_{target_t}.txt'

    data_vac = np.loadtxt(file_vac)
    data_kaw = np.loadtxt(file_kaw)

    plt.figure(figsize=(10, 6))
    
    plt.plot(data_vac[:, 0], data_vac[:, 1], 'b-o', markersize=4, label=f'Vacancy (t={target_t})')
    plt.plot(data_kaw[:, 0], data_kaw[:, 1], 'r--s', markersize=4, label=f'Kawasaki (t={target_t})')

    plt.axhline(0, color='k', linestyle='-', alpha=0.5)
    
    # 标注第一个零点 (即畴尺寸 R)
    # 简单寻找符号改变的位置
    def find_zero(r, c):
        for i in range(len(c)-1):
            if c[i] * c[i+1] <= 0: return r[i]
        return 0
    
    zero_vac = find_zero(data_vac[:, 0], data_vac[:, 1])
    zero_kaw = find_zero(data_kaw[:, 0], data_kaw[:, 1])
    
    plt.axvline(zero_vac, color='b', linestyle=':', alpha=0.5, label=f'R_vac ~ {zero_vac}')
    plt.axvline(zero_kaw, color='r', linestyle=':', alpha=0.5, label=f'R_kaw ~ {zero_kaw}')

    plt.xlabel(r'Distance $r$', fontsize=14)
    plt.ylabel(r'Pair Correlation $C(r)$', fontsize=14)
    plt.title(f'Spatial Correlation Comparison at t={target_t}', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)

    plt.savefig('../output/compare_4_correlation.png')
    print("已生成: ../output/compare_4_correlation.png")
    plt.show()

if __name__ == "__main__":
    # 输入target_t参数以指定时间点，或留空以自动寻找最大共有时间点
    plot_correlation_comparison()