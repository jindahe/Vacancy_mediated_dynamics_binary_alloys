import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.stats import linregress

def find_first_zero(r_values, c_values):
    """通过线性插值寻找 C(r) 的第一个零点"""
    for i in range(len(c_values) - 1):
        if c_values[i] * c_values[i+1] <= 0:
            # 线性插值公式: r0 = r1 - c1 * (r2 - r1) / (c2 - c1)
            r1, r2 = r_values[i], r_values[i+1]
            c1, c2 = c_values[i], c_values[i+1]
            return r1 - c1 * (r2 - r1) / (c2 - c1)
    return None

def analyze_simulation():
    # 1. 读取基于能量计算的 R (t_vs_R.txt)
    try:
        data_r = np.loadtxt('../output/t_vs_R.txt')
        times_e = data_r[:, 0]
        R_energy = data_r[:, 1]
    except FileNotFoundError:
        print("错误: 未找到 t_vs_R.txt")
        return

    # 2. 从 Cr_t_X.txt 文件中获取基于零点定义的 R
    cr_files = glob.glob('Cr_t_*.txt')
    # 按文件名中的数字排序
    cr_files.sort(key=lambda x: int(re.search(r't_(\d+)', x).group(1)))
    
    times_c = []
    R_zeros = []
    
    for f in cr_files:
        t = int(re.search(r't_(\d+)', f).group(1))
        if t == 0: continue # 初始随机状态通常没有零点
        
        data_c = np.loadtxt(f)
        r = data_c[:, 0]
        c = data_c[:, 1]
        
        z0 = find_first_zero(r, c)
        if z0:
            times_c.append(t)
            R_zeros.append(z0)

    # 3. 绘图分析
    plt.figure(figsize=(10, 6))
    
    # 绘制能量定义的 R
    plt.loglog(times_e[1:], R_energy[1:], 'o-', label='R from Energy (2/(<E>/N + 2))', alpha=0.7)
    
    # 绘制零点定义的 R
    if R_zeros:
        plt.loglog(times_c, R_zeros, 's--', label='R from C(r) Zero Crossing', alpha=0.7)

    # 4. 计算斜率 (针对后期数据进行拟合)
    # 选择 t > 100 的数据点进行线性回归
    fit_mask = times_e > 100
    if np.any(fit_mask):
        log_t = np.log10(times_e[fit_mask])
        log_r = np.log10(R_energy[fit_mask])
        slope, intercept, r_value, _, _ = linregress(log_t, log_r)
        
        # 绘制拟合线
        fit_line = 10**(slope * np.log10(times_e[fit_mask]) + intercept)
        plt.loglog(times_e[fit_mask], fit_line, 'k:', linewidth=2, 
                   label=f'Linear Fit (Slope={slope:.3f}, target=0.333)')
        print(f"拟合得到的幂律指数 (Scaling Exponent): {slope:.4f}")
        print(f"相关系数 (R-squared): {r_value**2:.4f}")

    plt.xlabel('Time (MCS)', fontsize=12)
    plt.ylabel('Domain Size R', fontsize=12)
    plt.title('Domain Growth Kinetics in Binary Alloys (Log-Log)', fontsize=14)
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    plt.savefig('../output/scaling_analysis.png')
    plt.show()

if __name__ == "__main__":
    analyze_simulation()