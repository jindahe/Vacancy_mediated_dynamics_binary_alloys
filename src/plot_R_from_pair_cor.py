import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.stats import linregress

def find_zero_crossing(r, c):
    for i in range(len(c)-1):
        if c[i] * c[i+1] <= 0:
            return r[i] - c[i] * (r[i+1] - r[i]) / (c[i+1] - c[i])
    return None

def plot_scaling_comparison():
    # 获取能量数据
    e_data = np.loadtxt('../output/t_vs_R.txt')
    t_e, R_e = e_data[:, 0], e_data[:, 1]

    # 获取关联函数零点数据
    cr_files = glob.glob('../output/Cr_t_*.txt')
    cr_files.sort(key=lambda x: int(re.search(r't_(\d+)', x).group(1)))
    
    t_c, R_c = [], []
    for f in cr_files:
        t_val = int(re.search(r't_(\d+)', f).group(1))
        data = np.loadtxt(f)
        z = find_zero_crossing(data[:, 0], data[:, 1])
        if z:
            t_c.append(t_val)
            R_c.append(z)

    plt.figure(figsize=(10, 6))
    plt.loglog(t_e[1:], R_e[1:], 'o-', label='R (Energy-based)')
    plt.loglog(t_c, R_c, 's--', label='R (C(r) zero-crossing)')

    # 线性拟合验证 1/3 指数
    if len(t_e) > 5:
        log_t = np.log10(t_e[t_e > 100])
        log_r = np.log10(R_e[t_e > 100])
        slope, intercept, _, _, _ = linregress(log_t, log_r)
        plt.loglog(t_e[t_e > 100], 10**(slope * log_t + intercept), 'k--', 
                   label=f'Fit slope: {slope:.3f} (Theory: 0.333)')

    plt.xlabel('log(t)')
    plt.ylabel('log(R)')
    plt.title(r'Scaling Law Verification $R \sim t^{1/3}$')
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.savefig('../output/final_scaling_verification.png')
    plt.show()

if __name__ == "__main__":
    plot_scaling_comparison()