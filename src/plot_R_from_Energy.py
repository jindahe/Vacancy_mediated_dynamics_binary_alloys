import numpy as np
import matplotlib.pyplot as plt

def plot_r_energy():
    try:
        data = np.loadtxt('../output/t_vs_R.txt')
    except:
        print("未找到 output/t_vs_R.txt")
        return

    t = data[:, 0]
    R = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.loglog(t[1:], R[1:], 'bo-', label='R from Energy')
    plt.xlabel('Time t (MCS)')
    plt.ylabel('Domain Size R')
    plt.title('Domain Growth $R(t)$ (Energy-based)')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    plt.savefig('../output/R_energy_scaling.png')
    plt.show()

if __name__ == "__main__":
    plot_r_energy()