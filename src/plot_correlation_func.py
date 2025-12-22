import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def plot_correlations():
    files = glob.glob('../output/Cr_t_*.txt')
    files.sort(key=lambda x: int(re.search(r't_(\d+)', x).group(1)))

    plt.figure(figsize=(8, 6))
    for f in files:
        t = re.search(r't_(\d+)', f).group(1)
        data = np.loadtxt(f)
        plt.plot(data[:, 0], data[:, 1], label=f't = {t}')

    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel('Distance r')
    plt.ylabel('C(r)')
    plt.title('Pair Correlation Function $C(r)$ over Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('../output/correlation_functions.png')
    plt.show()

if __name__ == "__main__":
    plot_correlations()