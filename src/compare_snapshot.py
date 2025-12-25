import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def plot_snapshot_comparison(target_t=None):
    # 如果不指定时间，默认找 output 文件夹里最大的那个时间点
    if target_t is None:
        files = glob.glob('../output/lattice_t_*.txt')
        if not files:
            print("未找到晶格文件。")
            return
        # 提取文件名中的数字并找到最大值
        times = [int(re.search(r't_(\d+)', f).group(1)) for f in files]
        target_t = max(times)

    # 构建文件名
    file_vac = f'../output/lattice_t_{target_t}.txt'
    file_kaw = f'../output/lattice_kawasaki_t_{target_t}.txt'

    try:
        # 读取数据
        lat_vac = np.loadtxt(file_vac)
        lat_kaw = np.loadtxt(file_kaw)
    except IOError:
        print(f"错误：找不到时刻 t={target_t} 的数据文件。请确保两个模拟都运行到了该时刻。")
        return

    # 绘图
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # 空位动力学快照
    axes[0].imshow(lat_vac, cmap='coolwarm', interpolation='nearest')
    axes[0].set_title(f'Vacancy Dynamics\n(t = {target_t} MCS)', fontsize=14)
    axes[0].axis('off')

    # Kawasaki 动力学快照
    axes[1].imshow(lat_kaw, cmap='coolwarm', interpolation='nearest')
    axes[1].set_title(f'Kawasaki Dynamics\n(t = {target_t} MCS)', fontsize=14)
    axes[1].axis('off')

    plt.suptitle(f"Microstructure Comparison at t={target_t}", fontsize=16)
    plt.tight_layout()
    plt.savefig(f'../output/compare_2_snapshots_t{target_t}.png')
    print(f"已生成: ../output/compare_2_snapshots_t{target_t}.png")
    plt.show()

if __name__ == "__main__":
    # 你可以手动指定一个存在的时间点，例如 4096 或 8192
    plot_snapshot_comparison(4096)