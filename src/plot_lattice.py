import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import re
import os

def generate_lattice_gif():
    # 1. 获取并按时间顺序排列所有晶格文件
    files = glob.glob('../output/lattice_t_*.txt')
    files.sort(key=lambda x: int(re.search(r't_(\d+)', x).group(1)))

    if not files:
        print("错误: 在 ../output/ 文件夹下未找到晶格数据文件。")
        return

    # 2. 设置绘图环境
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # 初始化第一帧
    initial_lattice = np.loadtxt(files[0])
    im = ax.imshow(initial_lattice, cmap='coolwarm', interpolation='nearest')
    title = ax.set_title(f"Lattice Evolution: t = {re.search(r't_(\d+)', files[0]).group(1)} MCS")
    ax.axis('off')

    def update(frame_file):
        # 读取新的一帧数据
        lattice = np.loadtxt(frame_file)
        t = re.search(r't_(\d+)', frame_file).group(1)
        
        # 更新图像内容
        im.set_array(lattice)
        title.set_text(f"Lattice Evolution: t = {t} MCS")
        return [im, title]

    # 3. 创建动画
    # frames 传入文件列表，interval 为帧间隔（毫秒）
    ani = animation.FuncAnimation(fig, update, frames=files, interval=500, blit=True)

    # 4. 保存为 GIF (需要安装 pillow 库)
    output_path = '../output/lattice_evolution.gif'
    try:
        ani.save(output_path, writer='pillow')
        print(f"成功生成动图: {output_path}")
    except Exception as e:
        print(f"保存失败: {e}")

    plt.show()

if __name__ == "__main__":
    generate_lattice_gif()