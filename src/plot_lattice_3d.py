import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import glob
import re
import os
import sys

# ===========================
# 参数设置
# ===========================
L = 64            # 必须与 C++ 代码中的 L 一致
INPUT_DIR = '../output_3d'
OUTPUT_GIF = '../output_3d/3d_evolution.gif'
FPS = 2           # 每秒帧数，数值越小动图越慢
SLICE_INDEX = L // 2  # 2D切片的位置（中间）
# ===========================

def generate_3d_gif():
    # 1. 获取并排序文件
    files = glob.glob(os.path.join(INPUT_DIR, 'lattice_3d_t_*.txt'))
    print(f"找到的文件列表: {files}")
    # 使用正则表达式提取时间步 t 进行数值排序
    files.sort(key=lambda x: int(re.search(r't_(\d+)\.txt', x).group(1)))

    if not files:
        print(f"错误: 在 {INPUT_DIR}/ 中未找到数据文件。请先运行 C++ 模拟。")
        sys.exit(1)

    print(f"找到 {len(files)} 个数据文件，开始生成动画...")

    # 2. 初始化画布 (只需设置一次)
    fig = plt.subplots(figsize=(12, 6))
    # 调整子图布局，留出标题空间
    fig[0].subplots_adjust(left=0.05, right=0.95, wspace=0.1, top=0.9)
    
    # ---- 初始化 2D 切片视图 (左图) ----
    ax2d = fig[0].add_subplot(121)
    # 先读取第一帧数据用于初始化图像对象
    first_data = np.loadtxt(files[0]).reshape((L, L, L))
    # 使用固定颜色范围 vmin/vmax，确保颜色映射在整个动画中保持一致
    im2d = ax2d.imshow(first_data[SLICE_INDEX, :, :], cmap='coolwarm', vmin=-1, vmax=1)
    ax2d.set_title(f"2D Mid-Slice (Z={SLICE_INDEX})")
    ax2d.axis('off')

    # ---- 初始化 3D 散点视图 (右图) ----
    ax3d = fig[0].add_subplot(122, projection='3d')
    ax3d.set_xlim(0, L)
    ax3d.set_ylim(0, L)
    ax3d.set_zlim(0, L)
    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlim(0, L)
    # 初始化一个空的散点图对象，之后在 update 中填充
    scatter_plot = ax3d.scatter([], [], [], c='red', s=2, alpha=0.15)

    # 全局标题
    main_title = fig[0].suptitle(f"3D Domain Evolution (t = 0 MCS)", fontsize=16)

    # 3. 定义动画更新函数
    def update(frame_idx):
        current_file = files[frame_idx]
        t_val = re.search(r't_(\d+)\.txt', current_file).group(1)
        print(f"Processing frame {frame_idx+1}/{len(files)}: t = {t_val} MCS...")

        # 读取并重构数据
        data_flat = np.loadtxt(current_file)
        lattice = data_flat.reshape((L, L, L))

        # ---- 更新 2D 视图 ----
        # 只需更新数据，无需重绘坐标轴
        im2d.set_array(lattice[SLICE_INDEX, :, :])

        # ---- 更新 3D 视图 ----
        # Matplotlib 3D 散点图更新比较复杂，最稳妥的方法是清空后重绘
        ax3d.clear()
        ax3d.set_xlim(0, L); ax3d.set_ylim(0, L); ax3d.set_zlim(0, L)
        ax3d.set_title("3D Morphology (A atoms only)")
        ax3d.set_xlabel('X'); ax3d.set_ylabel('Y'); ax3d.set_zlabel('Z')
        
        # 找到所有 A 原子 (值为 1) 的坐标
        x, y, z = np.where(lattice == 1)
        # 绘制新的散点
        ax3d.scatter(x, y, z, c='red', s=2, alpha=0.15)

        # 更新全局标题
        main_title.set_text(f"3D Domain Evolution (t = {t_val} MCS)")
        
        # 返回需要更新的艺术家对象 (对于 3D plot，这步不是严格必须的，但对 2D 是)
        return im2d, main_title

    # 4. 创建并保存动画
    # interval 计算公式：1000ms / FPS
    ani = animation.FuncAnimation(fig[0], update, frames=len(files), interval=1000/FPS, blit=False)
    
    print(f"正在保存 GIF 到 {OUTPUT_GIF} (这可能需要几秒钟)...")
    # 需要安装 pillow: pip install pillow
    ani.save(OUTPUT_GIF, writer='pillow')
    print("动画生成完成！")
    # plt.show() # 如果想在窗口中查看，取消注释

if __name__ == "__main__":
    # 检查依赖
    try:
        import PIL
    except ImportError:
        print("错误: 生成 GIF 需要安装 pillow 库。请运行: pip install pillow")
        sys.exit(1)
        
    generate_3d_gif()