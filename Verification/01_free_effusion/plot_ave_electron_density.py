#!/usr/bin/env python3
"""
计算 dsmc_flowfield_*.vtk 文件中平均密度和热能量的时间演化。

处理 Pantera 输出的二进制 VTK (UNSTRUCTURED_GRID, big-endian) 文件，
提取 species 场的密度和温度，计算面积加权平均密度和总热能量。

绘制 4 幅图: 平均电子密度 / 平均离子密度 / 电子总能量 / 离子总能量

使用方法:
    python plot_ave_electron_density.py

输出:
    - average_electron_density_v2.csv  (时间演化数据)
    - average_electron_density_v2.png  (演化曲线图，如果 matplotlib 可用)
"""

from pathlib import Path
import re
import csv
import sys

import numpy as np

# ===================== 用户参数 =====================

# VTK 文件所在目录
data_dir = Path(r"F:\Code\Results\dumps4\dumps")

# VTK 文件名模式
file_pattern = "dsmc_flowfield_*.vtk"

# Species 配置 (从 species 文件和 input 中得来)
species_list = ["e-", "H+"]  # 第一个是电子，第二个是离子

# 物理时间步长 (s) — 与 input 文件中 Timestep 一致
timestep = 2.5e-9

# Boltzmann 常数 (J/K)
kB = 1.380649e-23

# ===================================================

# 尝试导入 matplotlib
try:
    import matplotlib.pyplot as plt

    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False
    print("[信息] matplotlib 未安装，仅保存 CSV 数据文件。")


def natural_sort_key(path):
    """按文件名中的数字排序 (0, 100, 200, ...)"""
    nums = re.findall(r"\d+", path.stem)
    return int(nums[-1]) if nums else 0


def read_vtk_fields(vtk_path, field_names):
    """
    从二进制 VTK 文件中一次性读取多个场数据。

    参数:
        vtk_path: VTK 文件路径
        field_names: 需要提取的字段名列表

    返回:
        (fields, areas, time)
        fields: dict {name: ndarray(n_cells,)}
        areas: ndarray(n_cells,) — 单元面积
        time: float — 物理时间
    """
    with open(vtk_path, "rb") as f:
        raw = f.read()

    # ---------- 辅助函数 ----------
    def find_str(pat, start=0):
        return raw.find(pat.encode("ascii"), start)

    def parse_line(pos):
        nl = raw.find(b"\n", pos)
        return raw[pos:nl].decode("ascii").strip(), nl

    # ---------- 读取网格 ----------
    # POINTS
    p = find_str("POINTS ")
    hdr, nl = parse_line(p)
    n_pts = int(hdr.split()[1])
    pts_bin = nl + 1
    pts_end = pts_bin + n_pts * 3 * 8
    nodes = np.frombuffer(raw[pts_bin:pts_end], dtype=">f8").reshape(-1, 3).copy()

    # CELLS
    c = find_str("CELLS ", nl)
    hdr, nl = parse_line(c)
    n_cells = int(hdr.split()[1])
    n_int = int(hdr.split()[2])
    cells_bin = nl + 1
    cells_end = cells_bin + n_int * 4
    cell_ints = (
        np.frombuffer(raw[cells_bin:cells_end], dtype=">i4")
        .reshape(-1, 4)
        .copy()
    )

    # CELL_TYPES
    ct = find_str("CELL_TYPES ", nl)
    _, nl = parse_line(ct)

    # CELL_DATA
    cd = find_str("CELL_DATA ", nl)
    _, nl = parse_line(cd)

    # FIELD FieldData count
    fd = find_str("FIELD FieldData ", nl)
    hdr, nl = parse_line(fd)
    n_arrays = int(hdr.split()[2])

    # ---------- 扫描所有数组，收集需要的 ----------
    pos = nl + 1
    found = {}
    needed = set(field_names)

    for _ in range(n_arrays):
        nl = raw.find(b"\n", pos)
        hdr = raw[pos:nl].decode("ascii").strip()
        parts = hdr.split()
        name = parts[0]
        ndim = int(parts[1])
        ntuples = int(parts[2])
        dtype_str = parts[3]

        elem_size = 8 if dtype_str in ("double",) else 4
        data_start = nl + 1
        data_len = ndim * ntuples * elem_size

        if name in needed:
            dtype = ">f8" if dtype_str in ("double",) else ">i4"
            arr = np.frombuffer(
                raw[data_start : data_start + data_len], dtype=dtype
            ).copy()
            found[name] = arr
            needed.remove(name)

        pos = data_start + data_len + 1

    if needed:
        raise KeyError(
            f"{vtk_path.name} 中未找到字段: {needed}"
        )

    # ---------- 计算单元面积 (2D 三角形) ----------
    v1 = nodes[cell_ints[:, 2]] - nodes[cell_ints[:, 1]]
    v2 = nodes[cell_ints[:, 3]] - nodes[cell_ints[:, 1]]
    cross_z = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    areas = 0.5 * np.abs(cross_z)

    # ---------- 时间 ----------
    nums = re.findall(r"\d+", vtk_path.stem)
    time_val = (int(nums[-1]) if nums else 0) * timestep

    return found, areas, time_val


def area_weighted_mean(data, areas):
    """面积加权平均。"""
    valid = (areas > 0) & np.isfinite(data)
    if not valid.any():
        return 0.0
    return np.sum(data[valid] * areas[valid]) / np.sum(areas[valid])


def density_weighted_mean(data, density, areas):
    """密度加权平均 (用于温度)。"""
    valid = (areas > 0) & np.isfinite(data) & np.isfinite(density) & (density > 0)
    if not valid.any():
        return 0.0
    weights = density[valid] * areas[valid]
    return np.sum(data[valid] * weights) / np.sum(weights)


# ===================== 场名推导 =====================

# 定义需要读取的字段
fields_to_read = []
density_names = {}  # species -> density field name
temp_names = {}  # species -> temperature field name

for sp in species_list:
    dn = f"nrho_mean_{sp}"
    tn = f"Ttr_mean_{sp}"
    fields_to_read.append(dn)
    fields_to_read.append(tn)
    density_names[sp] = dn
    temp_names[sp] = tn

# ===================== 主流程 =====================

vtk_files = sorted(data_dir.glob(file_pattern), key=natural_sort_key)

if not vtk_files:
    print(f"错误: 在 {data_dir} 中未找到 {file_pattern}")
    sys.exit(1)

# CSV 列头
csv_fields = [
    "file",
    "timestep",
    "time_s",
]
for sp in species_list:
    csv_fields += [
        f"avg_density_{sp}_m-3",
        f"density_wavg_temp_{sp}_K",
        f"total_thermal_energy_{sp}_J",
    ]

print(f"找到 {len(vtk_files)} 个 VTK 文件，species: {species_list}")
print(f"需要读取的场: {fields_to_read}")
print(f"{'=' * 90}")

# 表格头
header = (
    f"{'#':>4s}  "
    f"{'t(us)':>10s}  "
    f"{'<ne>(m-3)':>14s}  "
    f"{'<Te>(K)':>12s}  "
    f"{'Ee(J)':>14s}  "
    f"{'<ni>(m-3)':>14s}  "
    f"{'<Ti>(K)':>12s}  "
    f"{'Ei(J)':>14s}"
)
print(header)
print("-" * 90)

results = []

for i, f in enumerate(vtk_files):
    fields, areas, t = read_vtk_fields(f, fields_to_read)

    row = {
        "file": f.name,
        "timestep": int(re.findall(r"\d+", f.stem)[-1]),
        "time_s": t,
    }

    # 打印字符串
    line_parts = [f"{i+1:4d}/{len(vtk_files):4d}", f"{t*1e6:>8.2f}"]

    for sp in species_list:
        density = fields[density_names[sp]]
        temperature = fields[temp_names[sp]]

        # 面积加权平均密度
        avg_n = area_weighted_mean(density, areas)

        # 密度加权平均温度
        avg_T = density_weighted_mean(temperature, density, areas)

        # 总热能量 (所有同类粒子的能量之和) 单位 J
        # 每个单元中粒子数 = density * area, 每个粒子能量 = 1.5 * kB * T
        total_energy_J = 1.5 * kB * np.sum(density * temperature * areas)

        row[f"avg_density_{sp}_m-3"] = avg_n
        row[f"density_wavg_temp_{sp}_K"] = avg_T
        row[f"total_thermal_energy_{sp}_J"] = total_energy_J

        line_parts.append(f"{avg_n:>12.4e}")
        line_parts.append(f"{avg_T:>10.2f}")
        line_parts.append(f"{total_energy_J:>12.4e}")

    results.append(row)
    print("  ".join(line_parts))

print("=" * 90)

# ===================== 保存 CSV =====================

csv_file = data_dir / "average_electron_density_v2.csv"
with open(csv_file, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=csv_fields)
    w.writeheader()
    w.writerows(results)

print(f"\n数据已保存: {csv_file}")

# ===================== 绘图 =====================

if HAVE_MPL:
    times = np.array([r["time_s"] for r in results])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # ---- 第一幅: 平均电子密度 ----
    ax = axes[0, 0]
    ne = np.array([r["avg_density_e-_m-3"] for r in results])
    ax.plot(times * 1e6, ne / 1e11, linewidth=1.5, color="#1f77b4")
    ax.set_xlabel("Time (μs)")
    ax.set_ylabel(r"$\langle n_e \rangle$ ($\times 10^{11}$ m$^{-3}$)")
    ax.set_title("Average Electron Density")
    ax.grid(True, alpha=0.3)

    # ---- 第二幅: 平均离子密度 ----
    ax = axes[0, 1]
    ni = np.array([r["avg_density_H+_m-3"] for r in results])
    ax.plot(times * 1e6, ni / 1e11, linewidth=1.5, color="#d62728")
    ax.set_xlabel("Time (μs)")
    ax.set_ylabel(r"$\langle n_i \rangle$ ($\times 10^{11}$ m$^{-3}$)")
    ax.set_title("Average Ion Density")
    ax.grid(True, alpha=0.3)

    # ---- 第三幅: 电子总能量 ----
    ax = axes[1, 0]
    Ee = np.array([r["total_thermal_energy_e-_J"] for r in results])
    ax.plot(times * 1e6, Ee, linewidth=1.5, color="#1f77b4")
    ax.set_xlabel("Time (μs)")
    ax.set_ylabel("Total Energy (J)")
    ax.set_title("Total Electron Thermal Energy")
    ax.grid(True, alpha=0.3)

    # ---- 第四幅: 离子总能量 ----
    ax = axes[1, 1]
    Ei = np.array([r["total_thermal_energy_H+_J"] for r in results])
    ax.plot(times * 1e6, Ei, linewidth=1.5, color="#d62728")
    ax.set_xlabel("Time (μs)")
    ax.set_ylabel("Total Energy (J)")
    ax.set_title("Total Ion Thermal Energy")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Average Electron / Ion Density and Total Energy Evolution",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    png_file = data_dir / "average_electron_density_v2.png"
    fig.savefig(png_file, dpi=150)
    print(f"图片已保存: {png_file}")

    plt.show()

print("完成!")
