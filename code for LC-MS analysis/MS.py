import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# ========== 手动设置显示的保留时间区间 ==========
start_time = 14.8   # 设置起始时间
end_time = 19.2     # 设置结束时间
# =============================================

# Excel 文件路径
file_path = '11.1-9-di-13-ac多个酶.xlsx'
xls = pd.ExcelFile(file_path)
sheet_names = xls.sheet_names

# 收集有效 sheet
valid_sheets = []
for sheet in sheet_names:
    df = pd.read_excel(xls, sheet_name=sheet)
    if 'Time' in df.columns and 'Intensity' in df.columns and len(df) >= 11:
        valid_sheets.append(sheet)

# 创建单一图
fig, ax = plt.subplots(figsize=(12, 8))

offset_increment = 300
colors = plt.cm.tab10(np.linspace(0, 1, len(valid_sheets)))
intensity_threshold = 1.0e5

for idx, sheet in enumerate(valid_sheets):
    df = pd.read_excel(xls, sheet_name=sheet)
    # ========== 仅保留手动设定的时间区间的数据 ==========
    df = df[(df['Time'] >= start_time) & (df['Time'] <= end_time)]
    if df.empty:  # 没有数据则跳过
        continue
    # =============================================
    x = df['Time']
    y = df['Intensity']

    y_smooth = savgol_filter(y, window_length=11, polyorder=3)
    y_smooth[y_smooth < 0] = 0

    y_filtered = y_smooth.copy()
    y_filtered[y_filtered < intensity_threshold] = 0

    max_intensity = np.max(y_filtered) if np.max(y_filtered) > 0 else 1
    y_norm = y_filtered / max_intensity * 200
    y_offset = y_norm + (idx * offset_increment)

    ax.plot(x, y_offset, label=sheet, color=colors[idx], linewidth=2)
    ax.text(x.max() + 0.5, idx * offset_increment + 50, sheet, 
            verticalalignment='center', fontsize=10, color=colors[idx])

# 设置X轴
ax.set_xlabel('Retention Time (min)', fontsize=12)
ax.set_xlim(start_time, end_time + 1)

ax.set_yticks([])
ax.set_ylabel('Intensity (a.u.)', fontsize=10)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('11.1-9-di-13-ac多个酶.pdf', dpi=300)
plt.show()
