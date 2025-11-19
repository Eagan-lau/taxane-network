import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

input_file = 'pLDDT_percentages2.xlsx'
sheets_dict = pd.read_excel(input_file, sheet_name=None)
pLDDT_ranges = ['0-50', '50-70', '70-90', '90-100']

group_names = list(sheets_dict.keys())
x = np.arange(len(group_names))  


all_ranges = []
for r in pLDDT_ranges:
    vals = []
    for sheet_name in group_names:
        df = sheets_dict[sheet_name]
        percent = df.set_index('pLDDT_range').reindex(pLDDT_ranges).loc[r, 'percent']
        vals.append(percent)
    all_ranges.append(vals)

width = 0.18  
plt.figure(figsize=(max(10, len(group_names)*0.6), 6))


for i, r in enumerate(pLDDT_ranges):
    plt.bar(x + i*width - 1.5*width, all_ranges[i], width, label=r)


high_conf_sum = np.array(all_ranges[2]) + np.array(all_ranges[3])
plt.plot(x, high_conf_sum, color='red', marker='o', linewidth=2, label='70-100% sum')


plt.xticks(x, group_names, rotation=60, ha='right')
plt.xlabel('Group (Sheet name)')
plt.ylabel('Percent')
plt.ylim(0, 100)
plt.title('pLDDT Range Distribution per Group\n(bar: all ranges; line: high-confidence sum)')
plt.legend(title='pLDDT range', loc='upper left')
plt.tight_layout()
plt.savefig('pLDDT_per_group_bar+line.png', dpi=300)
plt.show()
