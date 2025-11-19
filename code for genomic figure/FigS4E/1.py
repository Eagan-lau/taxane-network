import pandas as pd
import matplotlib.pyplot as plt


file_path = 'TAXUS_PDB.txt'  
values = []


with open(file_path, 'r') as file:
    for line in file:
        try:
            value = float(line.strip())
            values.append(value)
        except ValueError:
            pass  


df = pd.DataFrame(values, columns=['Value'])


total_count = len(df)



reference_number = 61375
percentage_of_total = (total_count / reference_number) * 100



bins = [50, 60, 70, 80, 90, 100]


binned_counts = pd.cut(df['Value'], bins=bins, right=False, include_lowest=True).value_counts(sort=False)


binned_percentage = (binned_counts / reference_number) * 100



for bin_range, count, perc in zip(binned_counts.index, binned_counts, binned_percentage):
    print(f"{bin_range}: {count} ({perc:.2f}%)")


fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [1, 0.6]})


axes[0].pie([percentage_of_total, 100 - percentage_of_total], labels=[f'Values ({percentage_of_total:.2f}%)', 'Other'], autopct='%1.1f%%', startangle=90, colors=['skyblue', 'lightgrey'], counterclock=False)
axes[0].set_title('Percentage of Values in TAXUS_PDB.txt relative to 61601', fontsize=16)


binned_counts.plot(kind='bar', color='skyblue', edgecolor='black', ax=axes[1])
for i, (count, perc) in enumerate(zip(binned_counts, binned_percentage)):
    axes[1].text(i, count + 5, f'{perc:.2f}%', ha='center', va='bottom', fontsize=10)

axes[1].set_title('Distribution of Values by Range in TAXUS_PDB.txt', fontsize=16)
axes[1].set_xlabel('Value Range', fontsize=14)
axes[1].set_ylabel('Count', fontsize=14)
axes[1].set_xticklabels(binned_counts.index, rotation=0, fontsize=12)  


axes[1].tick_params(axis='y', labelsize=12)


max_count = binned_counts.max()
axes[1].set_ylim(0, max_count + max_count * 0.2) 


plt.tight_layout()


plt.show()
