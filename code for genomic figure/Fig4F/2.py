import pandas as pd
import matplotlib.pyplot as plt


combined_file = 'combined_b_factors.txt'


b_factors = []


with open(combined_file, 'r') as file:
    for line in file:
        b_factor = float(line.strip())
        b_factors.append(b_factor)


df = pd.DataFrame(b_factors, columns=['B-factor'])


bins = [0, 50, 70, 90, 100]


binned_b_factors = pd.cut(df['B-factor'], bins=bins, right=False, include_lowest=True)
binned_counts = binned_b_factors.value_counts(sort=False)


binned_percentage = binned_counts / binned_counts.sum() * 100


plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)  
binned_counts.plot(kind='bar', color='skyblue', edgecolor='black')
plt.title('Residue Count by B-factor Range')
plt.xlabel('B-factor Range')
plt.ylabel('Residue Count')
plt.xticks(rotation=0)  


plt.subplot(1, 2, 2)  
binned_percentage.plot.pie(autopct='%1.1f%%', startangle=90, counterclock=False, labels=binned_percentage.index)
plt.title('B-factor Distribution')
plt.ylabel('')  


plt.tight_layout()
plt.show()
