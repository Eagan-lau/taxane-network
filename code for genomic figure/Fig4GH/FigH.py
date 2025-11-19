import matplotlib.pyplot as plt

# 1. 读取hap1_catalytic_activity.txt和hap2_catalytic_activity.txt的内容，并切割和合并成一个列表
def read_and_split(file_name, delimiter):
    with open(file_name, 'r') as file:
        content = file.read().strip()
        return content.split(delimiter)

hap1_list = read_and_split('hap1_catalytic_activity.txt', ';')
hap2_list = read_and_split('hap2_catalytic_activity.txt', ';')

list1 = hap1_list + hap2_list

# 2. 读取Q1_protein_ids.txt, Q2_protein_ids.txt, Q3_protein_ids.txt, Q4_protein_ids.txt，逐行读取，并根据'.'进行切割，取第一个元素
def read_and_split_by_dot(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()
        return [line.strip().split('.')[0] for line in lines]

Q1_list = read_and_split_by_dot('Q1_protein_ids.txt')
Q2_list = read_and_split_by_dot('Q2_protein_ids.txt')
Q3_list = read_and_split_by_dot('Q3_protein_ids.txt')
Q4_list = read_and_split_by_dot('Q4_protein_ids.txt')

# 3. 判断list1中的元素在哪个列表中，并进行标记
def classify_element(element, q1, q2, q3, q4):
    if element in q1:
        return 'Q1'
    elif element in q2:
        return 'Q2'
    elif element in q3:
        return 'Q3'
    elif element in q4:
        return 'Q4'
    else:
        return 'None'

classification = [classify_element(elem, Q1_list, Q2_list, Q3_list, Q4_list) for elem in list1]

# 4. 统计list1中元素在四个列表中的分布频率
distribution = {
    'Q1': classification.count('Q1'),
    'Q2': classification.count('Q2'),
    'Q3': classification.count('Q3'),
    'Q4': classification.count('Q4'),
    'None': classification.count('None')
}

print("Distribution in each list:")
for key, value in distribution.items():
    print(f"{key}: {value}")

# 计算百分比
total_elements = len(list1)
percent_distribution = {key: (value / total_elements) * 100 for key, value in distribution.items() if key != 'None'}

print("\nPercentage Distribution:")
for key, value in percent_distribution.items():
    print(f"{key}: {value:.2f}%")

# 进行百分比的绘图
plt.figure(figsize=(10, 6))
labels = list(percent_distribution.keys())
sizes = list(percent_distribution.values())
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']
explode = (0.1, 0, 0, 0)  # 仅 "Q1" 突出显示

plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
plt.title('Percentage Distribution of Elements in List1 among Q1, Q2, Q3, Q4')
plt.show()
