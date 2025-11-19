import numpy as np
from scipy.stats import entropy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Function to read data from Excel sheet
def read_spectrum_from_excel(file_path, sheet_name):
    df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)
    mz_values = df[0].values
    intensities = df[1].values
    data = [(mz, intensity) for mz, intensity in zip(mz_values, intensities)]
    return data

# Preprocessing function
def preprocess_spectrum(data, mz_tolerance=0.01):
    data = [(mz, intensity) for mz, intensity in data if intensity > 0]
    if not data:
        return np.array([]), np.array([])

    mz_values = np.array([mz for mz, _ in data])
    intensities = np.array([intensity for _, intensity in data])

    if intensities.sum() > 0:
        intensities = intensities / intensities.sum()

    return mz_values, intensities

# Align spectra within m/z tolerance
def align_spectra(mz1, int1, mz2, int2, mz_tolerance=0.01):
    aligned_int1 = []
    aligned_int2 = []
    all_mz = sorted(set(np.concatenate([mz1, mz2])))

    for mz in all_mz:
        int1_val = 0
        int2_val = 0
        for i, m in enumerate(mz1):
            if abs(m - mz) <= mz_tolerance:
                int1_val = int1[i]
                break
        for i, m in enumerate(mz2):
            if abs(m - mz) <= mz_tolerance:
                int2_val = int2[i]
                break
        aligned_int1.append(int1_val)
        aligned_int2.append(int2_val)

    return np.array(aligned_int1), np.array(aligned_int2)

# Compute spectral entropy
def spectral_entropy(intensities):
    if intensities.sum() == 0:
        return 0
    intensities = intensities + 1e-10
    intensities = intensities / intensities.sum()
    return entropy(intensities, base=2)

# Compute spectral entropy similarity
def spectral_entropy_similarity(int1, int2):
    int1 = np.maximum(int1, 0) + 1e-10
    int2 = np.maximum(int2, 0) + 1e-10
    int1 = int1 / int1.sum()
    int2 = int2 / int2.sum()

    entropy1 = spectral_entropy(int1)
    entropy2 = spectral_entropy(int2)

    combined = (int1 + int2) / 2
    entropy_combined = spectral_entropy(combined)

    similarity = 1 - (entropy_combined - (entropy1 + entropy2) / 2) / np.log2(len(int1))
    return max(0, min(1, similarity))

# Main processing
def main():
    # File path to the Excel file
    file_path = '243-16672.xlsx'

    # Read all sheets
    xls = pd.read_excel(file_path, sheet_name=None)
    sheet_names = list(xls.keys())
    spectra_data = {name: read_spectrum_from_excel(file_path, name) for name in sheet_names}

    # Initialize similarity matrix
    n = len(sheet_names)
    similarity_matrix = np.zeros((n, n))

    # Compute pairwise similarities
    for i in range(n):
        for j in range(i, n):  # Symmetric matrix, compute upper triangle
            mz1, int1 = preprocess_spectrum(spectra_data[sheet_names[i]])
            mz2, int2 = preprocess_spectrum(spectra_data[sheet_names[j]])
            aligned_int1, aligned_int2 = align_spectra(mz1, int1, mz2, int2)
            sim = spectral_entropy_similarity(aligned_int1, aligned_int2)
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim  # Mirror to lower triangle

    # Create DataFrame for similarity matrix
    sim_df = pd.DataFrame(similarity_matrix, index=sheet_names, columns=sheet_names)

    # Print similarity matrix
    print("Spectral Entropy Similarity Matrix:")
    print(sim_df.round(3))

    # --- Clustered heatmap with soft, high-impact journal style ---
    sns.set_theme(style="white")

    g = sns.clustermap(
        sim_df,
        method='average',
        metric='euclidean',
        cmap='crest',        
        annot=True,
        fmt='.2f',
        annot_kws={"size": 10, "weight": 'bold'},
        figsize=(9, 9),
        dendrogram_ratio=0.1,
        cbar_pos=(0.02, 0.8, 0.03, 0.18),
        linewidths=0,
        square=True,
        xticklabels=True,
        yticklabels=True
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=40, ha='right', fontsize=11, weight='bold')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=11, weight='bold')
    g.ax_heatmap.set_title(
        "Clustered Heatmap of Spectral Entropy Similarities",
        fontsize=16, fontweight='bold', pad=14
    )
    g.cax.set_ylabel("Similarity", fontsize=12, fontweight='bold')

    
    # plt.savefig('spectral_similarity_clustermap.tif', dpi=600, bbox_inches='tight')

    plt.show()

if __name__ == "__main__":
    main()
