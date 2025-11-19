# -*- coding: utf-8 -*-


import pandas as pd


INPUT_FILE = "pLDDT_percentages2.xlsx"
OUTPUT_FILE = "merged_pLDDT.xlsx"

def main():
    
    xls = pd.ExcelFile(INPUT_FILE)
    sheet_names = xls.sheet_names
    

    
    all_data = []
    for sheet in sheet_names:
        df = pd.read_excel(INPUT_FILE, sheet_name=sheet)
        df["Species"] = sheet  
        all_data.append(df)

    merged_df = pd.concat(all_data, ignore_index=True)

    
    merged_df.to_excel(OUTPUT_FILE, index=False)
    

if __name__ == "__main__":
    main()
