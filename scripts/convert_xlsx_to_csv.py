import pandas as pd
import sys
import os

script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))  # Get the script's directory
sys.path.append(script_dir)  # Add Visualization/ to sys.path

import GOATpy as gp

fd = gp.dh.GetFolder()

# Load the Excel file
excel_file = fd+"/R.xlsx"  # Change this to your actual file name

# Read the first sheet (or specify sheet_name="Sheet1" if needed)
df = pd.read_excel(excel_file, sheet_name=0, dtype=str)  # Read as strings to avoid formatting issues

# Save to CSV with proper comma separation
csv_file = fd+"/R.csv"  # Change output file name as needed
df.to_csv(csv_file, index=False, sep=",", encoding="utf-8")

print(f"Converted '{excel_file}' to '{csv_file}' with proper comma separation.")

# Load the Excel file
excel_file = fd+"/Z.xlsx"  # Change this to your actual file name

# Read the first sheet (or specify sheet_name="Sheet1" if needed)
df = pd.read_excel(excel_file, sheet_name=0, dtype=str)  # Read as strings to avoid formatting issues

# Save to CSV with proper comma separation
csv_file = fd+"/Z.csv"  # Change output file name as needed
df.to_csv(csv_file, index=False, sep=",", encoding="utf-8")

print(f"Converted '{excel_file}' to '{csv_file}' with proper comma separation.")

# Load the Excel file
excel_file = fd+"/psi.xlsx"  # Change this to your actual file name

# Read the first sheet (or specify sheet_name="Sheet1" if needed)
df = pd.read_excel(excel_file, sheet_name=0, dtype=str)  # Read as strings to avoid formatting issues

# Save to CSV with proper comma separation
csv_file = fd+"/psi.csv"  # Change output file name as needed
df.to_csv(csv_file, index=False, sep=",", encoding="utf-8")

print(f"Converted '{excel_file}' to '{csv_file}' with proper comma separation.")
