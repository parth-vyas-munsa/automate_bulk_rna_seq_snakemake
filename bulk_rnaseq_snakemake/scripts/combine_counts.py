import pandas as pd
import glob
import os

# 1. Find all count files
files = glob.glob("results/counts/*.txt")

# 2. Read the first file
df = pd.read_csv(files[0], sep="\t", comment="#")

# 3. Keep Geneid and the count column
sample = df[["Geneid", df.columns[-1]]]

# 4. Rename the count column
sample.columns = ["gene_id", os.path.basename(files[0]).replace(".txt", "")]

# 5. Add the other 29 samples
for file in files[1:]:

    df = pd.read_csv(file, sep="\t", comment="#")

    sample_name = os.path.basename(file).replace(".txt", "")

    current = df[["Geneid", df.columns[-1]]]

    current.columns = ["Geneid", sample_name]

    sample = sample.merge(current, on="gene_id")

# 6. Make Geneid the row name
sample = sample.set_index("gene_id")

# 7. Save the final matrix
sample.to_csv("config/count_matrix.csv")

print(sample)
print("Count matrix created!")
