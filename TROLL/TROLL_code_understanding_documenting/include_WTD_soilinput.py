import pandas as pd

original_file = "Paracou_input_pedology.txt"
modified_file = "Paracou_input_pedology_WTD.txt"

# read TXT as a dataframe 
df = pd.read_csv(original_file, sep="\t", dtype=str, keep_default_na=False)

#creates a new column
df["WTD"] = [1.0, 1.0, 1.0, 1.0, 1.0] 

# Inserindo a coluna "WTD" na posição desejada (por exemplo, posição 8)
# A posição deve ser um índice baseado em zero, então 8 significa que a coluna será inserida entre as colunas 8 e 9.
df.insert(8, "WTD", df.pop("WTD"))

# Salvando o novo arquivo
df.to_csv(modified_file, sep="\t", index=False)

print("Succed")