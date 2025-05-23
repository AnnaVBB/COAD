import pandas as pd

# Caminhos dos arquivos
manifesto_filtrado_path = "C:/Users/Usuário/Downloads/TCGA/manifest_coad_filtered.txt"  # Manifesto filtrado
manifesto_original_path = "C:/Users/Usuário/Downloads/TCGA/gdc_manifest_coad_one.txt"  # Manifesto original
output_manifesto_path = "C:/Users/Usuário/Downloads/TCGA/manifesto_coad_reb.txt"  # Caminho do arquivo de saída

# Carregar o manifesto filtrado e o manifesto original
manifesto_filtrado = pd.read_csv(manifesto_filtrado_path, sep='\t')
manifesto_original = pd.read_csv(manifesto_original_path, sep='\t')

# Realizar a união (merge) dos dados usando 'file_id' do manifesto filtrado e 'id' do manifesto original
merged_manifesto = pd.merge(manifesto_original, manifesto_filtrado, left_on='id', right_on='file_id', how='inner')

# Selecionar apenas as colunas relevantes do manifesto original
colunas_obrigatorias = ['id', 'filename', 'md5', 'size', 'state']
merged_manifesto = merged_manifesto[colunas_obrigatorias]

# Salvar o novo manifesto
if not merged_manifesto.empty:
    merged_manifesto.to_csv(output_manifesto_path, sep='\t', index=False)
    print(f"Manifesto unido salvo em: {output_manifesto_path}")
else:
    print("A união resultou em um DataFrame vazio. Verifique os valores de 'file_id' e 'id'.")