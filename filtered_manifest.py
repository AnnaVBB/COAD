import pandas as pd

def filter_manifest(input_file, output_file):
    # Carregar os dados do manifesto
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Arquivo carregado: {input_file}")
        print(f"Número de linhas carregadas: {len(df)}")
    except Exception as e:
        print(f"Erro ao carregar o arquivo: {e}")
        return

    # Verificar se as colunas necessárias estão presentes
    if 'submitter_id' not in df.columns or 'tissue_type' not in df.columns:
        print("As colunas 'submitter_id' e 'tissue_type' devem estar presentes no arquivo.")
        return

    # Contar as ocorrências de cada submitter_id
    counts = df['submitter_id'].value_counts()
    print(f"Contagens de submitter_id:\n{counts}")

    # Filtrar para manter apenas os submitter_ids que aparecem mais de uma vez
    filtered_df = df[df['submitter_id'].isin(counts[counts > 1].index)]
    print(f"Número de linhas após a filtragem: {len(filtered_df)}")

    # Ordenar pelo submitter_id
    filtered_df = filtered_df.sort_values(by='submitter_id')

    # Criar um DataFrame para armazenar as amostras válidas
    valid_samples = []

    # Agrupar pelo submitter_id e selecionar amostras
    grouped = filtered_df.groupby('submitter_id')
    for submitter_id, group in grouped:
        # Verificar se temos uma amostra Normal e uma Tumoral
        normal_sample = group[group['tissue_type'] == 'Normal']
        tumoral_sample = group[group['tissue_type'] == 'Tumor']  # Ajustado para "Tumor"

        if not normal_sample.empty and not tumoral_sample.empty:
            # Adiciona as amostras Normais e Tumorais ao resultado
            valid_samples.append(normal_sample.iloc[0])  # Seleciona a primeira amostra Normal
            valid_samples.append(tumoral_sample.iloc[0])  # Seleciona a primeira amostra Tumoral

    # Verifique se temos amostras válidas
    if valid_samples:
        # Criar um DataFrame a partir das amostras válidas
        result_df = pd.DataFrame(valid_samples)
        # Salvar o DataFrame filtrado em um novo arquivo
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"Arquivo filtrado salvo como '{output_file}'.")
    else:
        print("Nenhuma amostra válida encontrada para salvar.")

def main():
    input_file = 'gdc_manifest_coad_one.txt.map2submitterID.tsv'  # Altere para o nome correto do arquivo
    output_file = 'manifest_coad_filtered.txt'
    
    filter_manifest(input_file, output_file)

if __name__ == "__main__":
    main()
