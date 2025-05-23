# Análise estatística descritiva e inferencial do câncer de cólon (COAD)
Análise estatística de tecidos normais vs tumorais do câncer de cólon (COAD)

  Através da base de dados do portal NIH-GDC-TCGA foram coletados dados públicos de expressão RNA-seq de pacientes com câncer de colorretal (COAD) com o objetivo de se realizar uma análise estatística descritiva e inferencial comparando amostras teciduais normais versus tumorais. Queremos analisar quais são os genes mais expressos e como a idade, o gênero e a expressão gênica influenciam o estágio do câncer e o risco de mortalidade deste câncer. 
  Foram selecionadas 82 amostras que possuíam tanto um tecido normal, quanto um tumoral. As variáveis clínicas de interesse foram: project.project_id, case.case_id, cases.submitter_id, demographic.age_at_index, demographic.gender, demographic.vital_status, diagnoses.ajcc_pathologic_stage 
  Para filtragem e remodelação do manifesto foi utilizado o pacote Pandas do Python e o GDC-Client-Tools. Para as análises estatísticas foram utilizados pacotes do R como o DESeq2 (para realizar análises de expressão gênica diferencial com dados de RNA-Seq), Bioconductor, dplyr, tidyverse, pheatmap, RColorBrewer, ggplot2 e readr

Steps
1) Baixar o manifesto (Filtrar por RNA-seq e STAR-Counts) {GDC portal}
2) Baixar os dados clínicos e criar uma tabela com as colunas de interesse (clinical_data.R)
3) Baixar o manifesto com o tipo de amostra e barcode (generate_manifest.py)
  Executar no terminal (exemplo):
  python C:\\Users\\Usuário\\Downloads\\TCGA\\generate_manifest.py -i C:\\Users\\Usuário\\Downloads\\TCGA\\gdc_manifest_coad_one.txt
4) Filtrar o manifesto, pareando uma amostra normal com uma tumoral e excluindo as demais (filtered_manifest.py)
5) Excluir as colunas sample_type, id e barcode (rebuild_manifest.py)
6) Baixar os dados de expressão gênica 
   Executar no terminal:
  .\gdc-client download -m manifesto_coad_reb.txt
7) Matriz de Contagem RNA-Seq, de modo que os id's estejam alinhados ao manifesto (build_matrix.R)
  Deve gerar o arquivo TCGA_COAD_RNASEQ_ORDERED 
8) mapear os IDs do manifesto para os IDs da matriz de contagem, criar a matriz sample_conditions para análise diferencial (build_matrix.R)
9) Análises estatísticas inferenciais com DESeq (build_matrix.R)
10) Análises estatísticas descritivas e inferenciais com os dados clínicos (clinical_data.R)
11) Correlação entre a expressão gênica e os dados clínicos (clinical_data.R)
