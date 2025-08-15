############################################################
###### Script Kinship no SambaR e Plot dos Índices ######
###### Aula prática – Ornitologia no Antropoceno ####################
###### Autora: F. Bocalini – Agosto 2025 ###########################
############################################################

# Utilizando dados de T. melanonotus do artigo:
# Carcassola, M. V, F. Bocalini, L. F. Silveira, M. R. Francisco, e A. G. A. Migotto (2025).
# Unraveling the genetic parameters, social relationships, and conservation aspects of a group of Touit melanonotus (Brown-backed Parrotlet).
# Ornithological Applications:1–32.

############################################################
# 1️⃣ Definir o diretório de trabalho
############################################################
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/kinship")

############################################################
# 2️⃣ Carregar pacotes necessários
############################################################
library(vcfR)
library(SNPfiltR)
library(adegenet)
library(StAMPP)

############################################################
# 3️⃣ Carregar o script SambaR e instalar pacotes necessários
############################################################
#source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.08.txt")
source("SAMBAR_v1.10.txt")
getpackages()

############################################################
# 4️⃣ Ler dados de SNP
############################################################

# Ler arquivo VCF
my_vcf_touit <- read.vcfR("touit_snps75.vcf")

# Converter para objeto genlight
my_genlight_touit <- vcfR2genlight(my_vcf_touit)

# Ler tabela com informações das amostras
indmap <- read.table("pop_info_touit.csv", sep = ",", colClasses = "character", header = TRUE)

# Definir população de cada indivíduo
pop(my_genlight_touit) <- indmap[, "Pop"]

# Criar vetor das populações
pop_touit <- as.vector(pop(my_genlight_touit))

############################################################
# 5️⃣ Converter para formato SambaR
############################################################

genlight2sambar(
  genlight_object = "my_genlight_touit",
  do_confirm = TRUE,
  popvector = pop_touit,
  allow_edit = TRUE,
  colourvector= c("#29bf12")
)

############################################################
# 6️⃣ Filtragem do dataset
############################################################

# Permitir até 10% de dados ausentes, MAC mínimo = 2, aplicar filtro de heterozigosidade
filterdata(
  indmiss = 0.9,
  snpmiss = 0.9,
  min_mac = 2,
  dohefilter = TRUE,
  maxprop_hefilter = 0.5,
  do_distplot = FALSE
)

############################################################
# 7️⃣ Calcular índices de Kinship
############################################################

calckinship()

# Os resultados estarão salvos na pasta 'kinship' no arquivo 'pairwise_relatedness.txt'
# Abra a planilha no Excel e salve apenas as colunas de interesse

############################################################
# 8️⃣ Exemplo de plot das relações
############################################################

# Pairs with KING-robust scores around 0.25 or higher are considered full siblings (FS).
# Scores between 0.125 and 0.25 indicate second-degree relatives (SD).
# Scores between 0.0625 and 0.125 indicate third-degree relatives (TD).
# Scores below 0.0625 indicate unrelated pairs (U).
# To differentiate FS from parent-offspring (PO), the k0 statistic is used,
# representing the proportion of loci where the pair shares zero alleles identical by state.
# If k0 is below 0.005, the pair is considered parent-offspring (PO).

# Carregar pacotes necessarios
library(ggplot2)
library(dplyr)

# Ler o arquivo CSV (substitua pelo caminho correto do arquivo)
data <- read.csv("relatedness_touit.csv")

# Classifique o relacionamento dos pares de acordo com o kigrobust e k0
data <- data %>%
  mutate(relationship = case_when(
    kingrobust >= 0.2 & kingrobust <= 0.5 & k0 < 0.005 ~ "Parent Offspring",
    kingrobust >= 0.2 & kingrobust <= 0.5  ~ "Full-Siblings",
    kingrobust >= 0.1 & kingrobust <= 0.199  ~ "Second-Degree Relatives",
    kingrobust >= 0.05 & kingrobust <= 0.099 ~ "Third-Degree Relatives",
    kingrobust < 0.05 ~ "Unrelated"
  ))

# Criar um ggplot das relações dos pares:
ggplot(data, aes(x = k0, y = kingrobust, color = relationship)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "green", size = 0.6) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue", size = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "purple", size = 0.6) +
  labs(title = "", x = "k0", y = "KING-Robust", color = "") +
  theme_minimal() +
  theme(text = element_text(size = 10)) +
  scale_color_manual(values = c(
    "Parent Offspring" = "darkred",
    "Full-Siblings" = "darkgreen",
    "Second-Degree Relatives" = "blue",
    "Third-Degree Relatives" = "purple",
    "Unrelated" = "gray"
  ))

# Savar a planilha com a descrição dos relacionamento de cada par:
write.csv(data, "relatedness_with_relationships_touit.csv", row.names = FALSE)
