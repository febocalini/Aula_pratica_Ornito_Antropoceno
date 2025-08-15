############################################################
###### Script FST – Dados de SNPs ##########################
###### Aula prática – Ornitologia no Antropoceno ###########
###### Autora: F. Bocalini – Agosto 2025 ###################
############################################################

# O FST mede a diferenciação genética entre populações.
# Seu valor varia de 0 (nenhuma diferenciação) a 1 (diferenciação completa).

# Muitos estimadores diferentes de FST já foram propostos.
# Em dados reais, usamos com frequência estimadores mais complexos.
# Abaixo, vamos utilizar a função ‘stamppFst’ do pacote ‘StAMPP’
# para aplicar o estimador de FST de Weir & Cockerham (1984) 
# ao nosso conjunto real de dados de SNPs.
# Esse estimador leva em consideração o tamanho amostral
# e é um dos mais utilizados atualmente.

# Carregar pacotes necessários
library(vcfR)
library(adegenet)
#install.packages("StAMPP")
library(StAMPP)

# Definir o diretório de trabalho
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/FST")

# Ler os dados de SNP no formato VCF
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")

# Converter VCF para objeto genlight
my_genlight <- vcfR2genlight(my_vcf)

# Ler a tabela com as informações das amostras
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)

# Caso a tabela não esteja na ordem dos indivíduos do VCF, reordenar
indmap <- indmap[match(indNames(my_genlight), indmap$id), ]

# Definir a população de cada indivíduo
pop(my_genlight) <- indmap[, "pop"]

# Calcular FST usando a função 'stamppFst' do pacote 'StAMPP'
fst_estimates <- stamppFst(my_genlight, nboots = 1000, percent = 95)

# Ver resultados com quatro casas decimais
round(fst_estimates$Fsts, 4)

# Heatmap dos resultados - na pasta há uma função feita para plotar os resultados como heatmap
source("plot_fst_heatmap.R")

plot_fst_heatmap(fst_estimates$Fsts)
