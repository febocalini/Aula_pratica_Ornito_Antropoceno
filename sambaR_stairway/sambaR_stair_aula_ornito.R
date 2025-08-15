############################################################
###### Script SambaR - ??ndices de Diversidade e Stairway Plots ######
###### Aula pr??tica ??? Ornitologia no Antropoceno ####################
###### Autora: F. Bocalini ??? Agosto 2025 ###########################
############################################################

# 1?????? Definir o diret??rio de trabalho
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/sambaR_stairway")

# 2?????? Carregar pacotes necess??rios
library(vcfR)
library(SNPfiltR)
library(adegenet)
library(StAMPP)
#install.packages("phangorn")
library(phangorn) # Para c??lculos de dist??ncia gen??tica

# 3?????? Carregar o script SambaR e instalar pacotes necess??rios
#source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.08.txt")
source("SAMBAR_v1.10.txt")
getpackages()

############################################################
# 4?????? Ler dados de SNP
############################################################

# Ler arquivo VCF
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")

# Converter para objeto genlight
my_genlight <- vcfR2genlight(my_vcf)

# Ler tabela com informa????es das amostras
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)

# Reordenar tabela para coincidir com a ordem do VCF
indmap <- indmap[match(indNames(my_genlight), indmap$id), ]

# Definir popula????o de cada indiv??duo
pop(my_genlight) <- indmap[, "pop"]

# Criar vetor das popula????es
pop_xipho <- as.vector(pop(my_genlight))

############################################################
# 5?????? Converter para formato SambaR
############################################################

genlight2sambar(
  genlight_object = "my_genlight",
  do_confirm = TRUE,
  popvector = pop_xipho,
  allow_edit = TRUE,
  colourvector= c("#08bdbd","#29bf12", "#ff9914", "#f21b3f")
)

############################################################
# 6?????? Filtragem do dataset
############################################################

# Permitir at?? 10% de dados ausentes, MAC m??nimo = 1, aplicar filtro de heterozigosidade
filterdata(
  indmiss = 0.9,
  snpmiss = 0.9,
  min_mac = 1,
  dohefilter = TRUE,
  maxprop_hefilter = 0.5,
  do_distplot = FALSE
)

############################################################
# 7?????? An??lises de Afinidade Gen??tica, Dist??ncia & Diversidade
############################################################

# Calcular dist??ncias gen??ticas (opcional: excluir genpop)
calcdistance(dodistgenpop = FALSE)

# Calcular parentesco (kinship) - importante para estimar inbreeding
calckinship()

# Calcular diversidade e gerar SFS dobrado para Stairway Plot
calcdiversity(
  do_sfs = TRUE,
  do_continue = TRUE,
  nrsites = 410190,     # Multiplicar o n??mero de loci por 110
  do_venn = FALSE
)

############################################################
# 8?????? Visualiza????o do Stairway Plot (ap??s an??lise no Linux)
############################################################

# Copiar o resultado final para a pasta demography antes de rodar no R
run_plotstairway(
  mu_rate = "2.5*10^-8",
  Gtime = "4.4",
  x_range = c(10, 500000), 
  pop_names = c("AM", "Ceara", "PCE"), 
  my_colours = c("#29bf12", "#ff9914", "#f21b3f")
)

############################################################
# 9?????? Exportar resultados para Excel
############################################################
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/sambaR_stairway")
write.csv(inds, "inds.csv")  # Exportar m??tricas individuais
write.csv(snps, "snps.csv")  # Exportar dados em n??vel de SNP
