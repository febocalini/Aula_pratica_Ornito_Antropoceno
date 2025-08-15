############################################################
###### Script SambaR - Índices de Diversidade e Stairway Plots ######
###### Aula prática – Ornitologia no Antropoceno ####################
###### Autora: F. Bocalini – Agosto 2025 ###########################
############################################################

---
```r
# 1️⃣ Definir o diretório de trabalho
```

```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/sambaR_stairway")
```

```r
# 2️⃣ Carregar pacotes necessários
```

```r
library(vcfR)
```

```r
library(SNPfiltR)
```

```r
library(adegenet)
```

```r
library(StAMPP)
```

```r
#install.packages("phangorn")
```

```r
library(phangorn) # Para cálculos de distância genética
```

```r
# 3️⃣ Carregar o script SambaR e instalar pacotes necessários
```

```r
#source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.08.txt")
```

```r
source("SAMBAR_v1.10.txt")
```

```r
getpackages()
```

```r
############################################################
```

```r
# 4️⃣ Ler dados de SNP
```

```r
############################################################
```

```r
# Ler arquivo VCF
```

```r
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")
```

```r
# Converter para objeto genlight
```

```r
my_genlight <- vcfR2genlight(my_vcf)
```

```r
# Ler tabela com informações das amostras
```

```r
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)
```

```r
# Reordenar tabela para coincidir com a ordem do VCF
```

```r
indmap <- indmap[match(indNames(my_genlight), indmap$id), ]
```

```r
# Definir população de cada indivíduo
```

```r
pop(my_genlight) <- indmap[, "pop"]
```

```r
# Criar vetor das populações
```

```r
pop_xipho <- as.vector(pop(my_genlight))
```

```r
############################################################
```

```r
# 5️⃣ Converter para formato SambaR
```

```r
############################################################
```

```r
genlight2sambar(
```

```r
  genlight_object = "my_genlight",
```

```r
  do_confirm = TRUE,
```

```r
  popvector = pop_xipho,
```

```r
  allow_edit = TRUE,
```

```r
  colourvector= c("#08bdbd","#29bf12", "#ff9914", "#f21b3f")
```

```r
)
```

```r
############################################################
```

```r
# 6️⃣ Filtragem do dataset
```

```r
############################################################
```

```r
# Permitir até 10% de dados ausentes, MAC mínimo = 1, aplicar filtro de heterozigosidade
```

```r
filterdata(
```

```r
  indmiss = 0.9,
```

```r
  snpmiss = 0.9,
```

```r
  min_mac = 1,
```

```r
  dohefilter = TRUE,
```

```r
  maxprop_hefilter = 0.5,
```

```r
  do_distplot = FALSE
```

```r
)
```

```r
############################################################
```

```r
# 7️⃣ Análises de Afinidade Genética, Distância & Diversidade
```

```r
############################################################
```

```r
# Calcular distâncias genéticas (opcional: excluir genpop)
```

```r
calcdistance(dodistgenpop = FALSE)
```

```r
# Calcular parentesco (kinship) - importante para estimar inbreeding
```

```r
calckinship()
```

```r
# Calcular diversidade e gerar SFS dobrado para Stairway Plot
```

```r
calcdiversity(
```

```r
  do_sfs = TRUE,
```

```r
  do_continue = TRUE,
```

```r
  nrsites = 410190,     # Multiplicar o número de loci por 110
```

```r
  do_venn = FALSE
```

```r
)
```

```r
############################################################
```

```r
# 8️⃣ Visualização do Stairway Plot (após análise no Linux)
```

```r
############################################################
```

```r
# Copiar o resultado final para a pasta demography antes de rodar no R
```

```r
run_plotstairway(
```

```r
  mu_rate = "2.5*10^-8",
```

```r
  Gtime = "4.4",
```

```r
  x_range = c(10, 500000), 
```

```r
  pop_names = c("AM", "Ceara", "PCE"), 
```

```r
  my_colours = c("#29bf12", "#ff9914", "#f21b3f")
```

```r
)
```

```r
############################################################
```

```r
# 9️⃣ Exportar resultados para Excel
```

```r
############################################################
```

```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/sambaR_stairway")
```

```r
write.csv(inds, "inds.csv")  # Exportar métricas individuais
```

```r
write.csv(snps, "snps.csv")  # Exportar dados em nível de SNP
```

