############################################################
###### Script SambaR - Índices de Diversidade e Stairway Plots ######
###### Aula prática – Ornitologia no Antropoceno ####################
###### Autora: F. Bocalini – Agosto 2025 ###########################
############################################################

# 1️⃣ Definir o diretório de trabalho
```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/sambaR_stairway")

```
# 2️⃣ Carregar pacotes necessários
```r
library(vcfR)
library(SNPfiltR)
library(adegenet)
library(StAMPP)
#install.packages("phangorn")
library(phangorn) # Para cálculos de distância genética
```

# 3️⃣ Carregar o script SambaR e instalar pacotes necessários

```r
#source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.08.txt")
source("SAMBAR_v1.10.txt")
getpackages()

```

# 4️⃣ Ler dados de SNP


# Ler arquivo VCF
```r
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")

```
# Converter para objeto genlight
```r
my_genlight <- vcfR2genlight(my_vcf)

```
# Ler tabela com informações das amostras
```r
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)

```
# Reordenar tabela para coincidir com a ordem do VCF
```r
indmap <- indmap[match(indNames(my_genlight), indmap$id), ]

```
# Definir população de cada indivíduo
```r
pop(my_genlight) <- indmap[, "pop"]

```
# Criar vetor das populações
```r
pop_xipho <- as.vector(pop(my_genlight))

```

# 5️⃣ Converter para formato SambaR


```r
genlight2sambar(
  genlight_object = "my_genlight",
  do_confirm = TRUE,
  popvector = pop_xipho,
  allow_edit = TRUE,
  colourvector= c("#08bdbd","#29bf12", "#ff9914", "#f21b3f")
)

```
############################################################
# 6️⃣ Filtragem do dataset
############################################################

# Permitir até 10% de dados ausentes, MAC mínimo = 1, aplicar filtro de heterozigosidade
```r
filterdata(
  indmiss = 0.9,
  snpmiss = 0.9,
  min_mac = 1,
  dohefilter = TRUE,
  maxprop_hefilter = 0.5,
  do_distplot = FALSE
)

```
############################################################
# 7️⃣ Análises de Afinidade Genética, Distância & Diversidade
############################################################

# Calcular distâncias genéticas (opcional: excluir genpop)
```r
calcdistance(dodistgenpop = FALSE)

```
# Calcular parentesco (kinship) - importante para estimar inbreeding
```r
calckinship()

```
# Calcular diversidade e gerar SFS dobrado para Stairway Plot
```r
calcdiversity(
  do_sfs = TRUE,
  do_continue = TRUE,
  nrsites = 410190,     # Multiplicar o número de loci por 110
  do_venn = FALSE
)

```
############################################################
# 8️⃣ Visualização do Stairway Plot (após análise no Linux)
############################################################

# Copiar o resultado final para a pasta demography antes de rodar no R
```r
run_plotstairway(
  mu_rate = "2.5*10^-8",
  Gtime = "4.4",
  x_range = c(10, 500000), 
  pop_names = c("AM", "Ceara", "PCE"), 
  my_colours = c("#29bf12", "#ff9914", "#f21b3f")
)

```
############################################################
# 9️⃣ Exportar resultados para Excel
############################################################
```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/sambaR_stairway")
write.csv(inds, "inds.csv")  # Exportar métricas individuais
write.csv(snps, "snps.csv")  # Exportar dados em nível de SNP
```
