# Script PCA – Dados de SNPs
**Aula prática – Ornitologia no Antropoceno**  
**Autora:** F. Bocalini – Agosto 2025

---

## 1) Definir diretório de trabalho
```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/PCA")
```

## 2) Instalar e carregar pacotes necessários
```r
# install.packages("adegenet")
library(adegenet)

# install.packages("vcfR")
library(vcfR)

# install.packages("tidyr")
library(tidyr)
library(tibble)

# install.packages("tidyverse")
library(tidyverse)

# (opcional) ggrepel para rótulos nos gráficos
# install.packages("ggrepel")
# library(ggrepel)
```
---

## 3) Leitura e preparação dos dados
```r
# Ler o arquivo VCF no R
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")

# Transformar o VCF em objeto genind
my_genind <- vcfR2genind(my_vcf)
class(my_genind)  # Checar se a conversão foi feita corretamente

# Matriz de genótipos (homozigotos=2, heterozigotos=1)
tab(my_genind)[1:5, 1:10]

# Nomes de alelos por locus
my_genind@all.names %>% head()

# Fator com nomes dos loci
locFac(my_genind) %>% head()

# Nomes das amostras
ind_names <- indNames(my_genind)
ind_names
write.csv(ind_names, file = "ind_names.csv")  # Salvar nomes (ordem importa)

# Ploidia
ploidy(my_genind)

# Estatísticas básicas: número de loci e indivíduos
nLoc(my_genind)  # ex.: 3141
nInd(my_genind)  # ex.: 37

# Populações (carregar tabela e atribuir)
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)
pop(my_genind) <- indmap[, "pop"]
pop(my_genind)            # Checar
seppop(my_genind)         # Estatísticas por população
```
---

## 4) Análise de PCA
```r
# Frequências alélicas (substitui NA pela média)
allele_freq <- tab(my_genind, freq = TRUE, NA.method = "mean")

# PCA (retendo 10 eixos)
my_pca <- dudi.pca(allele_freq, scale = FALSE, scannf = FALSE, nf = 10)

# Eigenvalues
barplot(my_pca$eig, main = "PCA Eigenvalues")
barplot(my_pca$eig[1:5], main = "PCA Eigenvalues")

# Proporção de variância explicada
eig_perc <- my_pca$eig / sum(my_pca$eig)

# Data frame para gráficos
my_pca_df <- data.frame(
  individuo  = row.names(my_pca$li),
  PC1        = my_pca$li[, 1],
  PC2        = my_pca$li[, 2],
  PC3        = my_pca$li[, 3],
  populacao  = pop(my_genind)
)

# Paleta de cores por população
col <- c(
  AF     = "#08bdbd", 
  AM     = "#29bf12", 
  Ceara  = "#ff9914", 
  PCE    = "#f21b3f"
)
```
---

## 5) PCA no Base R
```r
pops  <- pop(my_genind)
l_pop <- levels(pops)

par(xpd = TRUE)
s.class(
  my_pca$li,
  factor(pops, levels = l_pop),
  xax = 1, yax = 2,
  col = transp(col, 0.7),
  axesell = FALSE, cstar = 0, cpoint = 2, grid = FALSE,
  label = NULL, cellipse = TRUE, pch = 19
)

legend("bottomleft", yjust = 1, bty = "n", cex = 0.6,
       title = "PCA", l_pop, fill = col)
```
---

## 6) PCA PC1 vs PC2 com ggplot2
```r
ggplot(data = my_pca_df, aes(x = PC1, y = PC2, color = indmap$pop)) +
  geom_point(size = 6, alpha = 0.9, stroke = 1.2, color = "black",
             aes(fill = indmap$pop), shape = 21) +
  scale_fill_manual(values = col) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(eig_perc[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_perc[2] * 100, 1), "%)"),
    fill = NULL
  ) +
  theme(
    legend.text = element_text(face = "bold", size = 14),
    axis.text   = element_text(face = "bold", size = 16),
    axis.title  = element_text(face = "bold", size = 18)
  )
```
---

## 7) PCA PC2 vs PC3 com ggplot2
```r
ggplot(data = my_pca_df, aes(x = PC2, y = PC3, color = indmap$pop)) +
  geom_point(size = 6, alpha = 0.9, stroke = 1.2, color = "black",
             aes(fill = indmap$pop), shape = 21) +
  scale_fill_manual(values = col) +
  theme_classic() +
  labs(
    x = paste0("PC2 (", round(eig_perc[2] * 100, 1), "%)"),
    y = paste0("PC3 (", round(eig_perc[3] * 100, 1), "%)"),
    fill = NULL
  ) +
  theme(
    legend.text = element_text(face = "bold.italic", size = 14),
    axis.text   = element_text(face = "bold", size = 16),
    axis.title  = element_text(face = "bold", size = 18)
  )
```
---

## 8) PCA com nomes das amostras (rótulos)
```r
library(ggrepel)

# Espaço extra nos eixos para rótulos
x_range   <- range(my_pca_df$PC1)
y_range   <- range(my_pca_df$PC2)
x_padding <- 0.1 * diff(x_range)
y_padding <- 0.1 * diff(y_range)

ggplot(my_pca_df, aes(x = PC1, y = PC2, fill = indmap$pop)) +
  geom_point(size = 6, alpha = 0.9, stroke = 1.2, color = "black", shape = 21) +
  geom_text_repel(
    aes(label = indmap$id),
    size = 1.5, color = "black",
    box.padding = 0.5, max.overlaps = Inf
  ) +
  scale_fill_manual(values = col) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(eig_perc[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_perc[2] * 100, 1), "%)"),
    fill = NULL
  ) +
  coord_cartesian(
    xlim = range(my_pca_df$PC1) * 1.1,
    ylim = range(my_pca_df$PC2) * 1.1
  )
```

---
**FIM**
