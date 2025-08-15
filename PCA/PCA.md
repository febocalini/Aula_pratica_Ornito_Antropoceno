# Script PCA – Dados de SNPs  
**Aula prática – Ornitologia no Antropoceno**  
Autora: F. Bocalini – Agosto 2025  

---

## Definir diretório de trabalho
```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/PCA")
Instalar e carregar pacotes necessários
r
Copiar
Editar
# install.packages("adegenet")
library(adegenet)

# install.packages("vcfR")
library(vcfR)

# install.packages("tidyr")
library(tidyr)
library(tibble)

# install.packages("tidyverse")
library(tidyverse)
Leitura e preparação dos dados
r
Copiar
Editar
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")
my_genind <- vcfR2genind(my_vcf)
class(my_genind)  # Checar conversão

tab(my_genind)[1:5, 1:10]
my_genind@all.names %>% head()
locFac(my_genind) %>% head()

ind_names <- indNames(my_genind)
ind_names
write.csv(ind_names, file = "ind_names.csv")

ploidy(my_genind)
nLoc(my_genind)
nInd(my_genind)

indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)
pop(my_genind) <- indmap[, "pop"]
pop(my_genind)
seppop(my_genind)
Análise de PCA
r
Copiar
Editar
allele_freq <- tab(my_genind, freq = TRUE, NA.method = "mean")
my_pca <- dudi.pca(allele_freq, scale = FALSE, scannf = FALSE, nf = 10)

barplot(my_pca$eig, main = "PCA Eigenvalues")
barplot(my_pca$eig[1:5], main = "PCA Eigenvalues")

eig_perc <- my_pca$eig / sum(my_pca$eig)

my_pca_df <- data.frame(
  individuo  = row.names(my_pca$li),
  PC1        = my_pca$li[, 1],
  PC2        = my_pca$li[, 2],
  PC3        = my_pca$li[, 3],
  populacao  = pop(my_genind)
)

col <- c(
  AF     = "#08bdbd", 
  AM     = "#29bf12", 
  Ceara  = "#ff9914", 
  PCE    = "#f21b3f"
)
PCA no Base R
r
Copiar
Editar
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
PCA PC1 vs PC2 com ggplot2
r
Copiar
Editar
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
PCA PC2 vs PC3 com ggplot2
r
Copiar
Editar
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
PCA com nomes das amostras
r
Copiar
Editar
library(ggrepel)

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
