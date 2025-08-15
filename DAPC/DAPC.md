# Script DAPC – Dados de SNPs
### Aula prática – Ornitologia no Antropoceno  
**Autora:** F. Bocalini – Agosto 2025

---

## Diferença entre PCA e DA
> **PCA** busca a direção de maior variância total, enquanto **DA** maximiza a separação entre grupos minimizando a variação dentro de cada grupo.

---

## Carregar pacotes
```r
library(vcfR)
library(adegenet)
library(ggplot2)
```

## Definir diretório de trabalho
```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/DAPC")
```

## Ler arquivo VCF e converter para `genind`
```r
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")
my_genind <- vcfR2genind(my_vcf)
class(my_genind)  # Checar conversão
```

## Ler tabela com informação das amostras
```r
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)
indmap <- indmap[match(indNames(my_genind), indmap$id), ]
pop(my_genind) <- indmap[, "pop"]
```

---

# Análise de DAPC

## Identificação de clusters
```r
my_grp <- find.clusters(my_genind, max.n.clust = 10)
my_grpk3 <- find.clusters(my_genind, n.pca = 40, n.clust = 3)
table(pop(my_genind), my_grpk3$grp)
```

## Criar e salvar tabela de resultados
```r
my_grp_results <- data.frame(
  pop = pop(my_genind),
  groups_k3 = my_grpk3$grp
)
write.csv(my_grp_results, file = "find_cluster_results.csv")
```

---

## Rodar DAPC
```r
dapc_grpk3 <- dapc(my_genind, pop = my_grpk3$grp, n.pca = 2, n.da = 2)
scatter(dapc_grpk3)
```

## Gráfico customizado com `scatter`
```r
col <- c("#08bdbd", "#29bf12", "#ff9914")
scatter(dapc_grpk3, 
        col=transp(col,0.9), 
        cellipse = TRUE,  
        pch = 19, 
        cex = 1.8, 
        legend = TRUE, 
        label = NULL, 
        posi.leg = "topright", 
        cleg = 0.7, 
        posi.da = "topleft")
```

---

## Gráfico usando `ggplot`
```r
col <- c(
  AF     = "#08bdbd", 
  AM     = "#29bf12", 
  Ceara  = "#ff9914", 
  PCE    = "#f21b3f"
)

plot_dapc_df <- data.frame(
  individual = row.names(dapc_grpk3$ind.coord),
  LD1 = dapc_grpk3$ind.coord[, 1],
  LD2 = dapc_grpk3$ind.coord[, 2],
  population = pop(my_genind),
  group = dapc_grpk3$grp
)

ggplot(data = plot_dapc_df, aes(x = LD1, y = LD2, fill = population)) +
  geom_point(shape = 21, color = "black", size = 6, alpha = 0.9) +
  scale_fill_manual(values = col) +
  theme_classic() +
  labs(
    x = "LD1",
    y = "LD2",
    fill = NULL
  ) +
  theme(
    legend.text = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 18)
  )
```
