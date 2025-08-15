############################################################
###### Script DAPC – Dados de SNPs ##########################
###### Aula prática – Ornitologia no Antropoceno ###########
###### Autora: F. Bocalini – Agosto 2025 ###################
############################################################


# Qual é a diferença entre Análise de Componentes Principais (PCA) e Análise Discriminante (DA)? 
# PCA busca a direção de maior variância total, enquanto DA maximiza a separação entre grupos 
# minimizando a variação dentro de cada grupo.

# Carregar os pacotes
library(vcfR)
library(adegenet)
library(ggplot2)

# Definir o diretório de trabalho
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/DAPC")

# Ler o arquivo VCF no R
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")

# Transformar o VCF em objeto genind
my_genind <- vcfR2genind(my_vcf)
class(my_genind)  # Checar se a conversão foi feita corretamente

# Ler a tabela com a informação das amostras
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)

# Caso a tabela não esteja na ordem dos indivíduos do VCF, reordenar
indmap <- indmap[match(indNames(my_genind), indmap$id), ]

# Atribuir populações
pop(my_genind) <- indmap[, "pop"]  

############################################################
### Análise de DAPC
############################################################

# A DAPC testa grupos definidos a priori (e.g., subespécies, populações geográficas)
# Quando os agrupamentos são incertos, podemos usar o algoritmo K-means para identificar
# o número de agrupamentos genéticos diagnósticos na amostra

# Ao usar o K-means, diferentes números possíveis de clusters são comparados usando o BIC.
# O número ótimo de clusters corresponde ao menor valor de BIC.

# Rodar 'find.clusters' interativamente
my_grp <- find.clusters(my_genind, max.n.clust = 10)

# Rodar 'find.clusters' especificando parâmetros, após escolher número ótimo de clusters
my_grpk3 <- find.clusters(my_genind, n.pca = 40, n.clust = 3)

# Ver quais indivíduos pertencem a cada agrupamento
table(pop(my_genind), my_grpk3$grp)

# Criar tabela com população e grupo de cada indivíduo
my_grp_results <- data.frame(
  pop = pop(my_genind),
  groups_k3 = my_grpk3$grp
)

my_grp_results

# Salvar tabela
write.csv(my_grp_results, file = "find_cluster_results.csv")

##### Rodar o DAPC nesses agrupamentos #####

# Evitar usar muitas PCs para não overfittar (inflar separação entre grupos)
# Seguindo o "critério k-1", para k = 3, manteremos 2 PCs

dapc_grpk3 <- dapc(my_genind, pop = my_grpk3$grp, n.pca = 2, n.da = 2)

# Visualizar DAPC usando a função scatter do pacote 'adegenet'
scatter(dapc_grpk3)

# Customizar o gráfico
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

###### Gráfico usando 'ggplot' ######
col <- c(
  AF     = "#08bdbd", 
  AM     = "#29bf12", 
  Ceara  = "#ff9914", 
  PCE    = "#f21b3f"
)
# Extrair informações relevantes para um novo dataframe
plot_dapc_df <- data.frame(
  individual = row.names(dapc_grpk3$ind.coord),
  LD1 = dapc_grpk3$ind.coord[, 1],
  LD2 = dapc_grpk3$ind.coord[, 2],
  population = pop(my_genind),
  group = dapc_grpk3$grp
)

# Gráfico customizado com ggplot
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
