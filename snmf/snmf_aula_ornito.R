############################################################
###### Script sNMF – Dados de SNPs #########################
###### Aula prática – Ornitologia no Antropoceno ###########
###### Autora: F. Bocalini – Agosto 2025 ###################
############################################################

# Carregar os pacotes
#install.packages("BiocManager")
# BiocManager::install("LEA")
library(LEA)
library(vcfR)

# Definir o diretório de trabalho
setwd("C:/brejos_analises_novas/snmf")

# Preparar dados

# Ler o arquivo VCF no R
my_vcf_un <- read.vcfR("xipho_95_unlinked_final.vcf.gz")

# Transformar VCF em genind
genind_un <- vcfR2genind(my_vcf_un)
class(genind_un)  # Checar se a transformação foi feita corretamente

# Estatísticas básicas
nLoc(genind_un)   # Número de locos
nInd(genind_un)   # Número de indivíduos

##### Transformar o genind em geno (formato de entrada no LEA) #####

write.table(genind_un, file = "xipho_un.nex")  # Exportar o objeto genind como tabela

# Transformar a tabela em matriz
table <- read.table("xipho_un.nex")
matrix <- as.matrix(table)

# Transformar a matriz em formato .geno
write.geno(matrix, output.file = "xipho_un.geno")
read.geno("xipho_un.geno")

##### Outra forma (prefiro a de cima) #####
vcf2geno(input.file = "xipho_95_unlinked_final.vcf",
         output.file = "xipho_un2.geno")

###### Rodar sNMF para cada valor de alpha ######

# O que é alpha?
# Alpha é um parâmetro de regularização não negativo que penaliza coeficientes de ancestralidade intermediários.
# Valores de alpha > 0 tendem a reduzir a variância das estimativas de coeficiente de ancestralidade
# e forçar estimativas irrelevantes a zero.
# Frichot et al. (2014) mostraram que esse parâmetro influencia mais em conjuntos de dados genômicos pequenos (< 10.000 SNPs).
# Em geral, valores grandes de alpha melhoram as estimativas quando o conjunto de dados é pequeno.
# Recomenda-se testar múltiplos valores de alpha e escolher aquele com menor valor de entropia cruzada (cross-entropy).
# Menores valores de entropia cruzada indicam melhores estimativas.

# OBS: normalmente rodamos o sNMF com 100 repetições ou mais,
# mas isso pode levar 20–40 min dependendo do computador.

# alpha = 10
project_xipho_total_a10 <- snmf("xipho_un.geno", K = 1:10, entropy = TRUE,
                                repetitions = 20, alpha = 10, CPU = 2, seed = 1,
                                project = "new")

# alpha = 50
project_xipho_total_a50 <- snmf("xipho_un.geno", K = 1:10, entropy = TRUE,
                                repetitions = 20, alpha = 50, CPU = 2, seed = 2,
                                project = "new")

# alpha = 100
project_xipho_total_a100 <- snmf("xipho_un.geno", K = 1:10, entropy = TRUE,
                                 repetitions = 20, alpha = 100, CPU = 2, seed = 3,
                                 project = "new")

# alpha = 500
project_xipho_total_a500 <- snmf("xipho_un.geno", K = 1:10, entropy = TRUE,
                                 repetitions = 20, alpha = 500, CPU = 2, seed = 4,
                                 project = "new")

##### Carregar os resultados #####

# Extrair valores de entropia cruzada
sum_total_a10   <- summary(project_xipho_total_a10)$crossEntropy
sum_total_a50   <- summary(project_xipho_total_a50)$crossEntropy
sum_total_a100  <- summary(project_xipho_total_a100)$crossEntropy
sum_total_a500  <- summary(project_xipho_total_a500)$crossEntropy

# Plot básico
plot(project_xipho_total_a10,   lwd = 5, col = "blue4",     pch = 19)
plot(project_xipho_total_a100,  lwd = 5, col = "blue",      pch = 19)
plot(project_xipho_total_a50,   lwd = 5, col = "cyan",      pch = 19)
plot(project_xipho_total_a500,  lwd = 5, col = "lightblue", pch = 19)

# Criar dataframe com entropias
k_snmf_total <- data.frame(
  alpha        = rep(c(10, 50, 100, 500), each = 10),
  k            = rep(1:10, 4),
  min.entropy  = c(sum_total_a10[1, ], sum_total_a50[1, ],
                   sum_total_a100[1, ], sum_total_a500[1, ]),
  mean.entropy = c(sum_total_a10[2, ], sum_total_a50[2, ],
                   sum_total_a100[2, ], sum_total_a500[2, ]),
  max.entropy  = c(sum_total_a10[3, ], sum_total_a50[3, ],
                   sum_total_a100[3, ], sum_total_a500[3, ])
)

write.csv(k_snmf_total, file = "kvalue_xipho_aula.csv")

alpha.labs <- c("\u03b1 = 10", "\u03b1 = 50", "\u03b1 = 100", "\u03b1 = 500")
names(alpha.labs) <- c(10, 50, 100, 500)

# Plot no ggplot2
ggplot(k_snmf_total, aes(x = k, y = min.entropy)) +
  geom_line(size = 0.9) +
  geom_point(size = 4) +
  facet_wrap(~alpha, scales = "free_y", nrow = 3,
             labeller = labeller(alpha = alpha.labs)) +
  xlab("Número de clusters (K)") +
  ylab("Entropia cruzada") +
  scale_x_continuous(breaks = 1:15) +
  theme_bw() +
  theme(
    strip.text.x  = element_text(size = 15, face = "bold"),
    axis.title    = element_text(size = 15, color = "black"),
    axis.text     = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Melhor K = 3, alpha = 10

#-----------------------------------#
### Resultados para alpha = 10, K = 3
#-----------------------------------#

# Identificar a melhor execução para o K escolhido
ce_a10_k3 <- cross.entropy(project_xipho_total_a10, K = 3)
bestrun_a10_k3 <- which.min(ce_a10_k3)

# Ver o valor mínimo de entropia cruzada para o K escolhido
min(bestrun_a10_k3)

# Obter a Q-matrix com proporções de ancestralidade da melhor execução
Qmatrix_a10_k3 <- Q(project_xipho_total_a50, K = 3, run = bestrun_a10_k3)

write.csv(Qmatrix_a10_k3, file = "Qmatrix_a10_k3_xipho.csv", row.names = TRUE)

# Plot da estrutura genética
col <- c("#29bf12", "#08bdbd", "#f21b3f", "#ff9914")

barchart(project_xipho_total_a10, K = 3, run = bestrun_a10_k3, col = col,
         xlab = "Indivíduos",
         ylab = "Coeficientes de mistura (admixture)",
         las = 2, border = TRUE, space = 0) -> bp

# Criar legenda dos indivíduos
indmap <- read.table("popmap_xipho.csv", sep = ",", header = TRUE)
head(indmap)

pop_xipho <- indmap[, "pop"]
ind <- seq(1:length(pop_xipho))
pops <- data.frame(ind, pop_xipho)

ordem <- c(bp$order)
ordem.nova <- factor(pops$ind, levels = ordem)
xlab <- pops[order(ordem.nova), ]
xlab <- xlab$pop_xipho

col <- c("#f21b3f", "#29bf12", "#08bdbd")

# Plot com nomes corretos dos indivíduos no eixo x
barchart(project_xipho_total_a10, K = 3, run = bestrun_a10_k3, col = col,
         ylab = "Coeficientes de mistura (admixture)",
         border = TRUE, las = 2, space = 0.1, names.arg = xlab,
         cex.names = 1, cex.axis = 1.2, cex.lab = 1.5, font = 2)

#### Como colocar o nome das amostras no eixo x?
#### Faça o gráfico para o segundo melhor valor de K que você obteve.
