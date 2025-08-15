############################################################
###### Script FST – Dados de SNPs ##########################
###### Aula prática – Ornitologia no Antropoceno ###########
###### Autora: F. Bocalini – Agosto 2025 ###################
############################################################

O **FST** mede a diferenciação genética entre populações.  
Seu valor varia de `0` (nenhuma diferenciação) a `1` (diferenciação completa).

Muitos estimadores diferentes de **FST** já foram propostos.  
Em dados reais, usamos com frequência estimadores mais complexos.  
Abaixo, vamos utilizar a função `stamppFst` do pacote **StAMPP** para aplicar  
o estimador de FST de **Weir & Cockerham (1984)** ao nosso conjunto real de dados de SNPs.  
Esse estimador leva em consideração o tamanho amostral e é um dos mais utilizados atualmente.

---
```r
# Carregar pacotes necessários
```

```r
library(vcfR)
```

```r
library(adegenet)
```

```r
#install.packages("StAMPP")
```

```r
library(StAMPP)
```

```r
# Definir o diretório de trabalho
```

```r
setwd("C:/Users/Fernanda Bocalini/Dropbox/pos-doc/aula_ornito_antropoceno/aula_pratica/FST")
```

```r
# Ler os dados de SNP no formato VCF
```

```r
my_vcf <- read.vcfR("xipho_95_total_final.vcf.gz")
```

```r
# Converter VCF para objeto genlight
```

```r
my_genlight <- vcfR2genlight(my_vcf)
```

```r
# Ler a tabela com as informações das amostras
```

```r
indmap <- read.table("popmap_xipho.csv", sep = ",", colClasses = "character", header = TRUE)
```

```r
# Caso a tabela não esteja na ordem dos indivíduos do VCF, reordenar
```

```r
indmap <- indmap[match(indNames(my_genlight), indmap['id']), ]
```

```r
# Definir a população de cada indivíduo
```

```r
pop(my_genlight) <- indmap[, "pop"]
```

```r
# Calcular FST usando a função 'stamppFst' do pacote 'StAMPP'
```

```r
fst_estimates <- stamppFst(my_genlight, nboots = 1000, percent = 95)
```

```r
# Ver resultados com quatro casas decimais
```

```r
round(fst_estimates$Fsts, 4)
```

```r
# Heatmap dos resultados - na pasta há uma função feita para plotar os resultados como heatmap
```

```r
source("plot_fst_heatmap.R")
```

```r
plot_fst_heatmap(fst_estimates$Fsts)
```

