#Script para obtener Diversidad alfa con Hill numbers

#Este parámetro se mide a nivel de especies (Leve7)
#Selecciono la carpeta de trabajo
setwd("~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity")

#Cargo librerías necesarias
library(readr)
library(hillR)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(devtools)

#Importo mi dataframe
#Los datos a importar deben verse así:

Df_batan_comm <- read.csv("batan.csv", header = TRUE, row.names = 1, fileEncoding = 'UTF-8-BOM')
View(Df_batan_comm)

q0 <- hill_taxa(Df_batan_comm, q=0, MARGIN = 1)
q1 <- hill_taxa(Df_batan_comm, q=1, MARGIN = 1)
q2 <- hill_taxa(Df_batan_comm, q=2, MARGIN = 1)

write.table(x=d0, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_q0.tsv")
write.table(x=d1, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_q1.tsv")
write.table(x=d2, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_q2.tsv")


#Para obtener una gamma, alpha, and beta diversity #Es importante quitar valores nulos 

d0p <- hill_taxa_parti(Df_batan_comm, q=0)
d1p <- hill_taxa_parti(Df_batan_comm, q=1)
write.table(x=d0p, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/gab_diversity_q0p.tsv")
write.table(x=d1p, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/gab_diversity_q1p.tsv")

#Mensaje que aparece si no eliminas valores nulos 
#Warning message:
#In hill_taxa_parti(Df_batan_comm, q = 0) :
#  Some species in comm data were not observed in any site,
#delete them...

#Un Ejemplo de plot para un archivo agrupado
qplot(q0 = d1, y = "", geom = "boxplot", col = I("darkblue"), fill = I("lightblue"), ylab = d1, main = "AlfaDiversity Batan")

qplot(q1 = d2, y = "", geom = "boxplot", col = I("")     )


qplot(q0 = Alfa_diversidad_q0, y = "", geom = "boxplot", col = I("darkblue"), fill = I("lightblue"), ylab = d1, main = "AlfaDiversity Batan")

###################################### BOXPLOTS ####################################
######################  ##############################  ######################################
##################################Análisis estadístico"#################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(devtools)
library(ggpubr)

Alfa_diversidad_q2 <- read.csv("div_q2.csv", fileEncoding = 'UTF-8-BOM')
View(Alfa_diversidad_q2)

##### El archivo que adjunto tiene la forma 
#   value tratamiento
#1    68     control
#2   351     control
#3    75     control
#4    27     residuo
#5    60     residuo
#6    63     residuo

Alfa_diversidad_q1$tratamiento <- factor(Alfa_diversidad_q1$tratamiento, levels = unique(Alfa_diversidad_q1$tratamiento))
levels(Alfa_diversidad_q1$tratamiento)
#Me debe salir 
#[1] "control" "residuo"

ggboxplot(Alfa_diversidad_q1, x = "tratamiento", y = "value", add = c("mean_se", "jitter"), color = "tratamiento", palette = c("#00AFBB", "#FC4E07"), order = c("control", "residuo"), ylab = "q=1", xlab = "tratamiento")


kruskal.test(value ~ tratamiento, data = Alfa_diversidad_q1)

#Me salió data:  value by tratamiento
#Kruskal-Wallis chi-squared = 12.402, df = 1, p-value =  0.0004289

#Como el valor p es menor que el nivel de significancia 0.05, podemos concluir que existen diferencias significativas entre los grupos de tratamiento.

# Mean plots
# ++++++++++++++++++++
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)

################# grafico de líneas

ggline(Alfa_diversidad_q1, x = "tratamiento", y = "value", 
       add = c("mean_se", "jitter"), color = "tratamiento", palette = c("#00AFBB", "#FC4E07"), 
       order = c("control", "residuo"),
       ylab = "q=1", xlab = "tratamiento")

#Como solo tengo dos conjuntos de datos que serán comparados esta prueba es suficiente para demostrar significancia estadística en el 
# efecto del tratamiento obserrvado sobre la diversidad de los suelos
#En caso de que sean más de 2 onjuntosde datos la Kruskal walis solo nis dirá qu hay diferencias significativas 
#Pero no nos dice entre cuales y ciales conjuntos de datos existe esa diferencia
#Por tano es necesario hacer una ptueba paeada una a una para identificar la significancia estadística entre grupos
#Aplicamos el Wilcox test
#Este  método de Wilcox por pares nos dice las significancia estadística en más de 2 grupos de información
#Me parece que también se puede aplicar la prueba de Tukey, pero esta información debe ser verificada


pairwise.wilcox.test(Alfa_diversidad_q1$value, Alfa_diversidad_q1$tratamiento,
                 p.adjust.method = "BH")

#Me salió
Pairwise comparisons using Wilcoxon rank sum exact test 

data:  Alfa_diversidad_q1$value and Alfa_diversidad_q1$tratamiento 

        control
residuo 0.00029

P value adjustment method: BH 

#########TukeyHSD(Alfa_diversidad_q1, tratamiento, ordered = FALSE, conf.level = 0.95) ##Verificar formula 

################################################################
#############################################
###PARA HACER EL ANÁLISIS DE LOS DIFERENTES TIPOS DE AGRICULTURA
##########################################################

Df_batan_comm <- read.csv("ASV_batan_CP-CA.csv", header = TRUE, row.names = 1, fileEncoding = 'UTF-8-BOM')
View(Df_batan_comm)

q0 <- hill_taxa(Df_batan_comm, q = 0, MARGIN = 1)
q1 <- hill_taxa(Df_batan_comm, q=1, MARGIN = 1)
q2 <- hill_taxa(Df_batan_comm, q=2, MARGIN = 1)
q3 <- hill_taxa(Df_batan_comm, q=3, MARGIN = 1)

write.table(x=q0, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_CP_CA_q0.tsv")
write.table(x=q1, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_CP_CA_q1.tsv")
write.table(x=q2, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_CP_CA_q2.tsv")
write.table(x=q3, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/diversity_CP_CA_q3.tsv")

#Para obtener una gamma, alpha, and beta diversity #Es importante quitar valores nulos 

d0p <- hill_taxa_parti(Df_batan_comm, q=0)
d1p <- hill_taxa_parti(Df_batan_comm, q=1)
write.table(x=d0p, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/gab_diversity_q0p.tsv")
write.table(x=d1p, quote = FALSE, sep = "\t", file = "~/MiProyecto/MiAnalisis/Taxonomy/MergeTabs/Level7/diversity/gab_diversity_q1p.tsv")

#####################################################################3#
#######################################################################3

Diversidad_q0_dia_cont_res <- read.csv("q0_CP-CA_residuo_control.csv", header = TRUE, fileEncoding = 'UTF-8-BOM')
View(Diversidad_q0_dia_cont_res)

Diversidad_q0_dia_cont_res$Agricultura <- factor(Diversidad_q0_dia_cont_res$Agricultura, levels = unique(Diversidad_q0_dia_cont_res$Agricultura))
levels(Diversidad_q0_dia_cont_res$Agricultura)

#Me debe salir 
#[1] "CP-R", "CP-C", "CA-R", "CA-C"

ggboxplot(Diversidad_q0_dia_cont_res, x = "Agricultura", y = "q0", add = c("mean_se", "jitter"), color = "Agricultura", palette = c("#00AFBB", "#000080", "#FC4E07", "#800000"), order = c("CP-R", "CP-C", "CA-R", "CA-C"), ylab = "q=0", xlab = "Agricultura")


kruskal.test(q0 ~ Agricultura, data = Diversidad_q0_dia_cont_res)

ggline(Diversidad_q0_dia_cont_res, x = "Agricultura", y = "q0", 
       add = c("mean_se", "jitter"), color = "Agricultura", palette = c("#00AFBB", "#000080", "#FC4E07", "#800000"), 
       order = c("CP-R", "CP-C", "CA-R", "CA-C"),
       ylab = "q=0", xlab = "Agricultura")


pairwise.wilcox.test(Diversidad_q0_dia_cont_res$q0, Diversidad_q0_dia_cont_res$Agricultura,
                     p.adjust.method = "BH")