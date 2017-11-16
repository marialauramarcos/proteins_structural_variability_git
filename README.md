### VariabilidadEstructuralProteica ###

# Descripción
En este respositorio se encuentran los directorios, programas y archivos necesarios para analizar la divergencia estructural proteica, comparando datos experimentales obtenidos de la base de alineamientos estructurales múltiples Homstrad y datos teóricos obtenidos simulando mutaciones con el modelo "Linearly Forced - Elastic Network Model" (LF-ENM).
Para correr los programas debe completarse el archivo input ("input_MainProgram.csv"). Luego, pueden ralizarse distintos análisis:

## Comparación proteínas experimentales con mutantes mútiples teóricas cuyas mutaciones fueron seleccionadas según el stress model:
Debe correrse el programa "MainProgram.R", que realiza los siguientes pasos:

1) Lee el input.
2) Analiza la familia de proteínas ingresada y genera archivos con la información extraída. Si estos archivos ya existen puede evitarse re-analizar la familia con una de las opciones del input.
3) Genera mutantes múltiples de la proteína de referencia de la familia usando el modelo LF-ENM. Para seleccionar a cada mutación puntual se utiliza el modelo "stress model" bajo diferentes regimenes de selección: nula, baja, media o alta. Se generan tantas mutaciones como corresponde al <% id> con la proteína de referencia. Las mutantes generadas se guardan en archivos. Si estos archivos ya existen puede evitarse volver a generarlos con una de las opciones del input.
4) Analiza proteínas teóricas y experimentales y calcula medidas de variabilidad estructural como dr.squarei y Pn.

Para generar reportes de las familias ingresadas, luego de correr ambos pogramas, debe correrse el programa "MainReport.R"", cuyo input es "input_MainReport.csv". Los archivos necesarios para correr este programa ya se encuentran en los directorios correspondientes.
Los resultados de la corrida se verán en "FIGURES_REPORTS/".

## Comparación de mutaciones teóricas en cada sitio de la proteína:
Debe correrse "GenerateMutantsAllSites". Este programa genera mutaciones en cada sitio de la proteína y compara los resultados.
Los resultados de la corrida se verán en "OUT/out_mutants_all_sites/". Todas las correlaciones que se muestran aquí son de Spearman.

## Comparación de mutaciones teóricas cuyas mutaciones fueron aceptadas según el score experimental o su inverso:
Debe correrse "GenerateMutantsScore". Este programa genera mutantes múltiples de la proteína de referencia, seleccionado a cada una según el score experimental obtenido de ConsurfDB, o su iverso. Se generan tantas mutaciones como corresponde al <% id> con la proteína de referencia. 
Los resultados de la corrida se verán en "OUT/out_mutants_scores/". Todas las correlaciones que se muestran aquí son de Spearman.

# Librerías necesarias
Las librerías necesarias para correr los programas son las siguientes:

-bio3d
-seqinr
-ggplot2
-plyr
-reshape2
-quantreg
-png

Para descargar e instalar las librerías usar el comando "install.packages()".

