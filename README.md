Se contruyó un modelo de probabilidad de ocurrencia para _Vicugna vicugna_ y sus tres linajes (_vicugna_, _mensalis_, híbrido) utilizando un ensamble de los algoritmos [Maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/), [Random Forest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf) y [Generalized Boosting Model (gbm)](https://cran.r-project.org/web/packages/gbm/gbm.pdf), usualmente llamado Boosted Regression Trees utilizando el paquete [biomod2](https://cran.r-project.org/web/packages/biomod2/biomod2.pdf) en [R](https://www.r-project.org/). 

biomod2 ofrece la posibilidad de ejecutar 10 técnicas de modelado de última generación para describir y modelar las relaciones entre una especie determinada y su entorno. Es un intento de definir el nicho ecológico de una especie en particular usando variables ambientales (temperatura, precipitación, ...) con el uso potencial de hacer, por ejemplo, proyecciones futuras bajo escenarios de cambio de uso de suelo y clima. Aunque se ha desarrollado principalmente para ecólogos que tienen como objetivo predecir la distribución de especies, biomod2 también se puede utilizar para modelar cualquier dato binomial (por ejemplo, gen, marcadores, ecosistema...) en función de cualquier variable explicativa.

## **Variable respuesta**
La variable respuesta está conformada a 1657 puntos de presencia de la especie, donde 792 corresponden al linaje _vicugna_, 465 al linaje _mensalis_ y 400 al linaje híbrido, obtenidos de .... Además, un set de  `5 * n` puntos de background, donde n es el número de puntos utilizados para la construcción de la instancia del modelo, es generado en cada instancia y distribuídos aleatóreamente dentro de la zona de estudio.

## **Variable predictoras**
Se utilizaron las [variables bioclimáticas](https://www.worldclim.org/data/bioclim.html) disponibles en [CHELSA](https://chelsa-climate.org/) como variables predictoras y posteriormente, mediante un [análisis de componentes principales (PCA)](https://www.sciencedirect.com/science/article/abs/pii/S0167947304002014?via%3Dihub) se determinó el set de variables más relevante para cada uno de los linajes. 

| Linaje | Variables bioclimáticas seleccionadas |
| ------------- | ------------- |
| mensalis  | bio3, bio7, bio9, bio10, bio13, bio15  |
| hybrid  | bio4, bio6, bio7, bio9, bio12, bio15 |
| vicugna | bio2, bio4, bio7, bio12, bio18, bio19 |

## **Construcción del modelo**

Comenzamos importando los paquetes a utilizar:
```
library(biomod2)
library(raster)
library(sf)
library(tidyverse)
library(doParallel)
```
Para luego leer los datos de ocurrencia y crear un stack de las variables bioclimáticas en la zona de estudio: 

```
data=read_sf('.../allvicugnas_nov2021.gpkg')              #ocurrencias
limite=read_sf('.../teow_vicugna1.shp')                   #límite zona estudio 
lista=list.files('.../bio/',pattern = '.tif')             #lista de archivos con variables bioclímáticas
bio=stack(lista)                                          #ráster stack con variables bioclimáticas
biof=stack(mask(crop(bio,st_zm(limite)),st_zm(limite)))   #recorte de stack a zona de estudio
```
