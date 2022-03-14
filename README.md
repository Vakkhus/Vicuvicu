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
```R
library(biomod2)
library(raster)
library(sf)
library(tidyverse)
library(doParallel)
```
Para luego leer los datos de ocurrencia y crear un stack de las variables bioclimáticas en la zona de estudio: 

```R
data=read_sf('.../allvicugnas_nov2021.gpkg')              #ocurrencias
limite=read_sf('.../teow_vicugna1.shp')                   #límite zona estudio 
lista=list.files('.../bio/',pattern = '.tif')             #lista de archivos con variables bioclímáticas
bio=stack(lista)                                          #ráster stack con variables bioclimáticas
biof=stack(mask(crop(bio,st_zm(limite)),st_zm(limite)))   #recorte de stack a zona de estudio
```
Construímos el modelo de manera independiente para cada linaje `i`, por lo que, para el linaje seleccionado el set de variables bioclimáticas correspondiente: 
```R
if (i=="mensalis"){bio=biof[[c("bio3","bio7","bio9","bio10","bio13","bio15")]]}
if (i=="hybrid"){bio=biof[[c("bio4","bio6","bio7","bio9","bio12","bio15")]]}
if (i=="vicugna"){bio=biof[[c("bio2","bio4","bio7","bio12","bio18","bio19")]]}
```
Luego extraemos de la base de datos, los puntos corresponientes al linaje seleccionado y le damos [formato](https://www.rdocumentation.org/packages/biomod2/versions/3.5.1/topics/BIOMOD_FormatingData) a los datos. `PA.nb.rep` corresponde al número de veces que se seleccionan las pseudoausencias y `PA.nb.absences` el número de puntos de pseudoausencias seleccionadas en cada ocación (`5 * n`, donde n es el número de putnos de ocurrencia).

```R
myRespName <- i
    data_resp=as.data.frame(data[which(data$Lineage==myRespName),])[,1:6]
    data_resp$val=1
    myResp <- as.numeric(data_resp$val)
    myRespXY <- data_resp[,c("Longitude","Latitude")]
    myExpl = bio
    gc()
    
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
    expl.var = myExpl,
    resp.xy = myRespXY,
    resp.name = myRespName,
    PA.nb.rep=5,
    PA.nb.absences=5*length(data_resp$val),
    PA.strategy='random')
```

Según el linaje seleccionado, se selecciona un porcentaje ds` de datos para fraccionar los datos de entrenamiento, de modo que, para cada linaje el set de datos de entrenamiento sea aproximadamente del mismo tamaño (80 puntos), equivalente al 20% de los datos del linaje con menos observaciones. 
```R
if (i=="mensalis"){ds=17}
if (i=="hybrid"){ds=20}
if (i=="vicugna"){ds=10}
```

Luego procedemos a la [construcción del modelo](https://www.rdocumentation.org/packages/biomod2/versions/3.5.1/topics/BIOMOD_Modeling), seleccionando el número de repeticiones `NbRunEval` para cada set de pseudoausencias creado, las métricas de evaluación `models.eval.meth` y los modelos `models` a utilizar: 
```R
myBiomodModelOut <- BIOMOD_Modeling(
  data=myBiomodData,
  models = c("GBM","RF","MAXENT.Phillips"),
  NbRunEval=20,
  DataSplit=ds,
  VarImport = 5,
  models.eval.meth = c("KAPPA", "TSS", "ROC","ACCURACY"),
  SaveObj = TRUE,
  do.full.models = TRUE)
```

Una vez construídos los modelos individuales (Maxent, RF, GBM), procedemos a generar un [ensamble de modelos](https://www.rdocumentation.org/packages/biomod2/versions/3.5.1/topics/BIOMOD_EnsembleModeling):

```R
myBiomodEM <- BIOMOD_EnsembleModeling(myBiomodModelOut,
  chosen.models = 'all',
  prob.mean=TRUE,
  prob.cv =TRUE,
  VarImport=5)
```

y finalmente una [proyección](https://www.rdocumentation.org/packages/biomod2/versions/3.5.1/topics/BIOMOD_Projection) de la ocurrencia utilizando el ensamble generado, con las variables bioclimáticas del escenario actual. 
```R
myBiomodProjection <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'curr',
  compress = TRUE,
  build.clamping.mask = FALSE,
  keep.in.memory=FALSE)
  
BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection, EM.output = myBiomodEM)
```
Cada proyección, para los modelos por separado y ensamblados es guardada en la carpeta seleccionada como directorio de trabajo, así como los archivos complementarios y los archivos que contienen cada uno de los modelos. 

Por defecto es posible modelar la distribución de una especie (o linaje) a la vez. Es posible aprovechar de mejor manera cada hilo del procesador haciendo uso de la librería [`snowfall`](https://cran.r-project.org/web/packages/snowfall/snowfall.pdf), con la que podemos modelar en paralelo la distribución de varias especies. Para ello, es necesario crear una función que a partir de un parámetro `i`, que corresponde a cada especie a modelar, realice el proceso completo.

```R
MyBiomodSF <- function(i){
#Proceso de modelado
}

library(snowfall)
sfInit(parallel=TRUE, cpus=6) #cpus = número de hilos/núcleos disponibles

## Paquetes a utilizar: 
sfLibrary('biomod2', character.only=TRUE)
sfLibrary('raster', character.only=TRUE)
sfLibrary('sf', character.only=TRUE)

## Seleccionar el nombre de las variables y exportar variables a utilizar, para que la función interna pueda "verlas":
sp.names=c("hybrid","vicugna","mensalis")
sfExportAll() 

## Ejecutar en paralelo:
mySFModelsOut <- sfLapply(sp.names, MyBiomodSF)

## stop:
sfStop(nostop=FALSE)
```

## **Evaluación del modelo**

#### **Exactitud**
Para evaluar el modelo hacemos uso de la función `load` para cargar el modelo a analizar y luego mediante manejo de los datos obtenemos las métricas para cada réplica, modelo o ensamble: 

```R
library(dplyr)
#setear el working directory en la carpeta donde estén las carpetas "mensalis", "vicugna" e "hybrid".
setwd(".../bio")

#se carga el modelo a evaluar
load('.../bio/vicugna/vicugna.1640727040.models.out')
#se carga el ensamble a evaluar
vicugna_en=load('.../bio/vicugna/vicugna.1640727040ensemble.models.out')
model_en <- get(vicugna_en)
```

utilizando @ en la variable de cada modelo, se exploran los diferentes ouputs, en este caso, las métricas exactitud. En el primer caso, la variable almacena los resultados de cada corrida o `run` para cada réplica de muestreo de pseudoausencias `PA`, además de un modelo `Full` que ensambla todas las réplicas para un mismo `PA`. Luego, se genera un data frame con las métricas de evaluación para cada modelo `Full` por separado, y se promedian: 

```R
#evaluación de los modelos por separado:
w=as_tibble(vicugna.1640727040.models.out@models.evaluation@val,rownames='names')

#Maxent:
test_eval_maxent=w %>% select(starts_with('Testing.data.MAXENT.Phillips.Full')) %>% add_column(stat=rownames(vicugna.1640727040.models.out@models.evaluation@val)) %>% mutate(mean = rowMeans(across(starts_with("Testing")))) %>% select(tail(names(.), 2))

#GBM:
test_eval_gbm=w %>% select(starts_with('Testing.data.GBM.Full')) %>% add_column(stat=rownames(vicugna.1640727040.models.out@models.evaluation@val)) %>% mutate(mean = rowMeans(across(starts_with("Testing")))) %>% select(tail(names(.), 2))

#Random Forest
test_eval_rf=w %>% select(starts_with('Testing.data.RF.Full')) %>% add_column(stat=rownames(vicugna.1640727040.models.out@models.evaluation@val)) %>% mutate(mean = rowMeans(across(starts_with("Testing")))) %>% select(tail(names(.), 2))
```
| vicugna |
| -------------  |
| Modelo | Kappa | TSS | ROC | Accuracy |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Maxent  | 0.844 | 0.876 | 0.974 | 0.958 |
| GBM  | 0.879 | 0.941 | 0.992 | 0.970 |
| Random Forest | 0.976 | 0.989 | 1.0 | 0.993 |

Mientras que en el segundo caso, la variable almacena las métricas del ensamble de modelos (Maxent, RF, GBM) para cada réplica de muestreo de pseudoausencias `PA`. Luego, se calcula el promedio de los indicadores: 

```R
#evaluación de los ensambles de modelos:

w=as_tibble(get_evaluations(model_en),rownames='names') %>% select(starts_with('vicugna_EMmeanByACCURACY_mergedAlgo_Full'))
test_eval_ensemble=data.frame(PA1=c(KAPPA=0,TSS=0,ROC=0))

for (i in 1:5){
  a=w[[i]][,1]
  test_eval_ensemble[paste('PA',i,sep='')]=a
}

test_eval_ensemble['mean']=rowMeans(test_eval_ensemble)
```

| eval | mean vicugna |  mean hybrid | mean mensalis  | 
| ------------- | ------------- | ------------- | ------------- |
|KAPPA | 0.9386| | |
|TSS   | 0.9608| | |
|ROC   | 0.9980| | |

#### **Importancia de variables**

De manera similar, podemos calcular la importancia de las variables para los modelos por separado, promediando las importancias de los modelos `Full`, o para el ensamble de modelos, promediando las importancia del mensamble para cada `PA`:

```R
#importancia de variables en modelos por separado

w=as_tibble(vicugna.1640727040.models.out@variables.importances@val,rownames='names')

var_imp_maxent=w %>% select(starts_with('MAXENT.Phillips.Full')) %>% add_column(var=rownames(vicugna.1640727040.models.out@variables.importances@val),.before = 1) %>% mutate(mean = rowMeans(across(starts_with("MAXENT")))) 

var_imp_gbm=w %>% select(starts_with('GBM.Full')) %>% add_column(var=rownames(vicugna.1640727040.models.out@variables.importances@val),.before = 1) %>% mutate(mean = rowMeans(across(starts_with("GBM")))) 

var_imp_rf=w %>% select(starts_with('RF.Full')) %>% add_column(var=rownames(vicugna.1640727040.models.out@variables.importances@val),.before = 1) %>% mutate(mean = rowMeans(across(starts_with("RF")))) 
```
var

#importancia de variables en ensambles

w=as_tibble(get_variables_importance(model_en),rownames='names') %>% select(contains(c('vicugna_EMcvByACCURACY_mergedAlgo_Full','names')))
varimp_ensemble=data.frame(PA1=rep(0,6));rownames(varimp_ensemble)=w$names

for (i in 1:5){
  o=seq(5*(i-1)+1,i*5)
  a=mutate(w[o], mean = rowMeans(across(starts_with("rand"))))[6]
  varimp_ensemble[paste('PA',i,sep='')]=a
}

varimp_ensemble['mean']=rowMeans(varimp_ensemble)
