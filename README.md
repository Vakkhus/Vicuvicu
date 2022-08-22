Se contruyó un modelo de probabilidad de ocurrencia para _Vicugna vicugna_ y sus tres linajes (_vicugna_, _mensalis_, híbrido) utilizando un ensamble de los algoritmos [Maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/), [Random Forest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf) y [Generalized Boosting Model (gbm)](https://cran.r-project.org/web/packages/gbm/gbm.pdf), usualmente llamado Boosted Regression Trees utilizando el paquete [biomod2](https://cran.r-project.org/web/packages/biomod2/biomod2.pdf) en [R](https://www.r-project.org/). 

biomod2 ofrece la posibilidad de ejecutar 10 técnicas de modelado de última generación para describir y modelar las relaciones entre una especie determinada y su entorno. Es un intento de definir el nicho ecológico de una especie en particular usando variables ambientales (temperatura, precipitación, ...) con el uso potencial de hacer, por ejemplo, proyecciones futuras bajo escenarios de cambio de uso de suelo y clima. Aunque se ha desarrollado principalmente para ecólogos que tienen como objetivo predecir la distribución de especies, biomod2 también se puede utilizar para modelar cualquier dato binomial (por ejemplo, gen, marcadores, ecosistema...) en función de cualquier variable explicativa.

## **Variable respuesta**
La variable respuesta está conformada a 1657 puntos de presencia de la especie, donde 792 corresponden al linaje _vicugna_, 465 al linaje _mensalis_ y 400 al linaje híbrido, obtenidos de .... Además, un set de  `5 * n` puntos de background, donde n es el número de puntos utilizados para la construcción de la instancia del modelo, es generado en cada instancia y distribuídos aleatóreamente dentro de la zona de estudio.

## **Variables predictoras**
Se utilizaron las [variables bioclimáticas](https://www.worldclim.org/data/bioclim.html) disponibles en [CHELSA](https://chelsa-climate.org/) como variables predictoras y posteriormente, mediante un [análisis de componentes principales (PCA)](https://www.sciencedirect.com/science/article/abs/pii/S0167947304002014?via%3Dihub) se determinó el set de variables más relevante para cada uno de los linajes. 

| Linaje | Variables bioclimáticas seleccionadas |
| ------------- | ------------- |
| mensalis  | bio3, bio7, bio9, bio10, bio13, bio15  |
| hybrid  | bio4, bio6, bio7, bio9, bio12, bio15 |
| vicugna | bio2, bio4, bio7, bio12, bio18, bio19 |

## **Área de estudio**


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

## **Exactitud**
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
#### vicugna
| Modelo | Kappa | TSS | ROC | Accuracy |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Maxent  | 0.844 | 0.876 | 0.974 | 0.958 |
| GBM  | 0.879 | 0.941 | 0.992 | 0.970 |
| Random Forest | 0.976 | 0.989 | 1.0 | 0.993 |

#### hybrid
| Modelo | Kappa | TSS | ROC | Accuracy |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Maxent  | 0.940 | 0.919 |0.967  |0.984  |
| GBM  | 0.972 |0.984 | 0.999 | 0.992 |
| Random Forest |0.987  |0.995  |1.0 |0.996  |

#### mensalis
| Modelo | Kappa | TSS | ROC | Accuracy |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Maxent  | 0.898 | 0.852 | 0.936 | 0.973  |
| GBM  | 0.919|0.952  |0.997  |0.977 | 
| Random Forest | 0.988  | 0.994 |1.0 | 0.996 |


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

| eval | vicugna |  hybrid | mensalis  | 
| ------------- | ------------- | ------------- | ------------- |
|KAPPA | 0.9386|0.5432 |0.7810 |
|TSS   | 0.9608|0.5478 |0.7734 |
|ROC   | 0.9980|0.8062 |0.8820 |

## **Importancia de variables**

De manera similar, podemos calcular la importancia de las variables para los modelos por separado, promediando las importancias de los modelos `Full`, o para el ensamble de modelos, promediando las importancia del mensamble para cada `PA`:

```R
#importancia de variables en modelos por separado

w=as_tibble(vicugna.1640727040.models.out@variables.importances@val,rownames='names')

var_imp_maxent=w %>% select(starts_with('MAXENT.Phillips.Full')) %>% add_column(var=rownames(vicugna.1640727040.models.out@variables.importances@val),.before = 1) %>% mutate(mean = rowMeans(across(starts_with("MAXENT")))) 

var_imp_gbm=w %>% select(starts_with('GBM.Full')) %>% add_column(var=rownames(vicugna.1640727040.models.out@variables.importances@val),.before = 1) %>% mutate(mean = rowMeans(across(starts_with("GBM")))) 

var_imp_rf=w %>% select(starts_with('RF.Full')) %>% add_column(var=rownames(vicugna.1640727040.models.out@variables.importances@val),.before = 1) %>% mutate(mean = rowMeans(across(starts_with("RF")))) 
```
#### vicugna
| var | Maxent |  GBM | Random Forest | 
| ------------- | ------------- | ------------- | ------------- |
|bio2  |0.4276 |0.5216  |0.2750|
|bio4  |0.2208 |0.1098  |0.1536|
|bio7  |0.1050 |0.0364  |0.1254|
|bio12 |0.5568 |0.2744  |0.3468|
|bio18 |0.0942 |0.0032  |0.0528|
|bio19 |0.0546 |0.0014  |0.0332|


#### hybrid
| var | Maxent |  GBM | Random Forest | 
| ------------- | ------------- | ------------- | ------------- |
|bio4 |0.0850 |0.0166 |0.0632 |
|bio6  |0.1832 |0.0112|0.0610 |
|bio7  |0.0876 |0.0168 |0.0474 |
|bio9 |0.2214 |0.1492 |0.0418 |
|bio12 |0.1502 |0.0068 |0.0576 |
|bio15 |0.1096 |0.6150 |0.1840 |

#### mensalis
| var | Maxent |  GBM | Random Forest | 
| ------------- | ------------- | ------------- | ------------- |
|bio3  |0.2006|0.0392|0.1460|
|bio7  |0.3048|0.1732|0.1120|
|bio9  |0.0228|0.0014|0.0600|
|bio10 |0.1780|0.5674|0.1080|
|bio13 |0.1328|0.0856|0.0494|
|bio15 |0.0804|0.0020|0.0688|

```R
#importancia de variables en ensambles

w=as_tibble(get_variables_importance(model_en),rownames='names') %>% select(contains(c('vicugna_EMcvByACCURACY_mergedAlgo_Full','names')))
varimp_ensemble=data.frame(PA1=rep(0,6));rownames(varimp_ensemble)=w$names

for (i in 1:5){
  o=seq(5*(i-1)+1,i*5)
  a=mutate(w[o], mean = rowMeans(across(starts_with("rand"))))[6]
  varimp_ensemble[paste('PA',i,sep='')]=a
}

varimp_ensemble['mean']=rowMeans(varimp_ensemble)
```
| var | vicugna|  hybrid | mensalis | 
| ------------- | ------------- | ------------- | ------------- |
|bio2  |0.4179 |-      |-      |
|bio3  |-      |-      | 0.2076|
|bio4  |0.1901 |0.1281 |-      |
|bio6  |-      |0.2368 |-      |
|bio7  |0.2251 |0.2089 |0.2194 |
|bio9  |-      |0.1395 |0.1367 |
|bio10 |-      |-      |0.4474 |
|bio12 |0.3482 |0.0722 |-      |
|bio13  |-     |-      |0.0993 |
|bio15  |-     |0.3538 |0.0895 |
|bio18 |0.1131 |-      |-      |
|bio19 |0.0737 |-      |-      |

## Curvas de respuesta

Para obtener las curvas de respuesta por modelo (cada línea representa un modelo `Full` y cada color un algoritmo): 

```R
setwd(".../bio")
#se carga el modelo a evaluar
load('.../bio/vicugna/vicugna.1640727040.models.out')

vicugna_mod=load('.../bio/vicugna/vicugna.1640727040.models.out')
my_model <- get(vicugna_mod)

model=BIOMOD_LoadModels(my_model, models = c('GBM','RF','MAXENT.Phillips'))
model=model[seq(21,315,21)]#se muestrans los modelos "full"

myRespPlot2D <- 
  response.plot2(models = model, Data = get_formal_data(my_model, 'expl.var'),show.variables = c("bio2",  "bio4",  "bio7",  "bio12", "bio18", "bio19"),fixed.var.metric = 'mean',
    col = c(rep("blue",5), rep("red",5),rep("green",5)),legend = TRUE)
```
#### vicugna
![enter image description here](https://github.com/Vakkhus/Vicuvicu/blob/main/Plots/Response_model_vicugna.png?raw=true)
#### hybrid
![enter image description here](https://github.com/Vakkhus/Vicuvicu/blob/main/Plots/Response_model_hybrid.png?raw=true)
#### mensalis
![enter image description here](https://github.com/Vakkhus/Vicuvicu/blob/main/Plots/Response_model_mensalis.png?raw=true)


Para obtener las curvas de respuesta por ensamble (cada línea y color representa una réplica de `PA` para el ensamble): 
```R
#se carga el modelo a evaluar
setwd("S:/UACh/vicuñas/chelsa_historico/bio")
load('S:/UACh/vicuñas/chelsa_historico/bio/vicugna/vicugna.1640727040ensemble.models.out')

#se carga el archivo con los ensambles
vicugna_mod=load('S:/UACh/vicuñas/chelsa_historico/bio/vicugna/vicugna.1640727040ensemble.models.out')
my_model <- get(vicugna_mod)
#se carga el archivo con los modelos por separado (solo se usa para obtener los valores de las variables)
vicugna_mod2=load('S:/UACh/vicuñas/chelsa_historico/bio/vicugna/vicugna.1640727040.models.out')
my_model2 <- get(vicugna_mod2)

#se seleccionan los ensambles "full" (combinación de corridas), promediados por Accuracy
model=c("vicugna_EMmeanByACCURACY_mergedAlgo_Full_PA1","vicugna_EMmeanByACCURACY_mergedAlgo_Full_PA2","vicugna_EMmeanByACCURACY_mergedAlgo_Full_PA3","vicugna_EMmeanByACCURACY_mergedAlgo_Full_PA4","vicugna_EMmeanByACCURACY_mergedAlgo_Full_PA5")
#para ver todos los ensambles disponibles: 
BIOMOD_LoadModels(my_model, models = c('PA1','PA2','PA3','PA4','PA5'))

#se grafican las curvas de respuesta
myRespPlot2D <- response.plot2(models = model, Data = get_formal_data(my_model2, 'expl.var'),show.variables = c("bio2",  "bio4",  "bio7",  "bio12", "bio18", "bio19"),fixed.var.metric = 'mean', col = c('#d9ed92', "#99d98c","#52b69a","#1a759f","#184e77"),legend = TRUE, save.file = 'pdf')
```
#### vicugna
![enter image description here](https://github.com/Vakkhus/Vicuvicu/blob/main/Plots/Response_ensemble_vicugna.png?raw=true)
#### hybrid
![enter image description here](https://github.com/Vakkhus/Vicuvicu/blob/main/Plots/Response_ensemble_hybrid.png?raw=true)
#### mensalis
![enter image description here](https://github.com/Vakkhus/Vicuvicu/blob/main/Plots/Response_ensemble_mensalis.png?raw=true)


## **Proyecciones**

#### mensalis 
<img width="500" alt="image" src="https://user-images.githubusercontent.com/43461660/159535805-c311ac4d-dc35-4446-8b36-e12a6fc2be79.png">

#### vicugna 
<img width="500" alt="image" src="https://user-images.githubusercontent.com/43461660/159535848-6c69a51b-70f6-43c7-baa5-deb8a6b18b07.png">

#### hybrid
<img width="500" alt="image" src="https://user-images.githubusercontent.com/43461660/159535881-543b8e87-1a15-48b7-9660-127745f49c95.png">

#### all lineages
<img width="500" alt="image" src="https://user-images.githubusercontent.com/43461660/159535938-bb0b06b8-f6d9-4bbb-99d1-e27b6945ef79.png">


# **Modelos con todas las variables bioclimáticas**

A raíz de la poca respuesta de la probabilidad de ocurrencia a las variables en las curvas, se construyó un modelo nuevo para cada linaje y que utiliza todas las variables bioclimáticas, para verificar si se están utilizando aquellas variables que efectivamente son relevantes para la especie. 

Las importancia de las variables para cada linaje es:

|       | vicugna     | mensalis    | hybrid      |
|-------|-------------|-------------|-------------|
| bio2  | *0.28074092 | 0.12229212  | *0.16986900 |
| bio6  | *0.22781432 | *0.17068072 | *0.14219548 |
| bio11 | *0.20332840 | 0.15484432  | 0.06883244  |
| bio1  | *0.17025972 | *0.17939124 | 0.04697016  |
| bio9  | *0.16242868 | *0.15700300 | 0.05752312  |
| bio7  | *0.15459416 | 0.09498436  | 0.04377308  |
| bio12 | 0.13059460  | 0.03437792  | *0.09249912 |
| bio16 | 0.11143008  | 0.05764960  | 0.04409472  |
| bio10 | 0.11005576  | *0.26365488 | 0.04860852  |
| bio4  | 0.10580292  | 0.09568504  | 0.05743544  |
| bio3  | 0.08931480  | 0.08082560  | 0.04712856  |
| bio5  | 0.06929360  | *0.18473520 | 0.04315840  |
| bio8  | 0.05987564  | *0.16932212 | 0.02510676  |
| bio13 | 0.04628144  | 0.05183364  | 0.04068596  |
| bio15 | 0.04247180  | 0.02703372  | *0.21201004 |
| bio18 | 0.04007872  | 0.05483144  | 0.02956224  |
| bio19 | 0.03891472  | 0.02570144  | *0.09635752 |
| bio17 | 0.03294220  | 0.02247824  | *0.11271572 |
| bio14 | 0.02865220  | 0.01588932  | 0.02812072  |

y las 6 variables más importantes son: 
| mensalis   PCA | mensalis FULL | \| | hybrid PCA | hybrid FULL | \| | vicugna PCA | vicugna FULL |
|----------------|---------------|----|------------|-------------|----|-------------|--------------|
| bio3           | bio1          | \| | bio4       | bio2        | \| | bio2        | bio2         |
| bio7           | bio6          | \| | bio6       | bio6        | \| | bio4        | bio6         |
| bio9           | bio9          | \| | bio7       | bio12       | \| | bio7        | bio11        |
| bio10          | bio10         | \| | bio9       | bio15       | \| | bio12       | bio1         |
| bio13          | bio5          | \| | bio12      | bio19       | \| | bio18       | bio9         |
| bio15          | bio8          | \| | bio15      | bio17       | \| | bio19       | bio7         |

Sin embargo, aún con las variables más importantes identificadas por el modelo, no hay cambios importantes en las curvas de respuesta, sobre todo para MAXENT.

<img width="568" alt="image" src="https://user-images.githubusercontent.com/43461660/159533211-38e2ca32-df67-40d2-93fe-6ff83bc4b235.png">

Si graficamos el modelo anterior para el linaje mensalis, pero esta vez para cada PA, cada modelo y cada réplica, el restulado es el siguiente: 

<img width="569" alt="image" src="https://user-images.githubusercontent.com/43461660/159518307-5964844c-a812-478d-b330-00b6c72390de.png">

# **Modelos con variables climáticas relevantes para la especie (08-2022)**

A partir de un análisis de componentes principales se estableció un set de variables relevantes para la especie, y se contruyeron modelos de distribución usando la metodología anterior e incorporando todas las variables de este set. Las variables seleccionadas fueron: "bio1", "bio2", "bio3", "bio4",  "bio6", "bio7", "bio8", "bio10", "bio12", "bio14", "bio15" y"bio19".

La importancia media de cada una de las variables entre los 5 `PA` se detalla a continuación:

|       | vicugna     | mensalis    | hybrid      |all lineages     |
|-------|-------------|-------------|-------------|-------------|
| bio1  |0.229|0.190|0.058 |0.277|
| bio2  |0.283|0.170|0.180|0.255|
| bio3  |0.112|0.113|0.054|0.119|
| bio4  |0.160|0.126|0.061|0.119|
| bio6  |0.292|0.185|0.142|0.257|
| bio7  |0.185|0.114|0.053|0.124|
| bio8  |0.117|0.158|0.043|0.128|
| bio10 |0.198|0.434|0.058|0.253|
| bio12 |0.125|0.068|0.042|0.118|
| bio14 |0.043|0.038|0.032|0.044|
| bio15 |0.055|0.038|0.250|0.057|
| bio19 |0.055|0.100|0.092|0.092|


![importance_bar](https://user-images.githubusercontent.com/43461660/185834385-7c1ad769-8ca2-45b9-9049-32e3de998357.png)

## Curvas de respuesta

### all lineages
#### by model
![response_model_all](https://user-images.githubusercontent.com/43461660/185834633-18bbe9d3-6066-40cc-9a9c-4b29773aa9ae.png)
#### by ensemble
![response_ensemble_all](https://user-images.githubusercontent.com/43461660/185834642-f681b8a6-6e83-4cd6-be0b-d583d35e46ce.png)

### vicugna
#### by model
![response_model_vicugna](https://user-images.githubusercontent.com/43461660/185834664-f83684bc-6210-4e87-aa26-e433f8ca8d38.png)
#### byensemble
![response_ensemble_vicugna](https://user-images.githubusercontent.com/43461660/185834669-e0975c8e-bb1b-46d4-b83f-44f0571ab9bb.png)

### mensalis
#### by model
![response_model_mensalis](https://user-images.githubusercontent.com/43461660/185834711-728805bb-bf64-4f88-b4b6-dd7eb305f35d.png)
#### by ensemble
![response_ensemble_mensalis](https://user-images.githubusercontent.com/43461660/185834722-7173afb1-3a2e-45ed-90f6-8b9afa10028e.png)

### hybrid
#### by model
![response_model_hybrid](https://user-images.githubusercontent.com/43461660/185834761-80e13480-0978-4ef9-ab7e-4af1365bddf3.png)
#### ensemble
![response_ensemble_hybrid](https://user-images.githubusercontent.com/43461660/185834781-e0b78ee3-4700-4045-8c02-5d851eee10f5.png)
