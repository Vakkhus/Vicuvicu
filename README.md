Se contruyó un modelo de probabilidad de ocurrencia para Vicugna vicugna y sus tres linajes (vicugna, mensalis, híbrido) utilizando un ensamble de los algoritmos [Maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/), [Random Forest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf) y [Generalized Boosting Model (gbm)](https://cran.r-project.org/web/packages/gbm/gbm.pdf), usualmente llamado Boosted Regression Trees utilizando el paquete [biomod2](https://cran.r-project.org/web/packages/biomod2/biomod2.pdf) en [R](https://www.r-project.org/). 

biomod2 ofrece la posibilidad de ejecutar 10 técnicas de modelado de última generación para describir y modelar las relaciones entre una especie determinada y su entorno. Es un intento de definir el nicho ecológico de una especie en particular usando variables ambientales (temperatura, precipitación, ...) con el uso potencial de hacer, por ejemplo, proyecciones futuras bajo escenarios de cambio de uso de suelo y clima. Aunque se ha desarrollado principalmente para ecólogos que tienen como objetivo predecir la distribución de especies, biomod2 también se puede utilizar para modelar cualquier dato binomial (por ejemplo, gen, marcadores, ecosistema...) en función de cualquier variable explicativa.

**Variable respuesta**
La variable
