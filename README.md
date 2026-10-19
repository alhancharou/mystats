# mystats

R package with statistical functions and probability distributions.

Paquete con funciones en **español** para tus alumnos.

## Instalación local (rápida)

1. Descomprime el `.zip`.
2. En R ejecuta:

```r
install.packages("devtools")        # si hace falta
devtools::install_local("misalumnos")
library(misalumnos)
```

## Funciones incluidas

- `media(x)`, `varianza(x)`, `desviacion(x)`
- `cuantiles(x, probs=c(.25,.5,.75))`
- `normalizar(x)`, `estandarizar(x)`
- `correlacion(x, y, metodo="pearson")`
- `regresion_lineal(x, y)`
- `grafica(x, y=NULL, tipo="p")`
- `grafica_residuos(modelo)`
- `datos_ejemplo(n=50, seed=NULL)`

## Ejemplo rápido

```r
library(misalumnos)
set.seed(123)
df <- datos_ejemplo(30)
media(df$y)
m <- regresion_lineal(df$x, df$y)
summary(m)
grafica(df$x, df$y)
grafica_residuos(m)
```