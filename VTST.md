# VTSTS.py

## Descrición

`VTSTS.py` é unha ferramenta desenvolvida para o tratamento de perfís VTST obtidos a partir de escaneos relaxados realizados con Gaussian. O script automatiza o cálculo de propiedades termodinámicas e cinéticas ao longo dunha coordenada de reacción e facilita a análise de procesos inicialmente clasificados como reaccións sen barreira.

## Funcionalidades principais

### Procesamento automático de escaneos

O programa identifica e procesa automaticamente todos os ficheiros do tipo:

```text
freq_1.out
freq_2.out
freq_3.out
...
```

presentes no directorio especificado polo usuario.

Isto permite analizar de forma eficiente escaneos compostos por un gran número de puntos.

### Selección de rexións específicas do escaneo

O usuario pode limitar a análise a un intervalo concreto mediante:

```python
POINT_RANGE
```

ou excluír puntos individuais mediante:

```python
POINTS_EXCLUDE
```

o que resulta útil cando determinados puntos presentan problemas de converxencia ou comportamentos anómalos.

### Cálculo de propiedades termodinámicas

Para cada punto do escaneo calcúlanse:

* funcións de partición;
* enerxías libres de Gibbs;
* perfís de enerxía libre dependentes da temperatura.

Os cálculos realízanse empregando as rutinas implementadas en Pilgrim.

### Estimación de constantes cinéticas

A partir dos perfís de enerxía libre obtéñense constantes cinéticas dependentes da temperatura mediante a aproximación GFEA (*Gibbs Free Energy Approximation*).

### Corrección automática de frecuencias baixas

O script permite substituír automaticamente frecuencias vibracionais moi baixas por un valor definido polo usuario:

```python
FREQ_FIX_FROM_POINT
FREQ_FIX_VALUE_CM
```

Esta funcionalidade resulta especialmente útil en rexións da coordenada de reacción onde aparecen modos vibracionais moi brandos que poden introducir erros numéricos nas funcións de partición.

### Seguimento da coordenada de reacción

É posible monitorizar automaticamente a evolución dunha distancia interatómica determinada mediante:

```python
SCAN_ATOMS
```

permitindo relacionar directamente:

```text
Distancia → ΔG → k(T)
```

ao longo do escaneo.

### Xeración automática de resultados

Os resultados expórtanse automaticamente a un ficheiro Excel que inclúe:

* táboas completas de resultados;
* perfís de enerxía libre;
* gráficos ΔG fronte á distancia;
* valores de constantes cinéticas para todas as temperaturas estudadas.

## Aplicación neste traballo

Esta ferramenta empregouse para estudar reaccións inicialmente clasificadas como procesos sen barreira durante a exploración automática de mecanismos. Os resultados permitiron identificar posibles máximos variacionais ao longo da coordenada de reacción, obter perfís VTST dependentes da temperatura e estimar constantes cinéticas para a súa incorporación aos modelos cinéticos desenvolvidos neste traballo.
