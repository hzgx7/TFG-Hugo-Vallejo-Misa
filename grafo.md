# grafo.py

## Descrición

`grafo.py` é unha ferramenta desenvolvida para a análise automática de redes de reacción xeradas mediante AutoMeKin e Pilgrim. O seu obxectivo principal é identificar rutas mecanisticamente relevantes dentro de redes reaccionais complexas e preparar automaticamente os ficheiros necesarios para posteriores simulacións cinéticas.

## Funcionalidades principais

### Construción automática do grafo reaccional

O script le os ficheiros:

* `RXNet.cg`
* `RXNet.barrless`

e constrúe internamente un grafo dirixido formado por mínimos (`MIN`), estados de transición (`TS`) e produtos (`PR`).

### Integración de reaccións con e sen barreira

As reaccións con barreira e sen barreira intégranse nunha única rede, permitindo analizar simultaneamente todos os procesos accesibles.

As reaccións sen barreira reciben identificadores artificiais mediante o parámetro:

```python
TS_BARRIERLESS_OFFSET = 800
```

o que facilita o seu tratamento como conexións equivalentes a estados de transición durante a análise da rede.

### Procura automática de camiños reaccionais

A partir dunha especie inicial definida polo usuario:

```python
start_path = ["MIN1"]
```

o script identifica automaticamente todos os camiños posibles cara aos produtos finais mediante un algoritmo de procura en profundidade (*Depth First Search*, DFS).

### Filtrado mecanístico e enerxético

A análise pode restrinxirse mediante:

```python
max_length
```

que limita o número máximo de mínimos intermedios presentes nun camiño,

e

```python
E_cutoff
```

que elimina rutas que superen unha enerxía máxima especificada polo usuario.

### Eliminación de redundancias

Cando existen varias conexións entre dúas especies, o script conserva automaticamente a ruta asociada ao estado de transición de menor enerxía, evitando duplicidades innecesarias.

### Identificación de pescozos de botella

Para cada camiño determínase a enerxía máxima alcanzada ao longo da ruta e identifícase automaticamente o camiño cun menor pescozo de botella enerxético.

Esta funcionalidade resulta especialmente útil para identificar as rutas potencialmente máis relevantes desde o punto de vista cinético.

### Estatísticas da rede

O programa xera información resumida sobre:

* número total de camiños identificados;
* número de reaccións consideradas;
* número de procesos sen barreira;
* distribución de camiños por produto;
* camiños energeticamente máis favorables.

### Xeración automática de ficheiros para Pilgrim

A partir das rutas seleccionadas, o script constrúe automaticamente o ficheiro:

```text
chem_C6H6/pif.chem
```

compatible con Pilgrim.

Ademais, copia ao directorio correspondente todas as estruturas de mínimos, estados de transición e produtos necesarias para a simulación cinética.

## Aplicación neste traballo

Esta ferramenta empregouse para reducir redes reaccionais inicialmente moi extensas a un conxunto de rutas mecanisticamente interpretables. A súa utilización permitiu identificar produtos relevantes, analizar a competencia entre diferentes canles reaccionais e construír modelos cinéticos adecuados para a simulación da pirólise unimolecular do benceno.
