# TFG-Hugo-Vallejo-Misa
# Scripts desenvolvidos

Durante este traballo desenvolvéronse varias ferramentas auxiliares para automatizar a análise de redes de reacción e o tratamento de procesos sen barreira. Estes scripts forman parte da metodoloxía empregada e permitiron reducir significativamente o traballo manual asociado ao estudo de redes reaccionais complexas.

## grafo.py

Ferramenta desenvolvida para a análise automática de redes de reacción xeradas mediante AutoMeKin/Pilgrim.

O script:

* Constrúe un grafo reaccional a partir dos ficheiros `RXNet.cg` e `RXNet.barrless`.
* Integra nunha única rede reaccións con e sen barreira.
* Identifica automaticamente camiños de reacción entre especies seleccionadas e produtos finais.
* Aplica filtros mecanísticos e enerxéticos para eliminar rutas pouco relevantes.
* Detecta pescozos de botella enerxéticos e camiños de menor barreira.
* Xera automaticamente ficheiros `pif.chem` compatibles con Pilgrim.
* Copia as estruturas necesarias para a posterior simulación cinética.

Esta ferramenta permitiu transformar redes inicialmente moi extensas en conxuntos de rutas mecanisticamente interpretables e aptas para cálculos cinéticos posteriores.

## VTSTS.py

Ferramenta desenvolvida para a análise de perfís VTST obtidos a partir de escaneos relaxados de Gaussian.

O script:

* Procesa automaticamente todos os puntos dun escaneo.
* Calcula funcións de partición e enerxías libres de Gibbs.
* Obtén perfís de enerxía libre dependentes da temperatura.
* Estima constantes cinéticas mediante a aproximación GFEA implementada en Pilgrim.
* Permite aplicar correccións automáticas a frecuencias vibracionais baixas cando resulta necesario.
* Xera táboas e gráficos de resultados en formato Excel.

Esta ferramenta empregouse para estudar reaccións inicialmente clasificadas como procesos sen barreira e avaliar a súa posible descrición mediante teoría variacional do estado de transición (VTST).

## Contribución metodolóxica

Os scripts desenvolvidos non constitúen unicamente ferramentas de postprocesado, senón que forman parte da metodoloxía empregada para a análise cinética da pirólise do benceno. En particular, `grafo.py` permitiu identificar automaticamente rutas mecanisticamente relevantes dentro de redes de reacción complexas, mentres que `VTSTS.py` se utilizou para caracterizar procesos sen barreira mediante cálculos VTST dependentes da temperatura. Ambas as ferramentas facilitaron a construción de modelos cinéticos consistentes e o tratamento sistemático de redes reaccionais de gran tamaño.

