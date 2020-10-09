# Modelo numérico de interior estelar

Modelo numérico que resuelve mediante el método de diferencias, el sistema de ecuaciones diferenciales acopladas que gobiernan el transporte de energía en el interior de una estrella.

<p align="center">
  <img src="https://github.com/DanielLapido/Modelo_estelar/blob/master/Figures/Ecuaciones.jpg" width=50% height=50%>
</p>

El modelo resuelve por separado las ecuaciones diferenciales para la región radiativa de la estrella y para la región convectiva y posteriormente ajusta correctamente las soluciones.
Se combina una integración numérica desde la superficie de la estrella hacia el centro
y otra integración desde el centro hasta la superficie. Posteriormente, ambas soluciones se unen
en la frontera entre la zona radiativa y la zona convectiva de la estrella.

Partiendo del valor estimado de la temperatura central, se busca un valor óptimo de la temperatura central
para el cual se minimizan las diferencias entre las dos soluciones.

El modelo resultante devuelve los valores de la temperatura, presión, masa, luminosidad y producción de energía
en función de la distancia al centro de la estrella.

Finalmente, se realiza una búsqueda de los valores óptimos de Radio total y Luminosidad total en una
malla entorno a los valores estimados inicialmente de Radio total y Luminosidad total.

La integración en la zona radiativa consta de tres fases:

* A.1.1 Envoltura radiativa: masa y luminosidad constantes
* A.1.2 Envoltura radiativa: masa variable y luminosidad constante
* A.1.3 Envoltura radiativa con masa y luminosidad variables

Algoritmo para las capas exteriores de la estrella. Masa y luminosidad más o menos constantes. 
<p align="center">
  <img src="https://github.com/DanielLapido/Modelo_estelar/blob/master/Figures/MasayLum_constantes.jpg" width=50% height=50%>
</p>

Los algoritmos para las fases A.1.2 y A.1.3 son similares a este.

Si se continúa la integración desde la superficie atravesando a continuación la zona convectiva, se obtiene una solución no válida para los parámetros físicos en el núcleo.
Por ello, se inicia una integración desde el núcleo hacia la superficie y se unen ambas soluciones.

El algoritmo para la zona convectiva es el siguiente:
<p align="center">
  <img src="https://github.com/DanielLapido/Modelo_estelar/blob/master/Figures/MasayLum_variables.jpg" width=50% height=50%>
</p>


Al unir las soluciones obtenidas mediante la integración hacia adentro de la estrella y hacia afuera, se observa una clara discontinuidad en la transición entre la zona convectiva y la zona radiativa. Por ello, se ajustan los parámetros estimados de Temperatura central, Luminosidad Total y Masa Total de la estrella.

Evolución de los parámetros sin realizar el ajuste:
<p align="center">
  <img src="https://github.com/DanielLapido/Modelo_estelar/blob/master/Figures/evolucion_sinajuste.jpeg" width=50% height=50%>
</p>

Evolución de los parámetros realizando el ajuste:
<p align="center">
  <img src="https://github.com/DanielLapido/Modelo_estelar/blob/master/Figures/evolucion_condobleajuste.jpeg" width=50% height=50%>
</p>
