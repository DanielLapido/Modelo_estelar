# Modelo numérico de interior estelar

Modelo numérico que resuelve mediante el método de diferencias, el sistema de ecuaciones diferenciales acopladas que gobiernan el interior de una estrella.

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
