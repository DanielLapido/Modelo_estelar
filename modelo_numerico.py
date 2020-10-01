# MODELO NUMÉRICO DE INTERIOR ESTELAR
# AUTOR: DANIEL LAPIDO MARTÍNEZ

# Modelo numérico que resuelve mediante el método de diferencias
# el sistema de ecuaciones diferenciales acopladas que gobiernan el interior de una estrella.

# El modelo resuelve por separado las ecuaciones diferenciales para la región radiativa de la estrella y
# para la región convectiva y posteriormente ajusta correctamente las soluciones.
# Se combina una integración numérica desde la superficie de la estrella hacia el centro
# y otra integración desde el centro hasta la superficie. Posteriormente, ambas soluciones se unen
# en la frontera entre la zona radiativa y la zona convectiva de la estrella.
# Partiendo del valor estimado de la temperatura central, se busca un valor óptimo de la temperatura central
# para el cual se minimizan las diferencias entre las dos soluciones.

# El modelo resultante devuelve los valores de la temperatura, presión, masa, luminosidad y producción de energía
# en función de la distancia al centro de la estrella.

# Finalmente, se realiza una búsqueda de los valores óptimos de Radio total y Luminosidad total en una
# malla entorno a los valores estimados inicialmente de Radio total y Luminosidad total.

from math import sqrt
from tabulate import tabulate
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
import copy

#VALORES INICIALES Y CONSTANTES
# Cambiando estas constantes se pueden resolver las ecuaciones para la estrella deseada.

Rtot = 12  # Radio total de la estrella
Ltot = 40  # Luminosidad total de la estrella
Tcc = 1.5  # Temperatura central de la estrella
Tc = copy.deepcopy(Tcc)  # Creamos copia de la temperatura central,
                         # pues necesitamos el valor original de Tcc más adelante

Mtot=5  #Masa total. Las unidades son 10^33 gramos
X=0.75        #Fracción de Hidrógeno
Y=0.20        #Fracción de Helio
Z=1-X-Y         #Fracción de Metales
Rini = 0.9*Rtot   #No tomamos radio total para evitar problemas de convergencia
h = -Rini/100   #Este es el paso de integración. Hemos tomado 100 capas.  Es negativo al integrar desde la #superficie hacia el centro

#Cálculo del peso molecular medio

#Supondremos que la estrella es homogénea en composicion química y que el material esta completamente ionizado.

mu=1/(2*X+3/4*Y+1/2*Z)
Cm =0.01523*mu
Cp = 8.084*mu
Ct = 0.01679*Z*(1+X)*mu**2
Ct_convec = 3.234*mu


#------------#

#Ritmo de generación de energía

#Debemos calcular la energía producida por el ciclo protó-protón y por el ciclo CN y quedarnos con el que sea mayor

def funcion_energia(T):
    """Esta función me calcula la energía producida por cada ciclo y
    me devuelve el valor correspondiente al ciclo que más energía genere.
    """

    #T viene dada en unidades de 10^7 K y para el cálculo se emplean unidades de 10^6 K:
    T=T*10

    #Calculamos primero energía del ciclo pp

    if T>=4 and T<6:
        epsilonpp=10**(-6.84)
        nupp=6
    elif T>=6 and T<9.5:
        epsilonpp = 10 ** (-6.04)
        nupp = 5
    elif T>=9.5 and T<12:
        epsilonpp = 10 ** (-5.56)
        nupp = 4.5
    elif T>=12 and T<16.5:
        epsilonpp = 10 ** (-5.02)
        nupp = 4
    elif T>=16.5 and T<=24:
        epsilonpp = 10 ** (-4.40)
        nupp = 3.5
    else:
        epsilonpp = 0
        nupp = 0

    energiapp=epsilonpp*X*X*T**(nupp)

    # Ahora calculamos la energía del ciclo CN

    if T>=12 and T<16:
        epsiloncn=10**(-22.2)
        nucn=20
    elif T>=16 and T<22.5:
        epsiloncn = 10 ** (-19.8)
        nucn = 18
    elif T>=22.5 and T<27.5:
        epsiloncn = 10 ** (-17.1)
        nucn = 16
    elif T>=27.5 and T<36:
        epsiloncn = 10 ** (-15.6)
        nucn = 15
    elif T>=36 and T<=50:
        epsiloncn = 10 ** (-12.5)
        nucn = 13
    else:
        epsiloncn = 0
        nucn = 0

    energiacn = epsiloncn * X * (1/3*Z) * (T) ** nucn

    if energiapp>energiacn:
        energia = energiapp
        epsilon = epsilonpp
        nu = nupp
        Xa = X
        Xb = X
        ciclo = "pp"

    elif energiacn>energiapp:
        energia = energiacn
        epsilon = epsiloncn
        nu = nucn
        Xa = X
        Xb = 1/3*Z
        ciclo= "CN"

    elif energiacn == 0 and energiapp == 0:
        energia = 0
        epsilon = 0
        nu = 0
        Xa = 0
        Xb = 0
        ciclo = "--"

    return(energia,epsilon,nu,Xa,Xb,ciclo)




# Creamos lista vacía con los valores del radio:
radio=[]

#Lista con los números de capa
I = list(range(0,101))

#PRIMERAS 3 CAPAS

#Calculamos temperatura y presión para las capas exteriores de la estrella

temperaturas=[]
presiones=[]
masa=[]
luminosidad=[]
derivadamasa = []
derivadaluminosidad =[]
derivadapresion = []
derivadatemperatura = []
fase=[]
E=[]
n =[] #lista del parametro n+1


def paso1():
    """"Paso 1. Esta función estima la presión y la temperatura en las tres primeras capas,
     aproximando la masa y la luminosidad por las totales
    """
    for i in range(0, 3):
        r = Rini + i*h
        A1 = 1.9022 * mu * Mtot
        A2 = 10.645 * sqrt(Mtot / (mu * Z * (1 + X) * Ltot))
        T = A1 * ((1 / r) - (1 / Rtot))  # Temperatura primeras capas
        P = A2 * T ** 4.25  # Presión primeras capas
        M = Mtot
        L = Ltot

        temperaturas.append(T)
        presiones.append(P)
        masa.append(M)
        luminosidad.append(L)
        radio.append(r)

        fm = 0  # En estas primeras capas la consideramos constante
        fl = 0  # En estas primeras capas la consideramos constante
        fp = -Cp * (P / T) * (Mtot / r ** 2)
        ft = -Ct * (P ** 2 / T ** 8.5) * (Ltot / r ** 2)

        derivadamasa.append(fm)
        derivadaluminosidad.append(fl)
        derivadapresion.append(fp)
        derivadatemperatura.append(ft)
        fase.append("INICIO")
        infoenergia = funcion_energia(T)
        E.append(infoenergia[5])
        n.append("NaN")  # No usamos aún el parámetro n+1



def paso2():
    """"Paso 2. Esta función estima la presión
       y la temperatura de una capa conociendo las tres capas anteriores
       """

    # Estimación de la presión en la capa i+1
    p_est = presiones[i] + h * derivadapresion[i] + 1 / 2 * (h * derivadapresion[i] - h * derivadapresion[i - 1]) + \
            5 / 12 * ((h * derivadapresion[i] - h * derivadapresion[i - 1]) - (
            h * derivadapresion[i - 1] - h * derivadapresion[i - 2]))

    # Estimación de la temperatura en la capa i+1
    t_est = temperaturas[i] + h * derivadatemperatura[i] + 1 / 2 * (
            h * derivadatemperatura[i] - h * derivadatemperatura[i - 1])

    return (p_est, t_est)

def paso4(p_est,t_est,m,r):
    """Paso 4. Esta función me calcula la presión y su derivada de la capa i+1 empleando la masa de la capa
    y las estimaciones de la presión y la temperatura.
    """

    fp = -Cp * (p_est / t_est) * (m / r ** 2)  # calcula la derivada de la presion en la capa i+1

    p_cal = presiones[i] + h * fp - 1 / 2 * (h * fp - h * derivadapresion[i])  # calcula la presion de la capa i+1

    return (fp,p_cal)



def paso7(p_cal,t_est,L,r):
    """Paso 7. Esta función me calcula la temperatura de la siguiente capa
    a partir de las derivadas f_i+1
    """
    ft = -Ct * (p_cal**2 / t_est**8.5) * (L / r ** 2)  # calcula la derivada de la temperatura en la capa i+1
    t_cal = temperaturas[i] + h * ft - 1 / 2 * (h * ft - h * derivadatemperatura[i])  # calcula la temperatura de la capa i+1

    return (ft, t_cal)



def paso3(p, t, m, derivadamasa, r):
    """Paso3. Esta función me calcula la masa conociendo p_cal, t_cal y la masa de la capa anterior.
    """
    fm = Cm * (p / t) * r**2  # calcula la derivada de la masa en la capa i+1
    m_cal = m + h * fm - 1 / 2 * (h * fm - h * derivadamasa[i])  # calcula la temperatura de la capa i+1
    return (fm, m_cal)


#EJECUCION DEL CODIGO

#INICIO PROGRAMA

paso1()
i = len(presiones)- 1  # Quiero que i tenga el valor de la tercera capa, es decir, i=2
m_cal = Mtot

print("FASE A.1.1")
#FASE RADIATIVA A.1.1.
#Consideramos masa y luminosidad constante. La masa empezará a variar mucho antes de que lo haga
#la luminosidad, así pues, emplearemos este algoritmo hasta que la masa deje de considerarse constante

while abs(Mtot-m_cal) / Mtot < 0.0001:

    r = Rini + (i+1)*h

    #Paso2
    estimacion = paso2()  # Vector de dos elementos que almacena la presión estimada y la temperatura estimada
    p_est = estimacion[0]
    t_est = estimacion[1]

    #Paso 4  En la fase A.1.1 empleamos como masa de la capa Mtot
    calculaP = paso4(p_est, t_est, Mtot, r)  # Vector que almacena la derivada de la presión y la presión calculada
    fp = calculaP[0]
    p_cal = calculaP[1]

    #Paso 5
    while abs(p_cal - p_est) / p_cal > 0.0001:      #Verifica si la presión calculada es suficientemente próxima a la estimada.
        p_est = p_cal                               #En caso de no serlo, vuelve al paso 4 y obtiene nuevos valores de la presión
        calculaP = paso4(p_est, t_est, Mtot, r)              #y su derivada hasta que se cumpla la condición
        fp = calculaP[0]
        p_cal = calculaP[1]

    #Paso 7
    calculaT = paso7(p_cal, t_est, Ltot, r)   # Vector que almacena la derivada de la temperatura y la temperatura calculada
    ft = calculaT[0]
    t_cal = calculaT[1]

    #Paso 8
    while abs(t_cal - t_est) / t_cal > 0.0001:
        t_est = t_cal
        calculaP = paso4(p_est, t_est, Mtot, r) #Volvemos a paso 4
        fp = calculaP[0]
        p_cal = calculaP[1]
        while abs(p_cal - p_est) / p_cal > 0.0001:
            p_est = p_cal
            calculaP = paso4(p_est, t_est,Mtot, r)
            fp = calculaP[0]
            p_cal = calculaP[1]
        calculaT = paso7(p_cal, t_est, Ltot, r) #Llegamos de nuevo a paso 7 y verifica si se cumple la condición
        ft = calculaT[0]
        t_cal = calculaT[1]

    #Paso 3. En la fase A.1.1 empleamos p_cal, t_cal y Mtot
    calculaM = paso3(p_cal, t_cal, Mtot, derivadamasa, r)
    fm = calculaM[0]
    m_cal = calculaM[1]

    #Si la masa no ha cambiado mucho, tomamos como válidos los valores de la capa y asumimos que la masa
    #de la capa es la masa total.
    if abs(Mtot-m_cal) / Mtot < 0.0001:
        presiones.append(p_cal)
        temperaturas.append(t_cal)
        masa.append(Mtot)                #La masa es la masa total
        luminosidad.append(Ltot)         #La luminosidad es la luminosidad total
        derivadapresion.append(fp)
        derivadatemperatura.append(ft)
        derivadamasa.append(fm)
        derivadaluminosidad.append(0)     #La luminosidad se considera constante
        fase.append("A.1.1.")
        infoenergia = funcion_energia(t_cal)
        E.append(infoenergia[5])
        n.append("NaN") #No usamos aún el parámetro n+1
        radio.append(r)



    i = i+1


print("--------------------------------------")


#FASE A.1.2.

#Ya no podemos considerar la masa constante. Emplearemos este algoritmo hasta que la
#luminosidad deje de considerarse constante


def paso6(p_cal , t, L, derivadaluminosidad, r):
    """Paso 6. Fase A.1.2. Esta función me calcula la luminosidad de la siguiente capa
    a partir de las derivadas f_i+1. Además emplea la función energia_funcion para calcular
    los parámetros:(energia,epsilon,nu,Xa,Xb)
    """
    infoenergia = funcion_energia(t)
    energia = infoenergia[0]
    epsilon = infoenergia[1]
    nu = infoenergia[2]
    Xa = infoenergia[3]
    Xb = infoenergia[4]
    Cl = 0.01845 * epsilon*Xa*Xb*10**(nu)*mu**2

    #Con todos estos parámetros podemos calcular la luminosidad y su derivada
    fL = Cl*p_cal**2*t**(nu-2)*r**2
    l_cal= L + h*fL - 1/2 * (h * fL - h * derivadaluminosidad[i]) -\
        1 / 12 * ((h * fL - h * derivadaluminosidad[i]) -
        (h * derivadaluminosidad[i] - h * derivadaluminosidad[i - 1]))

    return(fL,l_cal)


#EJECUCION DEL CODIGO
print("FASE A.1.2")

l_cal = Ltot
i = len(presiones) - 1

while abs(Ltot-l_cal) / Ltot < 0.0001:

    r = Rini + (i + 1) * h

    #Paso 2
    estimacion = paso2()  # Vector de dos elementos que almacena la presión estimada y la temperatura estimada
    p_est = estimacion[0]
    t_est = estimacion[1]

    #Paso3. En la fase A.1.2. empleamos p_est, t_est y la masa de la capa i
    calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
    fm = calculaM[0]
    m_cal = calculaM[1]

    #Paso4. En la fase A.1.2 empleamos la masa calculada en la capa i+1
    calculaP = paso4(p_est, t_est, m_cal, r)  # Vector que almacena la derivada de la presión y la presión calculada
    fp = calculaP[0]
    p_cal = calculaP[1]

    #Paso 5
    while abs(p_cal - p_est) / p_cal > 0.0001:                           #Verifica si la presión calculada es suficientemente próxima a la estimada.
        p_est = p_cal                                                    #En caso de no serlo, vuelve al paso 3 y obtiene nuevos valores de la masa,
        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)         #de la presión y sus derivadas hasta que se cumpla la condición
        fm = calculaM[0]
        m_cal = calculaM[1]
        calculaP = paso4(p_est, t_est, m_cal, r)
        fp = calculaP[0]
        p_cal = calculaP[1]

    #Paso 7. Utilizamos el valor de Ltot, pues consideramos que permanece constante
    calculaT = paso7(p_cal, t_est, Ltot, r)   # Vector que almacena la derivada de la temperatura y la temperatura calculada
    ft = calculaT[0]
    t_cal = calculaT[1]

    #Paso 8
    while abs(t_cal - t_est) / t_cal > 0.0001:
        t_est = t_cal
        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  #Volvemos al paso 3
        fm = calculaM[0]
        m_cal = calculaM[1]
        calculaP = paso4(p_est, t_est, m_cal, r) #Pasamos al 4
        fp = calculaP[0]
        p_cal = calculaP[1]
        while abs(p_cal - p_est) / p_cal > 0.0001:      #Verifica si p_cal próximo a p_est
            p_est = p_cal
            calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)    #Paso 3
            fm = calculaM[0]
            m_cal = calculaM[1]
            calculaP = paso4(p_est, t_est, m_cal, r)
            fp = calculaP[0]
            p_cal = calculaP[1]
        calculaT = paso7(p_cal, t_est, Ltot, r) #Llegamos de nuevo a paso 7 y verifica si se cumple la condición
        ft = calculaT[0]
        t_cal = calculaT[1]

    #Paso 6 Calculamos la luminosidad y su derivada en la capa i+1.
    #En la fase A.1.2. empleamos t_cal.
    # Como la luminosidad se considera constante, la luminosidad de la capa anterior es Ltot
    calculaL = paso6(p_cal, t_cal, Ltot, derivadaluminosidad, r)
    fl = calculaL[0]
    l_cal = calculaL[1]


    #Si la luminosidad no ha cambiado mucho, tomamos como válidos los valores de la capa y asumimos que la luminosidad
    #de la capa es la luminosidad total
    if abs(Ltot-l_cal) / Ltot < 0.0001:
        presiones.append(p_cal)
        temperaturas.append(t_cal)
        masa.append(m_cal)
        luminosidad.append(Ltot)         #La luminosidad es la luminosidad total
        derivadapresion.append(fp)
        derivadatemperatura.append(ft)
        derivadamasa.append(fm)
        derivadaluminosidad.append(fl)
        fase.append("A.1.2.")
        infoenergia = funcion_energia(t_cal)
        E.append(infoenergia[5])
        n.append("NaN")  # No usamos aún el parámetro n+1
        radio.append(r)

    i = i + 1

print("--------------------------------------")

#EJECUCIÓN DEL CÓDIGO
print("FASE A.1.3")


#Ahora tanto la luminosidad como la masa varían en cada capa.
i = len(presiones)-1
a = 3  #Esta variable es n+1. Hacemos que tome un valor mayor a 2.5 para que comience el bucle
while a > 2.5:

    r = Rini + (i + 1) * h

    # Paso 2
    estimacion = paso2()  # Vector de dos elementos que almacena la presión estimada y la temperatura estimada
    p_est = estimacion[0]
    t_est = estimacion[1]
    # Paso3. En la fase A.1.3. empleamos p_est, t_est y la masa de la capa i
    calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
    fm = calculaM[0]
    m_cal = calculaM[1]

    # Paso4. En la fase A.1.3 empleamos la masa calculada en la capa i+1
    calculaP = paso4(p_est, t_est, m_cal, r)  # Vector que almacena la derivada de la presión y la presión calculada
    fp = calculaP[0]
    p_cal = calculaP[1]

    # Paso 5
    while abs(p_cal - p_est) / p_cal > 0.0001:                    # Verifica si la presión calculada es suficientemente próxima a la estimada.
        p_est = p_cal                                             # En caso de no serlo, vuelve al paso 3 y obtiene nuevos valores de la masa,
        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # de la presión y sus derivadas hasta que se cumpla la condición
        fm = calculaM[0]
        m_cal = calculaM[1]
        calculaP = paso4(p_est, t_est, m_cal, r)
        fp = calculaP[0]
        p_cal = calculaP[1]

    #Paso 6. Una vez que p_cal = p_est, pasamos a calcular la luminosidad.
    # Esta vez empleamos T_est, pues t_cal se hallará en el paso 7.
    #La luminosidad varía y ya no podemos utilizar Ltot.
    calculaL = paso6(p_cal, t_est, luminosidad[i], derivadaluminosidad, r)
    fl = calculaL[0]
    l_cal = calculaL[1]

    #Paso 7. Dado que ahora la luminosidad no es constante, debemos emplear la luminosidad calculada en el
    #paso anterior para calcular la temperatura.
    calculaT = paso7(p_cal, t_est, l_cal, r)  # Vector que almacena la derivada de la temperatura y la temperatura calculada
    ft = calculaT[0]
    t_cal = calculaT[1]

    # Paso 8
    while abs(t_cal - t_est) / t_cal > 0.0001:
        t_est = t_cal
        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # Volvemos al paso 3
        fm = calculaM[0]
        m_cal = calculaM[1]
        calculaP = paso4(p_est, t_est, m_cal, r)  # Pasamos al 4
        fp = calculaP[0]
        p_cal = calculaP[1]
        while abs(p_cal - p_est) / p_cal > 0.0001:  # Verifica si p_cal próximo a p_est
            p_est = p_cal
            calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # Paso 3
            fm = calculaM[0]
            m_cal = calculaM[1]
            calculaP = paso4(p_est, t_est, m_cal, r)
            fp = calculaP[0]
            p_cal = calculaP[1]
        calculaL = paso6(p_cal, t_est,luminosidad[i], derivadaluminosidad, r)
        fl = calculaL[0]
        l_cal = calculaL[1]
        calculaT = paso7(p_cal, t_est, l_cal, r)  # Llegamos de nuevo a paso 7 y verifica si se cumple la condición
        ft = calculaT[0]
        t_cal = calculaT[1]

    #Finalmente, debemos comprobar si o la hipótesis inicial de transporte radiativo deja de ser válida,
    #vamos a calcular el párametro n+1 definido en el Novotny

    #Llamamos a=n+1

    a = t_cal*h*fp/(p_cal*h*ft)

    if a > 2.5:
        presiones.append(p_cal)
        temperaturas.append(t_cal)
        masa.append(m_cal)
        luminosidad.append(l_cal)
        derivadapresion.append(fp)
        derivadatemperatura.append(ft)
        derivadamasa.append(fm)
        derivadaluminosidad.append(fl)
        fase.append("A.1.3.")
        infoenergia = funcion_energia(t_cal)
        E.append(infoenergia[5])
        n.append(a)
        radio.append(r)

    else:    #Cálculos conducen a n+1 ≤ 2.5 y datos de la capa no son válidos
        t_convec = t_cal
        p_convec = p_cal
        n.append(a)   #añadimos el valor del parámetro n+1 de la primera capa convectiva

    r_transicion = Rini + i * h  # Guardamos el valor del radio para el que se produce la transición entre fase radiativa y convectiva
    transicion = [r_transicion, i]

    i = i + 1

print('El radio en el que se produce la transición es:', transicion[0], 'y ocurre en la capa:', transicion[1])

print("--------------------------------------")

#FASE A.2. NUCLEO CONVECTIVO

print("FASE A.2")
#Cuando la aplicación de la fase A.1.3 conduce a un cálculo del parámetro n+1 ≤ 2.5 sabemos que los cálculos
#realizados en dicha capa no son válidos y debemos repetirlos empleando el algoritmo A.2

#No obstante, los valores de presión y temperatura calculados en esa capa sí que nos sirven para
#estimar la constante K que nos relaciona estos dos parametros fíısicos en un polítropo.
#Asumiendo un índice adiabático γ = 5/3 (correspondiente a un gas perfecto monoatomico):

K = p_convec/t_convec**2.5   #Este valor calculado en la primera capa convectiva se asumirá constante

def paso2bis(temperaturas,derivadatemperatura):
    """Paso2bis. Esta función me estima el valor de la temperatura en las sucesivas capas
    de la fase convectiva. A diferencia del paso2, la presión no se estima en este paso,
    pues se calcula directamente a partir de la expresión del polítropo
    """

    # Estimación de la temperatura en la capa i+1 del núcleo convectivo
    t_est = temperaturas[i] + h * derivadatemperatura[i] + 1 / 2 * (
            h * derivadatemperatura[i] - h * derivadatemperatura[i - 1])

    return t_est

def paso7bis(M, temperaturas, derivadatemperatura, r):
    """Paso 7bis. Esta función me calcula la temperatura y su derivada de las sucesivas capas
    de la fase convectiva. La ecuación diferencial para la temperatura en el caso convectivo
    es diferente de la del radiativo.
    """
    if r > 0:
        ft = -Ct_convec * M / r ** 2  # calcula la derivada de la temperatura en la capa convectiva i+1
        t_cal = temperaturas[i] + h * ft - 1 / 2 * (
                h * ft - h * derivadatemperatura[i])  # calcula la temperatura de la capa convectiva i+1

    else: #Cuando r = 0, no podemos usar la ecuación anterior y asumimos t_cal = t_est
        t_cal = t_est
        ft = 0

    return (ft, t_cal)

#EJECUCIÓN DEL ALGORITMO A.2


i = len(presiones) - 1


while r >= 0:

    r = Rini + (i + 1) * h

    #Paso 2
    t_est = paso2bis(temperaturas, derivadatemperatura)
    # Estimamos la presión con expresión del polítropo
    p_est = K*t_est**2.5
    # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
    # expresión del polítropo.
    calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
    fm = calculaM[0]
    m_cal = calculaM[1]

    # No hay pasos 4 y 5 porque la presión ahora se calcula directamente con la expresión del polítropo.
    # Dado que el gradiente de temperatura en el caso convectivo no depende de la luminosidad, el paso 6 se
    # incorporará después del paso 8.

    # Paso 7bis
    calculaT = paso7bis(m_cal, temperaturas, derivadatemperatura, r)
    ft = calculaT[0]
    t_cal = calculaT[1]

    # Paso 8
    while abs(t_cal - t_est) / t_cal > 0.0001:
        t_est = t_cal
        p_est = K * t_est ** 2.5  #Volvemos a estimar la presión
        #Repetimos paso 3
        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
        fm = calculaM[0]
        m_cal = calculaM[1]
        #Repetimos paso 7bis
        calculaT = paso7bis(m_cal, temperaturas, derivadatemperatura, r)
        ft = calculaT[0]
        t_cal = calculaT[1]    #En este punto verifica si t_cal es prácticamente igual a t_est. De ser así sale del bucle

    #Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
    if r > 0:
        p_cal = K * t_cal ** 2.5

    else: #En el caso del centro de la estrella basta con utilizar p_cal = p_est y t_cal = t_est
        t_cal = paso2bis(temperaturas, derivadatemperatura)
        p_cal = K * t_cal ** 2.5
        # Paso 3
        calculaM = paso3(p_cal, t_cal, masa[i], derivadamasa, r)
        fm = calculaM[0]
        m_cal = calculaM[1]

    #Paso 6. La ecuación para el cálculo de la luminosidad es igual que en el caso radiativo, pero considerando la
    #expresión del polítropo
    calculaL = paso6(p_cal, t_cal, luminosidad[i], derivadaluminosidad, r)
    fl = calculaL[0]
    l_cal = calculaL[1]

    #Si los cálculos de la capa son válidos, los añadimos a nuestras listas
    if r >= 0:
        presiones.append(p_cal)
        temperaturas.append(t_cal)
        masa.append(m_cal)
        luminosidad.append(l_cal)
        derivadatemperatura.append(ft)
        derivadamasa.append(fm)
        derivadaluminosidad.append(fl)
        fase.append("CONVEC")
        infoenergia = funcion_energia(t_cal)
        E.append(infoenergia[5])
        n.append("NaN")  # No usamos el parámetro n+1
        radio.append(r)

    i = i + 1

datos={'E': E,'fase': fase,'i': I,'r':radio[0:101],'Presión': presiones, 'Temperatura': temperaturas,'Luminosidad':luminosidad,'Masa':masa,'n+1':n[0:101]}
print(tabulate(datos, headers='keys'))

#Al acercarnos al núcleo en las últimas capas, la presión y temperatura disminuyen. Además la luminosidad no se anula
#y la masa se hace negativa!

#La forma de arreglar este problema consiste en utilizar una temperatura central diferente a la estimada al
#final de la integracion desde arriba, con la esperanza de que en el límite entre la capa convectiva y la radiativa
#se obtengan unos parametros físicos similares a los calculados al integrar desde la superficie.


#VALORES EN LA FRONTERA (integración desde arriba hacia el centro de la estrella)
#Resulta útil calcular cuáles han sido los valores de los parametros físicos en la transición entre las regiones radiativa y convectiva de la estrella.
#Queremos interpolar linealmente los resultados obtenidos para obtener una predicción a un radio tal que el parámetro
#n+1 sea exactamente igual a 2.5.

print("n+1 vale", n[transicion[1]], "en la capa", transicion[1],"(radio = ", transicion[0], ")")
print("y n+1 vale", n[transicion[1] + 1], "en la capa", transicion[1] + 1, "(radio = ", radio[transicion[1] + 1], ")")
print("Realizamos interpolación lineal para hallar el radio en que n+1 = 2.5")

r_rad_conv = transicion[0] + (radio[transicion[1] + 1] - transicion[0])/(n[transicion[1] + 1] - n[transicion[1]])*\
           (2.5 - n[transicion[1]])

print("El radio interpolado vale:", r_rad_conv)

#Interpolamos linealmente el resto de parámetros físicos también:
p_rad_conv = presiones[transicion[1]] + \
             (presiones[transicion[1] + 1] - presiones[transicion[1]])/(n[transicion[1] + 1] - n[transicion[1]]) * \
             (2.5 - n[transicion[1]])
print("La presión interpolada vale:", p_rad_conv)

t_rad_conv = temperaturas[transicion[1]] + \
             (temperaturas[transicion[1] + 1] - temperaturas[transicion[1]])/(n[transicion[1] + 1] - n[transicion[1]])\
             * (2.5 - n[transicion[1]])
print("La temperatura interpolada vale:", t_rad_conv)

l_rad_conv = luminosidad[transicion[1]] + \
             (luminosidad[transicion[1] + 1] - luminosidad[transicion[1]])/(n[transicion[1] + 1] - n[transicion[1]])\
             * (2.5 - n[transicion[1]])
print("La luminosidad interpolada vale:", l_rad_conv)

m_rad_conv = masa[transicion[1]] + \
             (masa[transicion[1] + 1] - masa[transicion[1]])/(n[transicion[1] + 1] - n[transicion[1]]) * \
             (2.5 - n[transicion[1]])
print("La masa interpolada vale:", m_rad_conv)

print("--------------------------------------")

print("INTEGRACIÓN DESDE EL CENTRO")

h = Rini/100   #Este es el paso de integración. Hemos tomado 100 capas.  Es positivo al integrar desde el centro hacia la superficie

#PARA EL MODELO DE PRUEBA n+1 = 2.50093 en la capa i=81 (radio r=1.9665)y 2.379415
#en la capa i=82 (radio r=1.8630). Una simple interpolacion lineal predice que n+1=2.5 para un radio
#r_rad→r_conv =1.96571, y en el que podemos calcular los parámetros físicos por interpolación.

#PRIMERAS 3 CAPAS

#Creamos una función para el cálculo de los parámetros físicos en las capas centrales
#de una estrella, asumiendo transporte convectivo de energía

def parametros_nucleo(r,Tc):
    """Esta función me calcula los parámetros físicos en las capas centrales
    de una estrella, asumiendo transporte convectivo de energía
    """
    m = 0.005077*mu*K*Tc**1.5*r**3         #Cálculo de la masa
    t = Tc - 0.008207*mu**2*K*Tc**1.5*r**2   #Cálculo de la temperatura
    infoenergia = funcion_energia(t)
    energia = infoenergia[0]
    epsilon = infoenergia[1]
    nu = infoenergia[2]
    Xa = infoenergia[3]
    Xb = infoenergia[4]
    l = 0.006150*epsilon*Xa*Xb*10**nu*mu**2*K**2*Tc**(3+nu)*r**3     #Cálculo de la luminosidad
    p = K*t**2.5   #Cálculo de la presión
    dm = Cm*K*t**1.5*r**2  #Cálculo de derivada de masa
    Cl = 0.01845 * epsilon*Xa*Xb*10**(nu)*mu**2
    dl = Cl*K**2*t**(3+nu)*r**2 #Cálculo de derivada de luminosidad

    #Para las derivadas de la presión y temperatura hay que tener cuidado, pues ecuaciones divergen en 0
    if r == 0:
        dp = 0
        dt = 0
    else:
        dp = -Cp * K * t ** 1.5 * m / r ** 2
        dt = -Ct_convec* m / r ** 2

    return(m, l, t, p, dm, dl, dt, dp)

print("PRIMERAS 3 CAPAS")

masa_nucleo = []
luminosidad_nucleo = []
temperatura_nucleo = []
presion_nucleo = []
fase_nucleo = []
E_nucleo = []
I_nucleo = []
radio_nucleo = []
dm_nucleo =[]
dl_nucleo =[]
dt_nucleo =[]
dp_nucleo =[]

for i in range(0,3):
    r = 0+i*h
    parametros = parametros_nucleo(r, Tc)
    m = parametros[0]
    l = parametros[1]
    t = parametros[2]
    p = parametros[3]
    dm = parametros[4]
    dl = parametros[5]
    dt = parametros[6]
    dp = parametros[7]
    masa_nucleo.append(m)
    luminosidad_nucleo.append(l)
    temperatura_nucleo.append(t)
    presion_nucleo.append(p)
    info_energia = funcion_energia(t)
    E_nucleo.append(info_energia[5])
    fase_nucleo.append("CENTRO")
    radio_nucleo.append(r)
    dm_nucleo.append(dm)
    dl_nucleo.append(dl)
    dt_nucleo.append(dt)
    dp_nucleo.append(dp)
    I_nucleo.append(i)

datos_nucleo = {'E':E_nucleo,'fase':fase_nucleo,'i':I_nucleo,'r':radio_nucleo,'Presión': presion_nucleo, 'Temperatura': temperatura_nucleo,'Luminosidad':luminosidad_nucleo,'Masa':masa_nucleo}
print(tabulate(datos_nucleo, headers='keys'))


#CAPAS POSTERIORES

i = len(presion_nucleo)-1

while r <= transicion[0]:

    r = 0 + (i + 1) * h

    # Paso 1. Ya conocemos las presiones, temperaturas, masas, luminosidades y sus derivadas de las tres primeras capas

    # Paso 2bis. Estimamos la temperatura de la capa siguiente:
    t_est = paso2bis(temperatura_nucleo, dt_nucleo)

    # Estimamos la presión con expresión del polítropo
    p_est = K * t_est ** 2.5

    # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
    # expresión del polítropo.
    calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
    fm = calculaM[0]
    m_cal = calculaM[1]

    # Paso 7bis. A continuación podemos calcular la temperatura:
    calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
    ft = calculaT[0]
    t_cal = calculaT[1]

    # Paso 8. Comparamos la temperatura calculada con la estimada:
    while abs(t_cal - t_est) / t_cal > 0.0001:
        t_est = t_cal
        p_est = K * t_est ** 2.5  # Volvemos a estimar la presión
        # Repetimos paso 3
        calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
        fm = calculaM[0]
        m_cal = calculaM[1]
        # Repetimos paso 7bis
        calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
        ft = calculaT[0]
        t_cal = calculaT[1]  # En este punto verifica si t_cal es prácticamente igual a t_est. De ser así, sale del bucle

    # Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
    p_cal = K * t_cal ** 2.5

    # Paso 6. Conocida p_cal podemos hallar la luminosidad:
    calculaL = paso6(p_cal, t_cal, luminosidad_nucleo[i], dl_nucleo, r)
    fl = calculaL[0]
    l_cal = calculaL[1]

    if r <= transicion[0]:                 #Si los cálculos de la capa son válidos, los añadimos
        presion_nucleo.append(p_cal)
        temperatura_nucleo.append(t_cal)
        masa_nucleo.append(m_cal)
        luminosidad_nucleo.append(l_cal)
        dt_nucleo.append(ft)
        dm_nucleo.append(fm)
        dl_nucleo.append(fl)
        fase_nucleo.append("CONVEC")
        infoenergia = funcion_energia(t_cal)
        E_nucleo.append(infoenergia[5])
        radio_nucleo.append(r)
        I_nucleo.append(i + 1)
        transicion2 = [r, i+1]  # Guardamos el radio en el que comienza la fase radiativa

    i = i + 1



datos_nucleo = {'E':E_nucleo,'fase':fase_nucleo,'i':I_nucleo,'r':radio_nucleo,'Presión': presion_nucleo, 'Temperatura': temperatura_nucleo,'Luminosidad':luminosidad_nucleo,'Masa':masa_nucleo}
print(tabulate(datos_nucleo, headers='keys'))

#VALORES EN LA FRONTERA (integracion desde el centro hacia el exterior de la estrella)
#Al igual que hicimos al integrar desde la superficie, vamos a estimar los parametros físicos en la frontera
#entre la capa convectiva y la radiativa. Para ello volvemos a interpolar linealmente las soluciones que hemos
#encontrado.


r_conv_rad = radio_nucleo[transicion2[1]-1] + (transicion2[0] - radio_nucleo[transicion2[1]-1]) / \
             (n[transicion[1]] - n[transicion[1] + 1])*(2.5 - n[transicion[1] + 1])

print("El radio interpolado vale:", r_conv_rad)

#Interpolamos linealmente el resto de parámetros físicos también:
p_conv_rad = presion_nucleo[transicion2[1]-1] + \
             (presion_nucleo[transicion2[1]] - presion_nucleo[transicion2[1]-1])/ \
             (n[transicion[1]] - n[transicion[1] + 1])*(2.5 - n[transicion[1] + 1])

print("La presión interpolada vale:", p_conv_rad)

t_conv_rad = temperatura_nucleo[transicion2[1]-1] + \
             (temperatura_nucleo[transicion2[1]] - temperatura_nucleo[transicion2[1]-1])/ \
             (n[transicion[1]] - n[transicion[1] + 1])*(2.5 - n[transicion[1] + 1])
print("La temperatura interpolada vale:", t_conv_rad)

l_conv_rad = luminosidad_nucleo[transicion2[1]-1] + \
             (luminosidad_nucleo[transicion2[1]] - luminosidad_nucleo[transicion2[1]-1])/ \
             (n[transicion[1]] - n[transicion[1] + 1])*(2.5 - n[transicion[1] + 1])
print("La luminosidad interpolada vale:", l_conv_rad)

m_conv_rad = masa_nucleo[transicion2[1]-1] + \
             (masa_nucleo[transicion2[1]] - masa_nucleo[transicion2[1]-1])/ \
             (n[transicion[1]] - n[transicion[1] + 1])*(2.5 - n[transicion[1] + 1])
print("La masa interpolada vale:", m_conv_rad)


#Podemos ahora comparar las soluciones que hemos obtenido para r=rconv→rad=rrad→conv=1.96571, tanto
#al integrar desde arriba como desde abajo:

error_rel_p = (p_rad_conv - p_conv_rad)/p_rad_conv*100
error_rel_t = (t_rad_conv - t_conv_rad)/t_rad_conv*100
error_rel_l = (l_rad_conv - l_conv_rad)/l_rad_conv*100
error_rel_m = (m_rad_conv - m_conv_rad)/m_rad_conv*100
error_rel_tot = sqrt(error_rel_p**2 + error_rel_t**2 + error_rel_l**2 + error_rel_m**2)

errores = {'Error.relat P(%)': error_rel_p,'Error.relat T(%)': error_rel_t ,'Error.relat L(%)': error_rel_l,
           'Error.relat M(%)': error_rel_m}

for parametro, error in errores.items():
    print('{0:10} ==> {1:10f}'.format(parametro, error))
print("El error relativo total es:", error_rel_tot,"%")

print("El error relativo es demasiado grande.")
print("A continuación se calcula el valor de la temperatura central que minimiza este error:")

print("--------------------------------------------------------------------------------------")

#Vamos a hacer el mismo cálculo anterior para diferentes valores de Tc y veremos para cuál de ellos el error es menor:

Tcini = 1
Tcfin = 3

Tc = np.linspace(Tcini, Tcfin, 50) #Hacemos que tome 50 valores entre 1 y 3
Errores = []  # En esta lista almacenaremos los errores relativos totales para cada Tc
Lista_de_diccionarios = [] #Esta lista será una lista de diccionarios

iteraciones = 1  #Este contador mide el número de iteraciones que realizamos.
                 #En cada iteración obtenemos un valor más preciso de Tc


while iteraciones <= 3:    #Vamos a pedir que haga tres iteraciones. Podemos cambiar libremente este número

    intervalo = (Tcfin - Tcini)/50  # Este es el paso en Tc. (Nos será útil para afinar aún más el valor de Tc)

    for temperatura in Tc:

        masa_nucleo = []
        luminosidad_nucleo = []
        temperatura_nucleo = []
        presion_nucleo = []
        fase_nucleo = []
        E_nucleo = []
        I_nucleo = []
        radio_nucleo = []
        dm_nucleo = []
        dl_nucleo = []
        dt_nucleo = []
        dp_nucleo = []

        # PRIMERAS 3 CAPAS
        for i in range(0, 3):
            r = 0 + i * h
            parametros = parametros_nucleo(r, temperatura)
            m = parametros[0]
            l = parametros[1]
            t = parametros[2]
            p = parametros[3]
            dm = parametros[4]
            dl = parametros[5]
            dt = parametros[6]
            dp = parametros[7]
            masa_nucleo.append(m)
            luminosidad_nucleo.append(l)
            temperatura_nucleo.append(t)
            presion_nucleo.append(p)
            info_energia = funcion_energia(t)
            E_nucleo.append(info_energia[5])
            fase_nucleo.append("CENTRO")
            radio_nucleo.append(r)
            dm_nucleo.append(dm)
            dl_nucleo.append(dl)
            dt_nucleo.append(dt)
            dp_nucleo.append(dp)
            I_nucleo.append(i)

        # CAPAS POSTERIORES

        i = len(presion_nucleo) - 1

        while r <= transicion[0]:

            r = 0 + (i + 1) * h

            # Paso 1. Ya conocemos las presiones, temperaturas, masas, luminosidades y sus derivadas de las tres primeras capas

            # Paso 2bis. Estimamos la temperatura de la capa siguiente:
            t_est = paso2bis(temperatura_nucleo, dt_nucleo)

            # Estimamos la presión con expresión del polítropo
            p_est = K * t_est ** 2.5

            # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
            # expresión del polítropo.
            calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
            fm = calculaM[0]
            m_cal = calculaM[1]

            # Paso 7bis. A continuación podemos calcular la temperatura:
            calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
            ft = calculaT[0]
            t_cal = calculaT[1]

            # Paso 8. Comparamos la temperatura calculada con la estimada:
            while abs(t_cal - t_est) / t_cal > 0.0001:
                t_est = t_cal
                p_est = K * t_est ** 2.5  # Volvemos a estimar la presión
                # Repetimos paso 3
                calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
                fm = calculaM[0]
                m_cal = calculaM[1]
                # Repetimos paso 7bis
                calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
                ft = calculaT[0]
                t_cal = calculaT[
                    1]  # En este punto verifica si t_cal es prácticamente igual a t_est. De ser así, sale del bucle

            # Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
            p_cal = K * t_cal ** 2.5

            # Paso 6. Conocida p_cal podemos hallar la luminosidad:
            calculaL = paso6(p_cal, t_cal, luminosidad_nucleo[i], dl_nucleo, r)
            fl = calculaL[0]
            l_cal = calculaL[1]

            if r <= transicion[0]:  # Si los cálculos de la capa son válidos, los añadimos
                presion_nucleo.append(p_cal)
                temperatura_nucleo.append(t_cal)
                masa_nucleo.append(m_cal)
                luminosidad_nucleo.append(l_cal)
                dt_nucleo.append(ft)
                dm_nucleo.append(fm)
                dl_nucleo.append(fl)
                fase_nucleo.append("CONVEC")
                infoenergia = funcion_energia(t_cal)
                E_nucleo.append(infoenergia[5])
                radio_nucleo.append(r)
                I_nucleo.append(i + 1)
                transicion2 = [r, i + 1]  # Guardamos el radio en el que comienza la fase radiativa

            i = i + 1

        r_conv_rad = radio_nucleo[transicion2[1] - 1] + (transicion2[0] - radio_nucleo[transicion2[1] - 1]) / \
                     (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

        # Interpolamos linealmente el resto de parámetros físicos también:
        p_conv_rad = presion_nucleo[transicion2[1] - 1] + \
                     (presion_nucleo[transicion2[1]] - presion_nucleo[transicion2[1] - 1]) / \
                     (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

        t_conv_rad = temperatura_nucleo[transicion2[1] - 1] + \
                     (temperatura_nucleo[transicion2[1]] - temperatura_nucleo[transicion2[1] - 1]) / \
                     (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

        l_conv_rad = luminosidad_nucleo[transicion2[1] - 1] + \
                     (luminosidad_nucleo[transicion2[1]] - luminosidad_nucleo[transicion2[1] - 1]) / \
                     (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

        m_conv_rad = masa_nucleo[transicion2[1] - 1] + \
                     (masa_nucleo[transicion2[1]] - masa_nucleo[transicion2[1] - 1]) / \
                     (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

        # Calculamos los errores para esta temperatura central
        error_rel_p = (p_rad_conv - p_conv_rad) / p_rad_conv * 100
        error_rel_t = (t_rad_conv - t_conv_rad) / t_rad_conv * 100
        error_rel_l = (l_rad_conv - l_conv_rad) / l_rad_conv * 100
        error_rel_m = (m_rad_conv - m_conv_rad) / m_rad_conv * 100
        error_rel_tot = sqrt(error_rel_p ** 2 + error_rel_t ** 2 + error_rel_l ** 2 + error_rel_m ** 2)

        Errores.append(error_rel_tot)


    data = {"Tc": Tc, "Err_rel": Errores}
    df = pd.DataFrame(data)
    optimo = df[df.Err_rel == df.Err_rel.min()] #data frame con la Tc que hace mínimo el error
    Toptima = optimo.iloc[0, 0]  # Extraemos el valor óptimo de la temperatura central

    #Ahora ajustamos el intervalo para obtener un mejor valor de Tc:

    Tcini = Toptima - intervalo
    Tcfin = Toptima + intervalo

    Tc = np.linspace(Tcini, Tcfin, 50)  # Hacemos que tome 50 valores entre Tcini y Tcfin
    Errores = [] #Limpiamos la lista de errores para la siguiente iteracion

    iteraciones = iteraciones + 1

print("El valor óptimo de temperatura central es:", Toptima, ". El error relativo total es:", optimo.iloc[0, 1], "%")


#CÁLCULO DE LAS CAPAS CON LA TEMPERATURA CENTRAL ÓPTIMA.
#Volvemos a ejecutar el mismo código pero esta vez solo para la Toptima.
#Ahora sí guardamos todos los valores de cada capa

Tc = Toptima

#Limpiamos todas las listas y ejecutamos de nuevo el código
masa_nucleo = []
luminosidad_nucleo = []
temperatura_nucleo = []
presion_nucleo = []
fase_nucleo = []
E_nucleo = []
I_nucleo = []
radio_nucleo = []
dm_nucleo =[]
dl_nucleo =[]
dt_nucleo =[]
dp_nucleo =[]

for i in range(0,3):
    r = 0+i*h
    parametros = parametros_nucleo(r, Tc)
    m = parametros[0]
    l = parametros[1]
    t = parametros[2]
    p = parametros[3]
    dm = parametros[4]
    dl = parametros[5]
    dt = parametros[6]
    dp = parametros[7]
    masa_nucleo.append(m)
    luminosidad_nucleo.append(l)
    temperatura_nucleo.append(t)
    presion_nucleo.append(p)
    info_energia = funcion_energia(t)
    E_nucleo.append(info_energia[5])
    fase_nucleo.append("CENTRO")
    radio_nucleo.append(r)
    dm_nucleo.append(dm)
    dl_nucleo.append(dl)
    dt_nucleo.append(dt)
    dp_nucleo.append(dp)
    I_nucleo.append(i)


#CAPAS POSTERIORES

i = len(presion_nucleo)-1

while r <= transicion[0]:

    r = 0 + (i + 1) * h

    # Paso 1. Ya conocemos las presiones, temperaturas, masas, luminosidades y sus derivadas de las tres primeras capas

    # Paso 2bis. Estimamos la temperatura de la capa siguiente:
    t_est = paso2bis(temperatura_nucleo, dt_nucleo)

    # Estimamos la presión con expresión del polítropo
    p_est = K * t_est ** 2.5

    # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
    # expresión del polítropo.
    calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
    fm = calculaM[0]
    m_cal = calculaM[1]

    # Paso 7bis. A continuación podemos calcular la temperatura:
    calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
    ft = calculaT[0]
    t_cal = calculaT[1]

    # Paso 8. Comparamos la temperatura calculada con la estimada:
    while abs(t_cal - t_est) / t_cal > 0.0001:
        t_est = t_cal
        p_est = K * t_est ** 2.5  # Volvemos a estimar la presión
        # Repetimos paso 3
        calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
        fm = calculaM[0]
        m_cal = calculaM[1]
        # Repetimos paso 7bis
        calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
        ft = calculaT[0]
        t_cal = calculaT[1]  # En este punto verifica si t_cal es prácticamente igual a t_est. De ser así, sale del bucle

    # Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
    p_cal = K * t_cal ** 2.5

    # Paso 6. Conocida p_cal podemos hallar la luminosidad:
    calculaL = paso6(p_cal, t_cal, luminosidad_nucleo[i], dl_nucleo, r)
    fl = calculaL[0]
    l_cal = calculaL[1]

    if r <= transicion[0]:                 #Si los cálculos de la capa son válidos, los añadimos
        presion_nucleo.append(p_cal)
        temperatura_nucleo.append(t_cal)
        masa_nucleo.append(m_cal)
        luminosidad_nucleo.append(l_cal)
        dt_nucleo.append(ft)
        dm_nucleo.append(fm)
        dl_nucleo.append(fl)
        fase_nucleo.append("CONVEC")
        infoenergia = funcion_energia(t_cal)
        E_nucleo.append(infoenergia[5])
        radio_nucleo.append(r)
        I_nucleo.append(i + 1)
        transicion2 = [r, i+1]  # Guardamos el radio en el que comienza la fase radiativa

    i = i + 1



datos_nucleo = {'E':E_nucleo,'fase':fase_nucleo,'i':I_nucleo,'r':radio_nucleo,'Presión': presion_nucleo, 'Temperatura': temperatura_nucleo,'Luminosidad':luminosidad_nucleo,'Masa':masa_nucleo}

#UNIÓN DE SOLUCIONES

#Ahora vamos a unir las soluciones de la integración desde la superficie hasta el comienzo de la zona convectiva
# con la integración desde el centro hasta la misma zona de transición.

#En primer lugar vamos a convertir las listas en data frames con los valores de la integración desde la superficie
# hasta el comienzo de la zona convectiva. El diccionario que nos almacenaba estos valores lo habíamos llamado "datos".

arriba = pd.DataFrame(datos)[:transicion[1]+1] #Generamos dataframe con los datos hasta la capa de transición

#Convertimos a dataframe también el diccionario con los valores de la integración desde el centro

por_abajo = pd.DataFrame(datos_nucleo)[:100-transicion[1]] #Generamos dataframe con los datos hasta la capa de transición
abajo = por_abajo.iloc[::-1]
abajo.i = 100 - abajo.i

junto = [arriba, abajo]
modelo = pd.concat(junto, sort=False)
modelo.reset_index(drop=True, inplace=True)

#Configuramos ciertas prefrencias para la visualización de los valores:
pd.set_option("display.max_rows",150) #Así logramos que nos muestre todas las filas del data frame
pd.set_option("display.max_columns",10) #Así logramos que nos muestre todas las columnas del data frame
pd.set_option('display.expand_frame_repr', False) #Le pedimos que no nos trunque en otra línea las columnas

#MODELO COMPLETO

#El cálculo en la primera capa comienza para un Rini = 0.9 Rtot por motivos de convergencia.
#Vamos a añadir estas capas que faltan empleando las ecuacionesque usamos para las 3 primeras capas superficiales.

#Lo primero que vamos a hacer es calcular cuántas capas más necesitamos para tener el radio completo:

h = -Rini/100  #Volvemos a definir la variable h dado que la hemos cambiado en la integración desde el centro
#Sabemos que Rtot + i*h = Rini. Despejando para i:

i = -(int((Rini - Rtot)/h)) #Tomamos la parte entera. (No tiene sentido calcular media capa)

#Modificando ligeramente el paso1 obtenemos un algoritmo para calcular estas capas extra:
temp_extra = []
pres_extra = []
masa_extra = []
lum_extra = []
radio_extra = []
dm_extra = []
dl_extra = []
dt_extra = []
dp_extra = []
fase_extra = []
E_extra = []
n_extra = []
i_extra = []

j = -1  #Empezamos a añadir las capas extras como negativas
while j >= i:
    r = Rini + j*h #Ahora el radio va creciendo hacia afuera
    A1 = 1.9022 * mu * Mtot
    A2 = 10.645 * sqrt(Mtot / (mu * Z * (1 + X) * Ltot))
    T = A1 * ((1 / r) - (1 / Rtot))  # Temperatura primeras capas
    P = A2 * T ** 4.25  # Presión primeras capas
    M = Mtot
    L = Ltot

    temp_extra.append(T)
    pres_extra.append(P)
    masa_extra.append(M)
    lum_extra.append(L)
    radio_extra.append(r)

    fm = 0  # En estas primeras capas la consideramos constante
    fl = 0  # En estas primeras capas la consideramos constante
    fp = -Cp * (P / T) * (Mtot / r ** 2)
    ft = -Ct * (P ** 2 / T ** 8.5) * (Ltot / r ** 2)

    dm_extra.append(fm)
    dl_extra.append(fl)
    dp_extra.append(fp)
    dt_extra.append(ft)
    fase_extra.append("^^^^^^")
    infoenergia = funcion_energia(T)
    E_extra.append(infoenergia[5])
    n_extra.append("NaN")  # No usamos el parámetro n+1
    i_extra.append(j)

    j = j - 1


datos_extra = {'E': E_extra, 'fase': fase_extra, 'i': i_extra, 'r': radio_extra, 'Presión': pres_extra, 'Temperatura': temp_extra,'Luminosidad':lum_extra,'Masa':masa_extra,'n+1': n_extra}
extra = pd.DataFrame(datos_extra)
superficiales = extra.iloc[::-1] #Invertimos el data frame para juntarlo con el que ya tenemos

#Generamos el modelo completo
final = [superficiales, modelo]
modelo_completo = pd.concat(final, sort=False)
modelo_completo.reset_index(drop=True, inplace=True)
pd.options.display.float_format = '{:.7f}'.format  #Cambiamos el formato para que se vea como queremos.
pd.set_option('display.width', 1000)
print(modelo_completo)

modelo_completo.reset_index().to_csv('DatosExportados.csv', header=True, index=False) #Exportamos los datos en formato .csv

#JUGANDO CON EL RADIO Y LUMINOSIDAD TOTAL
# Tenemos un modelo completo de una estrella con parámetros constantes Mtot, X, y Y, y con unos valores iniciales para
# Rtot y Ltot. Mediante un ajuste fino hemos ajustado un valor de Tc que minimiza las diferencias entre la integración
# desde la superficie y desde el centro. Dado que la elección de Rtot y Ltot ha sido arbitraria, otra elección
# de estos parametros conducirá a modelos cuyo error relativo total en la transición entre las regiones radiativa
# y convectiva sea aún menor. Vamos a calcular una malla de modelos para diferentes valores de Rtot y Ltot, y
# veremos en cuál de ellos el error relativo total es mínimo.


deltaR = 0.5
deltaL = 5
Rtot_ini = Rtot - 2*deltaR
Rtot_fin = Rtot + 2*deltaR
R = np.linspace(Rtot_ini, Rtot_fin, 5) #Creamos 5 valores del radio total para realizar el cálculo
Ltot_ini = Ltot - 2*deltaL
Ltot_fin = Ltot + 2*deltaL
L = np.linspace(Ltot_ini, Ltot_fin, 5) #Creamos 5 valores de la luminosidad total para realizar el cálculo

#Para cada combinación de Rtot y Ltot, calcularemos el error relativo total. Y así podremos encontrar los valores óptimos.
#Sin embargo, la preción de estos valores está limitada a deltaR y deltaL. Por tanto, cada vez que encontremos
#los valores óptimos, disminuiremos los valores de deltaR y deltaL sucesivamente para obtener un resultado más fino.

k = 1 #Comenzamos la primera iteración:

while k <= 2: #Vamos a pedirle que nos disminuya el intervalo deltaR y deltaL en 3 ocasiones

    print("ITERACIÓN:", k)
    deltaR = (Rtot_fin - Rtot_ini)/5 #Nuevo deltaR que emplearemos para Rtot óptimo
    deltaL = (Ltot_fin - Ltot_ini)/5 #Nuevo deltaL que emplearemos para Ltot óptimo
    matriz = [] #En esta lista almacenaremos listas con los valores de cada fila
    for Ltot in L:
        filas = [] #Lista que almacena los valores de cada fila. Es decir, los valores para cada valor fijo de Ltot iterando sobre todos los Rtot
        for Rtot in R: #Ahora debemos ejecutar todos los cálculos hasta hallar el error relativo total para esta tupla de valores (Ltot, Rtot)

            Tc = Tcc  #Será necesario encontrar también el valor óptimo de Tc para estos (Ltot, Rtot)

            #En principio no modificamos Mtot, ni las fracciones de H y He
            Mtot = 5  # Masa total. Las unidades son 10^33 gramos
            X = 0.75  # Fracción de Hidrógeno
            Y = 0.20  # Fracción de Helio
            Z = 1 - X - Y  # Fracción de Metales
            Rini = 0.9 * Rtot  # No tomamos radio total para evitar problemas de convergencia
            h = -Rini / 100  # Este es el paso de integración. Hemos tomado 100 capas.  Es negativo al integrar desde la #superficie hacia el centro

            # Cálculo del peso molecular medio

            # Supondremos que la estrella es homogénea en composicion química y que el material esta completamente ionizado.

            mu = 1 / (2 * X + 3 / 4 * Y + 1 / 2 * Z)
            Cm = 0.01523 * mu
            Cp = 8.084 * mu
            Ct = 0.01679 * Z * (1 + X) * mu ** 2
            Ct_convec = 3.234 * mu

            # Creamos lista vacía con los valores del radio:
            radio = []

            # Lista con los números de capa
            I = list(range(0, 101))
            temperaturas = []
            presiones = []
            masa = []
            luminosidad = []
            derivadamasa = []
            derivadaluminosidad = []
            derivadapresion = []
            derivadatemperatura = []
            fase = []
            E = []
            n = []  # lista del parametro n+1

            # PRIMERAS 3 CAPAS
            # Calculamos temperatura y presión para las capas exteriores de la estrella
            paso1()
            i = len(presiones) - 1  # Quiero que i tenga el valor de la tercera capa, es decir, i=2
            m_cal = Mtot

            # FASE RADIATIVA A.1.1.
            while abs(Mtot - m_cal) / Mtot < 0.0001:

                r = Rini + (i + 1) * h

                # Paso2
                estimacion = paso2()  # Vector de dos elementos que almacena la presión estimada y la temperatura estimada
                p_est = estimacion[0]
                t_est = estimacion[1]

                # Paso 4  En la fase A.1.1 empleamos como masa de la capa Mtot
                calculaP = paso4(p_est, t_est, Mtot,
                                 r)  # Vector que almacena la derivada de la presión y la presión calculada
                fp = calculaP[0]
                p_cal = calculaP[1]

                # Paso 5
                while abs(
                        p_cal - p_est) / p_cal > 0.0001:  # Verifica si la presión calculada es suficientemente próxima a la estimada.
                    p_est = p_cal  # En caso de no serlo, vuelve al paso 4 y obtiene nuevos valores de la presión
                    calculaP = paso4(p_est, t_est, Mtot, r)  # y su derivada hasta que se cumpla la condición
                    fp = calculaP[0]
                    p_cal = calculaP[1]

                # Paso 7
                calculaT = paso7(p_cal, t_est, Ltot,
                                 r)  # Vector que almacena la derivada de la temperatura y la temperatura calculada
                ft = calculaT[0]
                t_cal = calculaT[1]

                # Paso 8
                while abs(t_cal - t_est) / t_cal > 0.0001:
                    t_est = t_cal
                    calculaP = paso4(p_est, t_est, Mtot, r)  # Volvemos a paso 4
                    fp = calculaP[0]
                    p_cal = calculaP[1]
                    while abs(p_cal - p_est) / p_cal > 0.0001:
                        p_est = p_cal
                        calculaP = paso4(p_est, t_est, Mtot, r)
                        fp = calculaP[0]
                        p_cal = calculaP[1]
                    calculaT = paso7(p_cal, t_est, Ltot,
                                     r)  # Llegamos de nuevo a paso 7 y verifica si se cumple la condición
                    ft = calculaT[0]
                    t_cal = calculaT[1]

                # Paso 3. En la fase A.1.1 empleamos p_cal, t_cal y Mtot
                calculaM = paso3(p_cal, t_cal, Mtot, derivadamasa, r)
                fm = calculaM[0]
                m_cal = calculaM[1]

                # Si la masa no ha cambiado mucho, tomamos como válidos los valores de la capa y asumimos que la masa
                # de la capa es la masa total.
                if abs(Mtot - m_cal) / Mtot < 0.0001:
                    presiones.append(p_cal)
                    temperaturas.append(t_cal)
                    masa.append(Mtot)  # La masa es la masa total
                    luminosidad.append(Ltot)  # La luminosidad es la luminosidad total
                    derivadapresion.append(fp)
                    derivadatemperatura.append(ft)
                    derivadamasa.append(fm)
                    derivadaluminosidad.append(0)  # La luminosidad se considera constante
                    fase.append("A.1.1.")
                    infoenergia = funcion_energia(t_cal)
                    E.append(infoenergia[5])
                    n.append("NaN")  # No usamos aún el parámetro n+1
                    radio.append(r)

                i = i + 1

            # FASE A.1.2.
            # Ya no podemos considerar la masa constante. Emplearemos este algoritmo hasta que la
            # luminosidad deje de considerarse constante

            l_cal = Ltot
            i = len(presiones) - 1

            while abs(Ltot - l_cal) / Ltot < 0.0001:

                r = Rini + (i + 1) * h

                # Paso 2
                estimacion = paso2()  # Vector de dos elementos que almacena la presión estimada y la temperatura estimada
                p_est = estimacion[0]
                t_est = estimacion[1]

                # Paso3. En la fase A.1.2. empleamos p_est, t_est y la masa de la capa i
                calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
                fm = calculaM[0]
                m_cal = calculaM[1]

                # Paso4. En la fase A.1.2 empleamos la masa calculada en la capa i+1
                calculaP = paso4(p_est, t_est, m_cal,
                                 r)  # Vector que almacena la derivada de la presión y la presión calculada
                fp = calculaP[0]
                p_cal = calculaP[1]

                # Paso 5
                while abs(
                        p_cal - p_est) / p_cal > 0.0001:  # Verifica si la presión calculada es suficientemente próxima a la estimada.
                    p_est = p_cal  # En caso de no serlo, vuelve al paso 3 y obtiene nuevos valores de la masa,
                    calculaM = paso3(p_est, t_est, masa[i], derivadamasa,
                                     r)  # de la presión y sus derivadas hasta que se cumpla la condición
                    fm = calculaM[0]
                    m_cal = calculaM[1]
                    calculaP = paso4(p_est, t_est, m_cal, r)
                    fp = calculaP[0]
                    p_cal = calculaP[1]

                # Paso 7. Utilizamos el valor de Ltot, pues consideramos que permanece constante
                calculaT = paso7(p_cal, t_est, Ltot,
                                 r)  # Vector que almacena la derivada de la temperatura y la temperatura calculada
                ft = calculaT[0]
                t_cal = calculaT[1]

                # Paso 8
                while abs(t_cal - t_est) / t_cal > 0.0001:
                    t_est = t_cal
                    calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # Volvemos al paso 3
                    fm = calculaM[0]
                    m_cal = calculaM[1]
                    calculaP = paso4(p_est, t_est, m_cal, r)  # Pasamos al 4
                    fp = calculaP[0]
                    p_cal = calculaP[1]
                    while abs(p_cal - p_est) / p_cal > 0.0001:  # Verifica si p_cal próximo a p_est
                        p_est = p_cal
                        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # Paso 3
                        fm = calculaM[0]
                        m_cal = calculaM[1]
                        calculaP = paso4(p_est, t_est, m_cal, r)
                        fp = calculaP[0]
                        p_cal = calculaP[1]
                    calculaT = paso7(p_cal, t_est, Ltot,
                                     r)  # Llegamos de nuevo a paso 7 y verifica si se cumple la condición
                    ft = calculaT[0]
                    t_cal = calculaT[1]

                # Paso 6 Calculamos la luminosidad y su derivada en la capa i+1.
                # En la fase A.1.2. empleamos t_cal.
                # Como la luminosidad se considera constante, la luminosidad de la capa anterior es Ltot
                calculaL = paso6(p_cal, t_cal, Ltot, derivadaluminosidad, r)
                fl = calculaL[0]
                l_cal = calculaL[1]

                # Si la luminosidad no ha cambiado mucho, tomamos como válidos los valores de la capa y asumimos que la luminosidad
                # de la capa es la luminosidad total
                if abs(Ltot - l_cal) / Ltot < 0.0001:
                    presiones.append(p_cal)
                    temperaturas.append(t_cal)
                    masa.append(m_cal)
                    luminosidad.append(Ltot)  # La luminosidad es la luminosidad total
                    derivadapresion.append(fp)
                    derivadatemperatura.append(ft)
                    derivadamasa.append(fm)
                    derivadaluminosidad.append(fl)
                    fase.append("A.1.2.")
                    infoenergia = funcion_energia(t_cal)
                    E.append(infoenergia[5])
                    n.append("NaN")  # No usamos aún el parámetro n+1
                    radio.append(r)

                i = i + 1

            #FASE A.1.3.

            # Ahora tanto la luminosidad como la masa varían en cada capa.
            i = len(presiones) - 1
            a = 3  # Esta variable es n+1. Hacemos que tome un valor mayor a 2.5 para que comience el bucle
            while a > 2.5:

                r = Rini + (i + 1) * h

                # Paso 2
                estimacion = paso2()  # Vector de dos elementos que almacena la presión estimada y la temperatura estimada
                p_est = estimacion[0]
                t_est = estimacion[1]
                # Paso3. En la fase A.1.3. empleamos p_est, t_est y la masa de la capa i
                calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
                fm = calculaM[0]
                m_cal = calculaM[1]

                # Paso4. En la fase A.1.3 empleamos la masa calculada en la capa i+1
                calculaP = paso4(p_est, t_est, m_cal,
                                 r)  # Vector que almacena la derivada de la presión y la presión calculada
                fp = calculaP[0]
                p_cal = calculaP[1]

                # Paso 5
                while abs(
                        p_cal - p_est) / p_cal > 0.0001:  # Verifica si la presión calculada es suficientemente próxima a la estimada.
                    p_est = p_cal  # En caso de no serlo, vuelve al paso 3 y obtiene nuevos valores de la masa,
                    calculaM = paso3(p_est, t_est, masa[i], derivadamasa,
                                     r)  # de la presión y sus derivadas hasta que se cumpla la condición
                    fm = calculaM[0]
                    m_cal = calculaM[1]
                    calculaP = paso4(p_est, t_est, m_cal, r)
                    fp = calculaP[0]
                    p_cal = calculaP[1]

                # Paso 6. Una vez que p_cal = p_est, pasamos a calcular la luminosidad.
                # Esta vez empleamos T_est, pues t_cal se hallará en el paso 7.
                # La luminosidad varía y ya no podemos utilizar Ltot.
                calculaL = paso6(p_cal, t_est, luminosidad[i], derivadaluminosidad, r)
                fl = calculaL[0]
                l_cal = calculaL[1]

                # Paso 7. Dado que ahora la luminosidad no es constante, debemos emplear la luminosidad calculada en el
                # paso anterior para calcular la temperatura.
                calculaT = paso7(p_cal, t_est, l_cal,
                                 r)  # Vector que almacena la derivada de la temperatura y la temperatura calculada
                ft = calculaT[0]
                t_cal = calculaT[1]

                # Paso 8
                while abs(t_cal - t_est) / t_cal > 0.0001:
                    t_est = t_cal
                    calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # Volvemos al paso 3
                    fm = calculaM[0]
                    m_cal = calculaM[1]
                    calculaP = paso4(p_est, t_est, m_cal, r)  # Pasamos al 4
                    fp = calculaP[0]
                    p_cal = calculaP[1]
                    while abs(p_cal - p_est) / p_cal > 0.0001:  # Verifica si p_cal próximo a p_est
                        p_est = p_cal
                        calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)  # Paso 3
                        fm = calculaM[0]
                        m_cal = calculaM[1]
                        calculaP = paso4(p_est, t_est, m_cal, r)
                        fp = calculaP[0]
                        p_cal = calculaP[1]
                    calculaL = paso6(p_cal, t_est, luminosidad[i], derivadaluminosidad, r)
                    fl = calculaL[0]
                    l_cal = calculaL[1]
                    calculaT = paso7(p_cal, t_est, l_cal,
                                     r)  # Llegamos de nuevo a paso 7 y verifica si se cumple la condición
                    ft = calculaT[0]
                    t_cal = calculaT[1]

                # Finalmente, debemos comprobar si o la hipótesis inicial de transporte radiativo deja de ser válida,
                # vamos a calcular el párametro n+1 definido en el Novotny

                # Llamamos a=n+1

                a = t_cal * h * fp / (p_cal * h * ft)

                if a > 2.5:
                    presiones.append(p_cal)
                    temperaturas.append(t_cal)
                    masa.append(m_cal)
                    luminosidad.append(l_cal)
                    derivadapresion.append(fp)
                    derivadatemperatura.append(ft)
                    derivadamasa.append(fm)
                    derivadaluminosidad.append(fl)
                    fase.append("A.1.3.")
                    infoenergia = funcion_energia(t_cal)
                    E.append(infoenergia[5])
                    n.append(a)
                    radio.append(r)

                else:  # Cálculos conducen a n+1 ≤ 2.5 y datos de la capa no son válidos
                    t_convec = t_cal
                    p_convec = p_cal
                    n.append(a)  # añadimos el valor del parámetro n+1 de la primera capa convectiva

                r_transicion = Rini + i * h  # Guardamos el valor del radio para el que se produce la transición entre fase radiativa y convectiva
                transicion = [r_transicion, i]

                i = i + 1

            # FASE A.2. NUCLEO CONVECTIVO
            # Cuando la aplicación de la fase A.1.3 conduce a un cálculo del parámetro n+1 ≤ 2.5 sabemos que los cálculos
            # realizados en dicha capa no son válidos y debemos repetirlos empleando el algoritmo A.2

            # No obstante, los valores de presión y temperatura calculados en esa capa sí que nos sirven para
            # estimar la constante K que nos relaciona estos dos parametros fíısicos en un polítropo.
            # Asumiendo un índice adiabático γ = 5/3 (correspondiente a un gas perfecto monoatomico):

            K = p_convec / t_convec ** 2.5  # Este valor calculado en la primera capa convectiva se asumirá constante

            i = len(presiones) - 1

            while r >= 0:

                r = Rini + (i + 1) * h

                # Paso 2
                t_est = paso2bis(temperaturas, derivadatemperatura)
                # Estimamos la presión con expresión del polítropo

                if t_est > 0: #Para radios muy pequeños cercanos a cero los errores se acumulan y podemos obtener resultados absurdos
                    p_est = K * t_est ** 2.5
                else:
                    p_est = 0

                # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
                # expresión del polítropo.
                calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
                fm = calculaM[0]
                m_cal = calculaM[1]

                # No hay pasos 4 y 5 porque la presión ahora se calcula directamente con la expresión del polítropo.
                # Dado que el gradiente de temperatura en el caso convectivo no depende de la luminosidad, el paso 6 se
                # incorporará después del paso 8.

                # Paso 7bis
                calculaT = paso7bis(m_cal, temperaturas, derivadatemperatura, r)
                ft = calculaT[0]
                t_cal = calculaT[1]

                # Paso 8
                while abs(t_cal - t_est) / t_cal > 0.0001:
                    t_est = t_cal
                    p_est = K * t_est ** 2.5  # Volvemos a estimar la presión
                    # Repetimos paso 3
                    calculaM = paso3(p_est, t_est, masa[i], derivadamasa, r)
                    fm = calculaM[0]
                    m_cal = calculaM[1]
                    # Repetimos paso 7bis
                    calculaT = paso7bis(m_cal, temperaturas, derivadatemperatura, r)
                    ft = calculaT[0]
                    t_cal = calculaT[
                        1]  # En este punto verifica si t_cal es prácticamente igual a t_est. De ser así sale del bucle

                # Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
                if r > 0 and t_cal > 0:
                    p_cal = K * t_cal ** 2.5

                else:  # En el caso del centro de la estrella basta con utilizar p_cal = p_est y t_cal = t_est
                    t_cal = paso2bis(temperaturas, derivadatemperatura)
                    if t_cal > 0:
                        p_cal = K * t_cal ** 2.5
                    else:
                        p_cal = 0

                    # Paso 3
                    calculaM = paso3(p_cal, t_cal, masa[i], derivadamasa, r)
                    fm = calculaM[0]
                    m_cal = calculaM[1]

                # Paso 6. La ecuación para el cálculo de la luminosidad es igual que en el caso radiativo, pero considerando la
                # expresión del polítropo
                calculaL = paso6(p_cal, t_cal, luminosidad[i], derivadaluminosidad, r)
                fl = calculaL[0]
                l_cal = calculaL[1]

                # Si los cálculos de la capa son válidos, los añadimos a nuestras listas
                if r >= 0:
                    presiones.append(p_cal)
                    temperaturas.append(t_cal)
                    masa.append(m_cal)
                    luminosidad.append(l_cal)
                    derivadatemperatura.append(ft)
                    derivadamasa.append(fm)
                    derivadaluminosidad.append(fl)
                    fase.append("CONVEC")
                    infoenergia = funcion_energia(t_cal)
                    E.append(infoenergia[5])
                    n.append("NaN")  # No usamos el parámetro n+1
                    radio.append(r)

                i = i + 1

            # VALORES EN LA FRONTERA (integración desde arriba hacia el centro de la estrella)
            #Realizamos interpolación lineal para hallar el radio en que n+1 = 2.5
            r_rad_conv = transicion[0] + (radio[transicion[1] + 1] - transicion[0]) / (
                    n[transicion[1] + 1] - n[transicion[1]]) * \
                         (2.5 - n[transicion[1]])
            # Interpolamos linealmente el resto de parámetros físicos también:
            p_rad_conv = presiones[transicion[1]] + \
                         (presiones[transicion[1] + 1] - presiones[transicion[1]]) / (
                                 n[transicion[1] + 1] - n[transicion[1]]) * \
                         (2.5 - n[transicion[1]])
            t_rad_conv = temperaturas[transicion[1]] + \
                         (temperaturas[transicion[1] + 1] - temperaturas[transicion[1]]) / (
                                 n[transicion[1] + 1] - n[transicion[1]]) \
                         * (2.5 - n[transicion[1]])
            l_rad_conv = luminosidad[transicion[1]] + \
                         (luminosidad[transicion[1] + 1] - luminosidad[transicion[1]]) / (
                                 n[transicion[1] + 1] - n[transicion[1]]) \
                         * (2.5 - n[transicion[1]])
            m_rad_conv = masa[transicion[1]] + \
                         (masa[transicion[1] + 1] - masa[transicion[1]]) / (n[transicion[1] + 1] - n[transicion[1]]) * \
                         (2.5 - n[transicion[1]])

            #INTEGRACIÓN DESDE EL CENTRO
            h = Rini / 100  # Paso de integración positivo al integrar desde el centro hacia la superficie

            #PRIMERAS 3 CAPAS
            masa_nucleo = []
            luminosidad_nucleo = []
            temperatura_nucleo = []
            presion_nucleo = []
            fase_nucleo = []
            E_nucleo = []
            I_nucleo = []
            radio_nucleo = []
            dm_nucleo = []
            dl_nucleo = []
            dt_nucleo = []
            dp_nucleo = []

            for i in range(0, 3):
                r = 0 + i * h
                parametros = parametros_nucleo(r, Tc)
                m = parametros[0]
                l = parametros[1]
                t = parametros[2]
                p = parametros[3]
                dm = parametros[4]
                dl = parametros[5]
                dt = parametros[6]
                dp = parametros[7]
                masa_nucleo.append(m)
                luminosidad_nucleo.append(l)
                temperatura_nucleo.append(t)
                presion_nucleo.append(p)
                info_energia = funcion_energia(t)
                E_nucleo.append(info_energia[5])
                fase_nucleo.append("CENTRO")
                radio_nucleo.append(r)
                dm_nucleo.append(dm)
                dl_nucleo.append(dl)
                dt_nucleo.append(dt)
                dp_nucleo.append(dp)
                I_nucleo.append(i)

            # CAPAS POSTERIORES

            i = len(presion_nucleo) - 1

            while r <= transicion[0]:

                r = 0 + (i + 1) * h

                # Paso 1. Ya conocemos las presiones, temperaturas, masas, luminosidades y sus derivadas de las tres primeras capas

                # Paso 2bis. Estimamos la temperatura de la capa siguiente:
                t_est = paso2bis(temperatura_nucleo, dt_nucleo)

                # Estimamos la presión con expresión del polítropo
                p_est = K * t_est ** 2.5

                # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
                # expresión del polítropo.
                calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
                fm = calculaM[0]
                m_cal = calculaM[1]

                # Paso 7bis. A continuación podemos calcular la temperatura:
                calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
                ft = calculaT[0]
                t_cal = calculaT[1]

                # Paso 8. Comparamos la temperatura calculada con la estimada:
                while abs(t_cal - t_est) / t_cal > 0.0001:
                    t_est = t_cal
                    p_est = K * t_est ** 2.5  # Volvemos a estimar la presión
                    # Repetimos paso 3
                    calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
                    fm = calculaM[0]
                    m_cal = calculaM[1]
                    # Repetimos paso 7bis
                    calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
                    ft = calculaT[0]
                    t_cal = calculaT[
                        1]  # En este punto verifica si t_cal es prácticamente igual a t_est. De ser así, sale del bucle

                # Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
                p_cal = K * t_cal ** 2.5

                # Paso 6. Conocida p_cal podemos hallar la luminosidad:
                calculaL = paso6(p_cal, t_cal, luminosidad_nucleo[i], dl_nucleo, r)
                fl = calculaL[0]
                l_cal = calculaL[1]

                if r <= transicion[0]:  # Si los cálculos de la capa son válidos, los añadimos
                    presion_nucleo.append(p_cal)
                    temperatura_nucleo.append(t_cal)
                    masa_nucleo.append(m_cal)
                    luminosidad_nucleo.append(l_cal)
                    dt_nucleo.append(ft)
                    dm_nucleo.append(fm)
                    dl_nucleo.append(fl)
                    fase_nucleo.append("CONVEC")
                    infoenergia = funcion_energia(t_cal)
                    E_nucleo.append(infoenergia[5])
                    radio_nucleo.append(r)
                    I_nucleo.append(i + 1)
                    transicion2 = [r, i + 1]  # Guardamos el radio en el que comienza la fase radiativa

                i = i + 1

            # VALORES EN LA FRONTERA (integracion desde el centro hacia el exterior de la estrella)
            r_conv_rad = radio_nucleo[transicion2[1] - 1] + (transicion2[0] - radio_nucleo[transicion2[1] - 1]) / \
                         (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])
            p_conv_rad = presion_nucleo[transicion2[1] - 1] + \
                         (presion_nucleo[transicion2[1]] - presion_nucleo[transicion2[1] - 1]) / \
                         (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])
            t_conv_rad = temperatura_nucleo[transicion2[1] - 1] + \
                         (temperatura_nucleo[transicion2[1]] - temperatura_nucleo[transicion2[1] - 1]) / \
                         (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])
            l_conv_rad = luminosidad_nucleo[transicion2[1] - 1] + \
                         (luminosidad_nucleo[transicion2[1]] - luminosidad_nucleo[transicion2[1] - 1]) / \
                         (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])
            m_conv_rad = masa_nucleo[transicion2[1] - 1] + \
                         (masa_nucleo[transicion2[1]] - masa_nucleo[transicion2[1] - 1]) / \
                         (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

            # Podemos ahora comparar las soluciones que hemos obtenido para r=rconv→rad=rrad→conv
            error_rel_p = (p_rad_conv - p_conv_rad) / p_rad_conv * 100
            error_rel_t = (t_rad_conv - t_conv_rad) / t_rad_conv * 100
            error_rel_l = (l_rad_conv - l_conv_rad) / l_rad_conv * 100
            error_rel_m = (m_rad_conv - m_conv_rad) / m_rad_conv * 100
            error_rel_tot = sqrt(error_rel_p ** 2 + error_rel_t ** 2 + error_rel_l ** 2 + error_rel_m ** 2)

            errores = {'Error.relat P(%)': error_rel_p, 'Error.relat T(%)': error_rel_t, 'Error.relat L(%)': error_rel_l, 'Error.relat M(%)': error_rel_m}

            Tcini = 1
            Tcfin = 3

            Tc = np.linspace(Tcini, Tcfin, 50)  # Hacemos que tome 50 valores entre 1 y 3
            Errores = []  # En esta lista almacenaremos los errores relativos totales para cada Tc
            Lista_de_diccionarios = []  # Esta lista será una lista de diccionarios

            iteraciones = 1  # Este contador mide el número de iteraciones que realizamos.
            # En cada iteración obtenemos un valor más preciso de Tc

            while iteraciones <= 3:  # Vamos a pedir que haga tres iteraciones. Podemos cambiar libremente este número

                intervalo = (Tcfin - Tcini) / 50  # Este es el paso en Tc. (Nos será útil para afinar aún más el valor de Tc)

                for temperatura in Tc:

                    masa_nucleo = []
                    luminosidad_nucleo = []
                    temperatura_nucleo = []
                    presion_nucleo = []
                    fase_nucleo = []
                    E_nucleo = []
                    I_nucleo = []
                    radio_nucleo = []
                    dm_nucleo = []
                    dl_nucleo = []
                    dt_nucleo = []
                    dp_nucleo = []

                    # PRIMERAS 3 CAPAS
                    for i in range(0, 3):
                        r = 0 + i * h
                        parametros = parametros_nucleo(r, temperatura)
                        m = parametros[0]
                        l = parametros[1]
                        t = parametros[2]
                        p = parametros[3]
                        dm = parametros[4]
                        dl = parametros[5]
                        dt = parametros[6]
                        dp = parametros[7]
                        masa_nucleo.append(m)
                        luminosidad_nucleo.append(l)
                        temperatura_nucleo.append(t)
                        presion_nucleo.append(p)
                        info_energia = funcion_energia(t)
                        E_nucleo.append(info_energia[5])
                        fase_nucleo.append("CENTRO")
                        radio_nucleo.append(r)
                        dm_nucleo.append(dm)
                        dl_nucleo.append(dl)
                        dt_nucleo.append(dt)
                        dp_nucleo.append(dp)
                        I_nucleo.append(i)

                    # CAPAS POSTERIORES

                    i = len(presion_nucleo) - 1

                    while r <= transicion[0]:

                        r = 0 + (i + 1) * h

                        # Paso 1. Ya conocemos las presiones, temperaturas, masas, luminosidades y sus derivadas de las tres primeras capas

                        # Paso 2bis. Estimamos la temperatura de la capa siguiente:
                        t_est = paso2bis(temperatura_nucleo, dt_nucleo)

                        # Estimamos la presión con expresión del polítropo
                        p_est = K * t_est ** 2.5

                        # Paso 3. En la fase A.2. la ecuación de la masa es igual que en el caso radiativo, pero considerando la
                        # expresión del polítropo.
                        calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
                        fm = calculaM[0]
                        m_cal = calculaM[1]

                        # Paso 7bis. A continuación podemos calcular la temperatura:
                        calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
                        ft = calculaT[0]
                        t_cal = calculaT[1]

                        # Paso 8. Comparamos la temperatura calculada con la estimada:
                        while abs(t_cal - t_est) / t_cal > 0.0001:
                            t_est = t_cal
                            p_est = K * t_est ** 2.5  # Volvemos a estimar la presión
                            # Repetimos paso 3
                            calculaM = paso3(p_est, t_est, masa_nucleo[i], dm_nucleo, r)
                            fm = calculaM[0]
                            m_cal = calculaM[1]
                            # Repetimos paso 7bis
                            calculaT = paso7bis(m_cal, temperatura_nucleo, dt_nucleo, r)
                            ft = calculaT[0]
                            t_cal = calculaT[
                                1]  # En este punto verifica si t_cal es prácticamente igual a t_est. De ser así, sale del bucle

                        # Ahora que disponemos de t_cal, calculamos el valor de la presión con la expresión del polítropo.
                        p_cal = K * t_cal ** 2.5

                        # Paso 6. Conocida p_cal podemos hallar la luminosidad:
                        calculaL = paso6(p_cal, t_cal, luminosidad_nucleo[i], dl_nucleo, r)
                        fl = calculaL[0]
                        l_cal = calculaL[1]

                        if r <= transicion[0]:  # Si los cálculos de la capa son válidos, los añadimos
                            presion_nucleo.append(p_cal)
                            temperatura_nucleo.append(t_cal)
                            masa_nucleo.append(m_cal)
                            luminosidad_nucleo.append(l_cal)
                            dt_nucleo.append(ft)
                            dm_nucleo.append(fm)
                            dl_nucleo.append(fl)
                            fase_nucleo.append("CONVEC")
                            infoenergia = funcion_energia(t_cal)
                            E_nucleo.append(infoenergia[5])
                            radio_nucleo.append(r)
                            I_nucleo.append(i + 1)
                            transicion2 = [r, i + 1]  # Guardamos el radio en el que comienza la fase radiativa

                        i = i + 1

                    r_conv_rad = radio_nucleo[transicion2[1] - 1] + (
                                transicion2[0] - radio_nucleo[transicion2[1] - 1]) / \
                                 (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

                    # Interpolamos linealmente el resto de parámetros físicos también:
                    p_conv_rad = presion_nucleo[transicion2[1] - 1] + \
                                 (presion_nucleo[transicion2[1]] - presion_nucleo[transicion2[1] - 1]) / \
                                 (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

                    t_conv_rad = temperatura_nucleo[transicion2[1] - 1] + \
                                 (temperatura_nucleo[transicion2[1]] - temperatura_nucleo[transicion2[1] - 1]) / \
                                 (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

                    l_conv_rad = luminosidad_nucleo[transicion2[1] - 1] + \
                                 (luminosidad_nucleo[transicion2[1]] - luminosidad_nucleo[transicion2[1] - 1]) / \
                                 (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

                    m_conv_rad = masa_nucleo[transicion2[1] - 1] + \
                                 (masa_nucleo[transicion2[1]] - masa_nucleo[transicion2[1] - 1]) / \
                                 (n[transicion[1]] - n[transicion[1] + 1]) * (2.5 - n[transicion[1] + 1])

                    # Calculamos los errores para esta temperatura central
                    error_rel_p = (p_rad_conv - p_conv_rad) / p_rad_conv * 100
                    error_rel_t = (t_rad_conv - t_conv_rad) / t_rad_conv * 100
                    error_rel_l = (l_rad_conv - l_conv_rad) / l_rad_conv * 100
                    error_rel_m = (m_rad_conv - m_conv_rad) / m_rad_conv * 100
                    error_rel_tot = sqrt(error_rel_p ** 2 + error_rel_t ** 2 + error_rel_l ** 2 + error_rel_m ** 2)

                    Errores.append(error_rel_tot)

                data = {"Tc": Tc, "Err_rel": Errores}
                df = pd.DataFrame(data)
                optimo = df[df.Err_rel == df.Err_rel.min()]  # data frame con la Tc que hace mínimo el error
                Toptima = optimo.iloc[0, 0]  # Extraemos el valor óptimo de la temperatura central

                # Ahora ajustamos el intervalo para obtener un mejor valor de Tc:

                Tcini = Toptima - intervalo
                Tcfin = Toptima + intervalo

                Tc = np.linspace(Tcini, Tcfin, 50)  # Hacemos que tome 50 valores entre Tcini y Tcfin
                Errores = []  # Limpiamos la lista de errores para la siguiente iteracion

                iteraciones = iteraciones + 1

            filas.append(optimo.iloc[0, 1])
        matriz.append(filas)

    row_labels = L
    column_labels = R
    df = pd.DataFrame(matriz, columns=column_labels, index=row_labels)
    print(df)
    ax = plt.axes()
    mapa_de_color = sns.heatmap(df, vmax=30)
    ax.set_ylabel('Luminosidad total')
    ax.set_xlabel('Radio total')
    ax.set_title('Mapa de color de los errores relativos')
    plt.show()

    #Ahora tenemos un data frame que contiene los errores relativos para cada Rtot y Ltot.
    #Buscamos la celda con el menor relativo.
    row_index = np.unravel_index(np.argmin(df.values), df.shape)[0] #Indice de la fila donde se encuentra el valor mínimo
    column_index = np.unravel_index(np.argmin(df.values), df.shape)[1] #Indice de la columna donde se encuentra el valor mínimo

    #Valores óptimos de Ltot y Rtot:
    Loptimo = df.index[row_index]
    Roptimo = df.columns[column_index]
    print("Loptimo es:",Loptimo,"Roptimo es:", Roptimo)

    #Ajustamos el intervalo para obtener mejores valores:
    Rtot_ini = Roptimo - deltaR
    Rtot_fin = Roptimo + deltaR
    Ltot_ini = Loptimo - deltaL
    Ltot_fin = Loptimo + deltaL

    R = np.linspace(Rtot_ini, Rtot_fin, 5)
    L = np.linspace(Ltot_ini, Ltot_fin, 5)

    print("Fin iteración:", k)
    print("------------------")

    k = k + 1

ax = plt.axes()
mapa_de_color = sns.heatmap(df,  vmax=30)
ax.set_ylabel('Luminosidad total')
ax.set_xlabel('Radio total')
ax.set_title('Mapa de color de los errores relativos')
plt.show()

#-----------------------------------------------------------------------------#
#---------------------------FIN DEL PROGRAMA----------------------------