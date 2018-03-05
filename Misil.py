# -*- coding: utf-8 -*-
"""
@author: Team REOS
"""
from math import exp, radians, cos, pi, sin, degrees, atan, log10

'''
Este programa servirá para establecer la trayectoria a describir por el
avión.  Contiene la maniobra de giro, seguida por otra de ascenso.  Al final se
procederá a exportar las características y cómo varían a lo largo de la
trayectoria para poder obtener sus gráficas.
'''

'''Condiciones g0itatorias y constantes atmosféricas'''
G = 6.673e-11  # Constante de g0itación.
MT = 5.972e24  # Masa terrestre.
MU = G * MT
RT = 6378136.3  # Radio terrestre.
R_AIR = 287  # Constante de los gases ideales.
g0 = 9.80665  # Aceleración gravitatoria inicial
RHO_SL = 101325 / (R_AIR * 288.15)  # Densidad del aire a nivel del mar.
GAMMA = 1.4  # Coeficiente de dilatación adiabática.

##Condiciones de viscosidad
beta_visc = 0.000001458
S_visc = 110.4

'''
-----------------------ATMÓSFERA ESTÁNDAR INTERNACIONAL-----------------------
Esta subrutina nos permitirá obtener los valores de presión, temperatura y
densidad del aire en funcióin de la altura. Están sacados de la ISA.  Más
adelante, en los siguientes bucles, se llamará a las funciones de P, T y rho
que variarán con respecto a h.
'''

#Estas son las alturas estipuladas según la normativa.

H_ISA1 = 11000
H_ISA2 = 20000
H_ISA3 = 32000
H_ISA4 = 47000
H_ISA5 = 51000
H_ISA6 = 71000
H_ISA7 = 84852

  

#Ahora se programan las variables termodinámicas, en función de la altura,
# y se relacionarán con los valores de T y alfa para cada altura estipulada.

def temperature(alt):
    '''Cálculo de la temperatura en función de la altura dada por el modelo ISA
    '''
    if alt < H_ISA1:
        h_0 = 0
        t_0 = 288.15
        alfa_isa = -.0065
    elif alt < H_ISA2:
        h_0 = H_ISA1
        t_0 = 216.65
        alfa_isa = 0
    elif alt < H_ISA3:
        h_0 = H_ISA2
        t_0 = 216.65
        alfa_isa = .001
    elif alt < H_ISA4:
        h_0 = H_ISA3
        t_0 = 228.65
        alfa_isa = .0028
    elif alt < H_ISA5:
        h_0 = H_ISA4
        t_0 = 270.65
        alfa_isa = 0
    elif alt < H_ISA6:
        h_0 = H_ISA5
        t_0 = 270.65
        alfa_isa = -.0028
    elif alt < H_ISA7:
        h_0 = H_ISA6
        t_0 = 214.65
        alfa_isa = -.002
    else:
        h_0 = H_ISA7
        t_0 = 214.65 - .002 * (H_ISA7 - H_ISA6)
        alfa_isa = 0
    return t_0 + alfa_isa * (alt - h_0)

def density(alt):
    '''Cálculo de la densidad en función de la altura dada por el modelo ISA
    '''
    rho0 = RHO_SL
    t_isa = temperature(alt)
    h_0 = 0
    t_0 = 288.15
    alfa_isa = -.0065
    if alt < H_ISA1:
    
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA1) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA1
    if alt < H_ISA2:
        
        return rho0 * exp(-g0 * (alt - h_0) / (R_AIR * t_isa))
    rho0 = rho0 * exp(-g0 * (H_ISA2 - h_0) / (R_AIR * temperature(H_ISA2)))
    h_0 = H_ISA2
    t_0 = 216.65
    alfa_isa = .001
    if alt < H_ISA3:
        
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA3) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA3
    t_0 = 228.65
    alfa_isa = .0028
    if alt < H_ISA4:
        
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA4) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA4
    if alt < H_ISA5:
        
        return rho0 * exp(-g0 * (alt - h_0)/(R_AIR * t_isa))
    rho0 = rho0 * exp(-g0 * (H_ISA5 - h_0)/(R_AIR * temperature(H_ISA5)))
    h_0 = H_ISA5
    t_0 = 270.65
    alfa_isa = -.0028
    if alt < H_ISA6:
        
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA6) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA6
    t_0 = 214.65
    alfa_isa = -.002
    if alt < H_ISA7:
       
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA7) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA7
    return rho0 * exp(-g0 * (alt - h_0) / (R_AIR * t_isa))

def pressure(alt):
    '''Cálculo de la presión en función de la altura dada por el modelo ISA
    '''
    return density(alt) * R_AIR * temperature(alt)

def viscosity(alt):	
    '''Cálculo de la viscosidad en función de la altura dada por el modelo ISA
    '''							
    if alt < H_ISA1:
	    h_0 = 0
	    t_0 = 288.15
	    alfa_isa = -0.0065
	    t = t_0 + alfa_isa * (alt- h_0)
     
    elif alt< H_ISA2:
      
	    h_0 = H_ISA1
	    t_0 = 216.65   
	    alfa_isa = 0
	    t = t_0 + alfa_isa * (alt- h_0)  
		 
    elif alt< H_ISA3:
   
	    h_0 = H_ISA2
	    t_0 = 216.65
	    alfa_isa = 0.001
		   
	    t = t_0 + alfa_isa * (alt- h_0) 
                    
    elif alt< H_ISA4:
          
        h_0 = H_ISA3
        t_0 = 228.65                          
        alfa_isa = 0.0028
     
        t = t_0 + alfa_isa * (alt- h_0)                           
                          
    elif alt< H_ISA5:
        
        h_0 = H_ISA4
        t_0 = 270.65                                                               
        alfa_isa = 0                                 
        t = t_0 + alfa_isa * (alt- h_0)                                  
    
    elif alt< H_ISA6:
       
        h_0 = H_ISA5
        t_0 = 270.65                                     
        alfa_isa = -0.0028
    
        t = t_0 + alfa_isa * (alt- h_0)                                         
    
    elif alt< H_ISA7:
    
        h_0 = H_ISA6
        t_0 = 214.65                                               
        alfa_isa = -0.002
                                               
        t = t_0 + alfa_isa * (alt- h_0)                                         
        
    elif alt> H_ISA7: 
                                                                                             
	    h_0 = H_ISA7
	    t_0 = 214.65 - 0.002 * (H_ISA7 - H_ISA6)													                                                    
	    alfa_isa = 0
	    t = t_0 + alfa_isa * (alt- h_0)                                                    
		
    return (beta_visc * t ** (3 / 2)) / (t + S_visc)   

'''
-------------------CARACTERÍSTICAS GEOMÉTRICAS DEL VEHÍCULO-------------------
'''

TH_SL = 100000  # Empuje a nivel del mar (máximo).
N = 3.5  # Factor de carga máximo.
W = 14273 * g0  # Peso del avión (N).
S_W = 49.2  # Superficie alar (m2).
B = 11.7  # Envergadura (m).
AR = B**2 / S_W  # Alargamiento = 2,78.
FLECHA = radians(52)  # Flecha BA.
FLECHA2 = radians(41.4)  # Flecha 1/4.
LF = 19.2  # Longitud del fuselaje.
BF = 2.87  # Longitud del fuselaje.
#Se desrecian las pérdidas por consumo de combustible en el peso del avión.

K = .4  # EL perfil del F-4 es el NACA0006.4-64 por tanto la K es 0,4.
ESTRECHAMIENTO = .26  # Estrechamiento.
ESPESOR = .064  # Espesor relativo máximo.

NE = 2  # Número de motores.
DM = 0  # Diámetro de cada motor.
LM = 0  # Longitud de cada motor.
#La longitud y el diámetro de cada motor aparecen nulos porque no se ha
# refelejado todavía en los cálculos la importancia de la geometría de los
# motores.

MASS = 14273  # Masa de carga

#A continuación, se definen las superficies de mando del avión, que nos
# servirán para, más adelante, calcular el coeficiente de resistencia parásita.
S_LEX = 0  # Área del Lex. [¡!]
S_H = 6.39  # Área de la superficie de mando horizontal
S_V = 5.035  # Área de la superficie de mando vertical

#Se necsitarán para más adelante los valores de Mach crítico y de divergencia,
# los cuales son función de la flecha y del espesor relativo máximo.  Estos
# valores marcarán los límites de los dominios que nos servirán para calcular
# el coeficicente de resitencia inducida con efecto de compresibilidad.
M_C = 1 - (.065 * (cos(FLECHA))**.6) * (100 * ESPESOR)**.637814  # Mach crítico
M_D = 1.08 * M_C  # Mach de divergencia
M_D098 = .98 * M_D


'''
-------------------CARACTERÍSTICAS GEOMÉTRICAS DEL MISIL-------------------
'''

diametro_m=0.5           #diametro misil
longitud_cono = 0.9      #longitud cono misil
longitud_misil = 3       #longitud total misil


Sw_aleta = 0.07875                           #Superficie de una aleta del AIM (tomado como ref.)
Craiz_aleta = 0.24                           #cuerda raiz de aleta
Cmedia_aleta = 0.18                          #cuerda media de aleta
espesor_aleta = 0.0065                       #espesor de aleta
tao_aleta = espesor_aleta/Cmedia_aleta       #tao aleta
num_aletas = 4                               #numero aletas
 
Swtotal_aletas=Sw_aleta*num_aletas           #Superficie total de aletas
Sref_aletas = Swtotal_aletas/2               #Superficie de referencia aletas


Sup_cono=pi*(diametro_m/2)*(longitud_cono**2+(diametro_m/2)**2)**(1/2)
Sup_total=2*pi*(diametro_m/2)*(longitud_misil-longitud_cono)
Sref_misil=pi*(diametro_m/2)**2
Sgases=pi*((diametro_m)*0.9/2)**2 #Area salida de los gases (consideramos el área de salida de la tobera)
Ratio_areas=(Sref_misil-Sgases)/Sref_misil #Relación de áreas 

angulo_cono = atan(0.5*diametro_m/longitud_cono)*(180/pi)  #angulo de la generatriz del cono

#gasto=30                      #Empuje constante del misil
#masa_misil=1000
#masa_propulsante=680
#Isp = 254.84*9.81
#Empuje_misil= gasto*Isp
#t_combustion= masa_propulsante/gasto                             #masa del misil, es necesario cambiarlo tambien abajo

'''
---------------------------CURVAS DE EMPUJE DEL F-4---------------------------
Previamente a este código, se han obtenido unas funciones de coeficientes
relativas a las curvas de empuje de la aeronave.  Partiendo del empuje máximo
del F-4 (el empuje máximo a nivel del mar) y con el dato de la altura, el cual
nos dará el valor de densidad, podremos obtener el valor del empuje equivalente
a cada altura en la que nos estemos moviendo.
'''

def thrust(mach, den):
    '''Cálculo del empuje de la aeronave.  Esta función se obtiene a partir de
    las gráficas del empuje del motor GE F404-400, cf. "Thrust Data for
    Performance Calculations" en M. Saarlas, "Aircraft Performance", p. 273.
    Se tiene en cuenta que la aeronave cuenta con dos motores GE F404-400.
    '''
    d_th = den / density(0)
    i = (.050618013228 + .11323534299 * d_th + 7.8263530571 * d_th**2
         - 15.012158645 * d_th**3)
    a_th = 1.1062543547 * d_th**1.276913816
    c_th = d_th * .862301392 + 1.937299323
    z_th = -.347382668*d_th + 1.71160358
    return TH_SL * (a_th + i * exp(-c_th * (mach - z_th)**2))

'''
--------------------------COEFICIENTES AERODINÁMICOS--------------------------
'''

def cl_alfa(mach):
    '''Cálculo de la pendiente del coeficiente de sustentación, que es lineal
    respecto del ángulo de ataque.  Este coeficiente varía con respecto al Mach
    y, con el ángulo de ataque, será posible obtener el coeficiente de
    sustentación.  Como el perfil es simétrico, el término independiente de la
    recta es nulo (CL0 = 0), por lo que el coeficiente de sustentación es
    proporcional al ángulo de ataque y cl_alfa(mach) da la constante de
    proporcionalidad en función del número de Mach.
    '''
    if mach <= .8:
        return (-3.4102 * mach**5 + 1.9918 * mach**4 + 2.9597 * mach**3
                - 1.6251 * mach**2 + .4172 * mach + 2.7915)
    if 0.8 < mach < 1.05:
        return 4.1216 * mach**3 - 13.250 * mach**2 + 14.343 * mach - 1.8055
    if mach >= 1.05:
        return .7534 * mach**3 - 4.2913 * mach**2 + 6.5935 * mach + .3476

def angulo_ataque(alfa_posible, mach):
    '''En función de si el ángulo de ataque obtenido es inferior o superior al
    de pérdida, la función angulo_ataque devolverá un ángulo u otro.  Cuando se
    supera el ángulo de ataque de entrada en pérdida, la función angulo_ataque
    devuelve el valor del ángulo de entrada en pérdida para así no volar en
    pérdida.
    '''
    if mach < .4:
        angulo_perdida = radians(16)
    else:
        angulo_perdida = radians(.8262 * mach**3 - .37724 * mach**2 - 6.4264
                                 * mach + 18.05)
    if alfa_posible < angulo_perdida:
        return alfa_posible
    return angulo_perdida

def cd0(mach):
    '''CD0 (coeficiente de resistencia parásita).  El estudio de este
    coeficiente se da en dos partes: una incompresible y otra compresible.  Es
    función del número de Mach y de distintos coeficientes referidos a las
    partes del avión.
    '''
    #Coeficientes de las distintas partes del avión, son coeficientes
    # experimentales que se introducen en la fórmula para calcular CD0.
    x_ala = .003
    x_fuselaje = .0024
    x_gondolas = .006
    x_med = (x_fuselaje + x_gondolas) / 2
    x_cola = .0025
    #Áreas totales de las partes del avión mojdadas por la corriente de aire y
    # que, por tanto, influyen en la resistencia.
    a_ala = 2 * S_W
    a_fuselaje = .75 * pi * LF * BF
    a_gondolas = pi * LM * DM * NE
    a_cola = 2 * (S_H + S_V)
    incremento = 1.1 # Incremento por interferencias e imperfecciones (10%)
    #No se tiene en cuenta el efecto de la compresibilidad; pero su cálculo es
    # necesario para el cálculo del coeficiente de resistencia parásita con
    # compresibilidad.
    cd0_inc = ((x_ala * a_ala + x_med * a_fuselaje + x_gondolas * a_gondolas
                + x_cola * a_cola) * incremento) / S_W
    n_comp = 3 / (1 + (1 / AR)) # Coeficiente de la fórmula con compresibilidad.
    #Cálculo de CD0 con efectos de compresibilidad.
    #Se emplea Md98, el 98% del Mach de divergencia ya que la fórmula tiende
    # asintóticamente a infinito en rangos cercanos al Mach de divergencia.
    if mach < M_D098 and mach < M_C:
        cd0_compresible = cd0_inc / ((1 - (mach / M_D)**2)**.25)
        cd_incremento = 0
        return cd0_compresible + cd_incremento
    #Para valores superiores al Mach crítico aparece un incremento por
    # resistencia de onda.
    if mach < M_D098 and mach >= M_C:
        cd0_compresible = cd0_inc / ((1 - (mach / M_D)**2)**.25)
        cd_incremento = (K / 1000) * ((10 * (mach - M_C)) / ((1 / cos(FLECHA))
                                                             - M_C))**n_comp
        return cd0_compresible + cd_incremento
    if 1.2 >= mach > M_D098:
        return (-379.32512053 * mach**5 + 1994.1499524 * mach**4
                - 4177.4704011 * mach**3 + 4358.3944768 * mach**2
                - 2264.4097020 * mach + 468.71342687)
    if mach > 1.2:
        return .031

def k(mach):
    '''Coeficiente de resistencia inducida que multiplica al coeficiente
    de sustentación.
    '''
    fos = .005 * (1 + 1.5 * (ESTRECHAMIENTO - .6)**2)
    #Este valor es una función lambda que aparece dentro del factor de Oswald.
    e_mach = 1 / ((1 + .12 * mach**2) * (1 + (.1 * (3 * NE + 1)) / (4 + AR)
                                         + (.142 + fos * AR * (10
                                                               * ESPESOR)**.33)
                                         / (cos(FLECHA2)**2)))  # Factor de Oswald
    return 1 / (e_mach * pi * AR)

def cd_inducida(k_d, c_l):
    '''Coeficiente de resistencia inducida.
    '''
    return k_d * c_l**2

def resistencia(vel, dens, c_d):
    '''Fuerza aerodinámica de resistencia total.
    '''
    return .5 * dens * S_W * c_d * vel**2

def sustentacion(vel, dens, c_l):
    '''Fuerza aerodinámica de sustentación total.
    '''
    return .5 * dens * S_W * c_l * vel**2


def Cdll(Ml):

    def coef_resistencia_base_misil(Ml):
    
    	if Ml < 0.8:
    		x0 = 0
    		x1 = 0
    		x2 = 0
    		x3 = 0
    		x4 = 0
    	elif Ml < 1:
    		x0 = -1.548523
    		x1 = 6.05972764
    		x2 = -7.30548391
    		x3 = 2.96129532
    		x4 = 0
    	elif Ml < 1.1:
    		x0 = 5.79090984*10**3
    		x1 = -2.19843314*10**4
    		x2 = 3.12774812*10**4
    		x3 = -1.97644892*10**4
    		x4 = 4.68059822*10**3
    	elif Ml < 1.5:
    		x0 = -4.11856506
    		x1 = 1.42267421*10**1
    		x2 = -1.69678524*10**1
    		x3 = 8.771665
    		x4 = -1.67398037
    	elif Ml < 2.2:
    		x0 = 3.0748*10**-1
    		x1 = -1.3258*10**-1
    		x2 = 2.8812*10**-2
    		x3 = 0
    		x4 = 0
    	elif Ml <=3.5:                            
    		x0 = 1.8481*10**-1
    		x1 = -2.2895*10**-2
    		x2 = 5.1876*10**-3
    		x3 = -4.0742*10**-4
    		x4 = 0
    	elif Ml >3.5:                            
    		x0 = 0.15
    		x1 = 0
    		x2 = 0
    		x3 = 0
    		x4 = 0
    		    		
    	return x4*Ml**4 + x3*Ml**3 + x2*Ml**2 + x1*Ml + x0
     
    CD_base_misil = coef_resistencia_base_misil(Ml)*Ratio_areas
    
  
    #CALCULO DEL COEFICIENTE DE FRICCION
    ##COEFICIENTE DE FRICCIÓN DEL CONO
    ###CALCULO DEL REYNOLDS1 
    
    Re_cono=rho*vl*longitud_cono/Mu_Visc #REYNOLDS 2
    
    def cfcono_misil(Re_cono):
    		
        #LAMINAR
    	if Re_cono < 1e6 :
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    	    cfi_cono=0.664*Re_cono**(-1/2)
    		
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL MEDIO
    
    	    cf_cono=2*cfi_cono
    
    			#CALCULO COEFICIENTE DE FRICCION COMPRESIBLE
    
    	    cfm_cono=cf_cono*(1/(1+0.17*Ml**2))**0.1295
    
    			#CALCULO COEFICIENTE DE FRICCION DEL CONO
    	    
    		
    		#TURBULENTO
    	else:	
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfi_cono=0.288*((log10(Re_cono))**(-2.45))
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL COMPRESIBLE
    
    		cf_cono=cfi_cono*1.597*((log10(Re_cono))**(-0.15))
    
    			#CALCULO COEFICIENTE DE FRICCION MEDIO
    
    		cfm_cono=cf_cono*(1/(1+(GAMMA-1)/2*Ml**2)**0.467)
    
    			#CALCULO COEFICIENTE DE FRICCION DEL CONO
    
    		
    		
    	return cfm_cono*Sup_cono/Sref_misil
    
    
    ##COEFICIENTE DE FRICCIÓN DEL CILINDRO
    
    
    Re_cil=rho*vl*(longitud_misil-longitud_cono)/Mu_Visc  #REYNOLDS 2
    
    
    def cfcil(Re_cil):
    	
    		
    		#LAMINAR
    	if Re_cil < 1e6 :
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfi_cil=0.664*Re_cil**(-1/2)
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL MEDIO
    
    		cf_cil=2*cfi_cil
    
    			#CALCULO COEFICIENTE DE FRICCION COMPRESIBLE
    
    		cfm_cil=cf_cil*(1/(1+0.17*Ml**2))**0.1295
    
    			
    		
    		
    		
                #TURBULENTO
    	else:			
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfi_cil=0.288*((log10(Re_cil))**(-2.45))
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL COMPRESIBLE
    
    		cf_cil=cfi_cil*1.597*((log10(Re_cil))**(-0.15))
    
    			#CALCULO COEFICIENTE DE FRICCION MEDIO
    
    		cfm_cil=cf_cil*(1/(1+(GAMMA-1)/2*Ml**2)**0.467)
    
    			
    
    		
    		
    	return cfm_cil*Sup_total/Sref_misil
    
    
    #CALCULO DEL COEFICIENTE DE FRICCION TOTAL REFERIDO A LA SUPERFICIE TRANSVERSAL
    CDFriccion_cono = cfcono_misil(Re_cono)
    CDFriccion_cil = cfcil(Re_cil)
    CDFriccion=CDFriccion_cono+CDFriccion_cil
    
    #CALCULO DEL COEFICIENTE DE ONDA
    def cd_onda(Ml,angulo_cono):
    	if Ml >= 1 :
    		
    		return (0.083+0.096/(Ml**2))*(angulo_cono/10)**1.69
       #REGIMEN SUBSONICO
    
    	elif Ml <1 :
    		
    		return 	((60/((longitud_cono/diametro_m)**3))+0.0025*(longitud_cono/diametro_m))*CDFriccion
    
    CD_onda = cd_onda(Ml,angulo_cono)
    	
    	
    ####RESISTENCIA ALETAS
    ######################
    
    #COEFICIENTE DE ONDA
    def cd_onda_aletas(Ml,angulo_cono):
    	if Ml >= 1 :
    		
    		return ((4*tao_aleta**2)/((Ml**2-1)**0.5))*Swtotal_aletas/Sref_misil
       #REGIMEN SUBSONICO
    
    	elif Ml <1 :
    		
    		return 	0
    
    CD_onda_aletas = cd_onda_aletas(Ml,angulo_cono)
        
    Re_aletas=rho*vl*Craiz_aleta/Mu_Visc #REYNOLDS aletas
    
    #coeficiente fricción aletas
    def cf_aletas(Re_aletas):
    		
        #LAMINAR
    	if Re_aletas < 1e6 :
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    	    cfialetas=0.664*Re_aletas**(-1/2)
    		
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL MEDIO
    
    	    cf1aletas=2*cfialetas
    
    			#CALCULO COEFICIENTE DE FRICCION COMPRESIBLE
    
    	    cfmaletas=cf1aletas*(1/(1+0.17*Ml**2))**0.1295
    
    			
    	    
    		
    		#TURBULENTO
    	else:	
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfialetas=0.288*(log10(Re_aletas))**(-2.45)
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL COMPRESIBLE
    
    		cf1aletas=cfialetas*1.597*((log10(Re_aletas))**(-0.15))
    
    			#CALCULO COEFICIENTE DE FRICCION MEDIO
    
    		cfmaletas=cf1aletas*(1/(1+(GAMMA-1)/2*Ml**2)**0.467)
    
    			
    
    		
    		
    	return cfmaletas*Swtotal_aletas/Sref_misil
    
    
    CDFriccion_aletas = cf_aletas(Re_aletas)
    
    return CD_base_misil + CDFriccion + CD_onda + CD_onda_aletas + CDFriccion_aletas


#A partir de aquí, se abre un bucle en función del ángulo theta ("beta" en el
# código).  Este ángulo sirve para determinar el ángulo final de la maniobra de
# giro; es decir, es el ángulo con el que se quiera que comience el ascenso
# tras el giro vertical.  Dado que se contempla un abanico de posibilidadaes,
# se decide estudiar in rango desde 10 grados hasta 90 grados.  Este bucle debe
# durar todo el programa para conseguir que exporte los distintos ficheros
# correspondientes a cada ángulo.

############-----Parametro a variar--------##########
f = open('2 etapas', 'w')  # Fichero de escritura sin extensión.
f.write('Gasto (kg/s) \tTiempo \tALTURA lanzamiento (m)\tVELOCIDAD (m/s)')
f.write('\tTHETA (deg)')
f.write('\tMasa  \tVelocidad perdida \n')




for gasto in range(1,5,2): 
    '''beta es el ángulo de asiento al final del giro y lo variaremos en el
    bucle para obtener los resultados en función de sus distintos valores.  Se
    creará un fichero para cada ángulo que contenga los tiempos, las alturas,
    las posiciones, las velocidades, los números de Mach, los ángulos de
    ataque, los ángulos de asiento de la velocidad, los ángulos de asiento, las
    fuerzas de sustentación, las energías mecánicas, las energías cinéticas,
    las energías potenciales, los empujes, las fuerzas de resistencia y los
    factores de carga correspondientes a cada instante de la maniobra.
    '''
    beta = 89
    beta = radians(beta)
    gasto_texto=str(gasto)       
    '''Inicialización de variables y diferenciales para la maniobra del misil'''
    
    vl = 467.22 #Inicialización de la velocidad
    yl = 13518.77  #Inicialización de la altitud 
    thetal = 30.95*pi/180 #Inicialización del ángulo de asiento
    thetalgrados=thetal*(180/pi) #Conversión de radianes a grados del ángulo de asiento


    
    vxl=vl*cos(thetal) #Inicialización de la componente horizontal de velocidad  
    vyl=vl*sin(thetal) #Inicialización de la componente verical de velocidad 
    tl = 0 #Inicialización temporal 
    xl = 0 #Inicialización de la posición en el eje x 
    sl = 0 #Inicialización del arco recorrido
    dvxl=0 #Inicialización del diferencial de la componente horizontal de velocidad 
    dvyl=0 #Inicialización del diferencial de la componente vertical de velocidad 
    dsl=0 #Inicialización del diferncial del arco recorrido
    dxl=0 #Inicialización del diferencial de la posición
    dyl=0 #Inicialización del diferencial de la altitud
    dtl=0.01 #Inicialización del diferencial de tiempo
    dthetal = 0 #Inicialización del diferencial del ángulo de asiento                                
    thetalgrados=thetal*(180/pi)
    
    
    rho = density(yl)  # Densidad inicial del aire (kg/m3).
    p = pressure(yl)  # Presión inicial del aire (Pa).
    T = temperature(yl)  # Temperatura inicial del aire (K).
    Mu_Visc = viscosity(yl) # Viscosidad 
    #A la altura inicial el avión vuela en vuelo estacionario.   
    dt = 0.1  # Diferencial de tiempo (s).
    
    Ml=vl/((GAMMA*R_AIR*T)**0.5)
    Cdl=Cdll(Ml)
    D_misil=0.5*rho*Cdl*Sref_misil*vl**2
    
                          #ETAPAS
                          #Menos la masa misil todos lo parametros pueden variar
    masa_misil=782.83 #ver linea 886
    masa_propulsante1=599.60
    masa_propulsante2=117.37
    masa_estructura1=29.98 #♣Hay que inventarlo
    masa_estructura2=5.86 #Hay que introducirlo
    Isp = 280*9.81
    Empuje_misil= gasto*Isp
    gasto2= gasto #♥Poner el que toque
    t_combustion= masa_propulsante1/gasto
    t_combustion_2_etapa=masa_propulsante2/gasto2 #Hay que calcularlo segun lo que surja
    retardo_encendido = 0 #de momento dejarlo asi
    t_fin_combustion2= t_combustion + t_combustion_2_etapa + retardo_encendido
     
    
    Velocidad_perdida=0

        
        
        
    '''BUCLE TIRO BALÍSTICO DEL MISIL'''
        
        
    while tl<= t_combustion and thetal>0:
            
        '''Con este bucle se calculan todas las variables correspondientes al tiro balístico del misil SIN EMPUJE . 
        La condición de parada del bucle es que el ángulo de asiento del misil deje de ser positivo. Esto es cuando
        el misil se encuentre completamente en posición horizontal. La particularidad de este caso
        es que el ángulo de lanzamiento del misil (ángulo de asiento) es el ángulo de asiento de la VELOCIDAD del avión.
        '''

        
        tl=tl+dtl #Evolución temporal (s)
        xl=xl+dxl #Posición horizontal (m)
        yl=yl+dyl #Altitud (m) 
        sl=sl+dsl #Arco recorrido (m)
        
        thetal=thetal+dthetal #Ángulo de asiento
        thetalgrados=thetal*(180/pi) #Conversión de radianes a grados
        
        
        vxl=vxl+dvxl #Componente horizontal de la velocidad (m/s)
        vyl=vyl+dvyl #Componente vertical de la velocidad (m/s)
        vl=(vxl**2+vyl**2)**0.5 #Módulo de la velocidad (m/s)
        
        r = RT + yl  # Distancia al centro de la Tierra (m).
        g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
        
        rho = density(yl)  # Densidad (kg/m3).
        T = temperature(yl)  # Temperatura (K).
        Mu_Visc = viscosity(yl) # Viscosidad 
        
        Ml=vl/((GAMMA*R_AIR*T)**0.5) #Mach de vuelo
        
        
        Ratio_areas=(Sref_misil-Sgases)/Sref_misil #Relación de áreas 
                        
              
        Cdl=Cdll(Ml)                     
        D_misil=0.5*rho*Cdl*Sref_misil*vl**2 #Fuerza de resistencia (N)
        Dx=D_misil*cos(thetal) #Componente horizontal de la fuerza de resistencia (N)
        Dy=D_misil*sin(thetal) #Componente vertical de la fuerza de resistencia (N)
        
        dvxl_perdidas=-(Dx/masa_misil)*dtl #Diferencial de la componente horizontal de la velocidad (m/s)
        dvyl_perdidas=-g0*dtl-(Dy/masa_misil)*dtl #Diferencial de la componente vertical de la velocidad (m/s)

        
        dvxl=dvxl_perdidas+Empuje_misil*cos(thetal)*dtl/masa_misil
        dvyl=dvyl_perdidas+Empuje_misil*sin(thetal)*dtl/masa_misil
        masa_misil=masa_misil-gasto*dtl
            
#            #COMBUSITON SEGUNDA ETAPA    
#            if t_fin_combustion2>=tl>t_combustion:
#                
#                masa_misil = masa_misil - masa_estructura1
#                dvxl=dvxl_perdidas+Empuje_misil*cos(thetal)*dtl/masa_misil
#                dvyl=dvyl_perdidas+Empuje_misil*sin(thetal)*dtl/masa_misil
#                masa_misil=masa_misil-gasto2*dtl

        dthetal=dtl*g0*cos(thetal)/(-vl) #Diferencial del ángulo de asiento
        dxl=vxl*dtl #Diferencial de la posición (m)
        dyl=vyl*dtl #Diferencial de la altitud (m)
        dsl=vl*dtl #Diferencial del arco recorrido (m)
        
        Velocidad_perdida = Velocidad_perdida + (dvxl_perdidas**2 + dvyl_perdidas**2)**0.5
    
    masa_misil = masa_misil - masa_estructura1
    print(gasto)
    print('1etapa velocidad',vl)
    print('1etapa altura',yl)
    print('1etapa',masa_misil)
    print('1etapa',thetalgrados)
    while t_fin_combustion2>=tl>t_combustion and thetal>0 and yl<500000:
    
        '''Con este bucle se calculan todas las variables correspondientes al tiro balístico del misil SIN EMPUJE . 
        La condición de parada del bucle es que el ángulo de asiento del misil deje de ser positivo. Esto es cuando
        el misil se encuentre completamente en posición horizontal. La particularidad de este caso
        es que el ángulo de lanzamiento del misil (ángulo de asiento) es el ángulo de asiento de la VELOCIDAD del avión.
        '''

        
        tl=tl+dtl #Evolución temporal (s)
        xl=xl+dxl #Posición horizontal (m)
        yl=yl+dyl #Altitud (m) 
        sl=sl+dsl #Arco recorrido (m)
        
        thetal=thetal+dthetal #Ángulo de asiento
        thetalgrados=thetal*(180/pi) #Conversión de radianes a grados
        
        
        vxl=vxl+dvxl #Componente horizontal de la velocidad (m/s)
        vyl=vyl+dvyl #Componente vertical de la velocidad (m/s)
        vl=(vxl**2+vyl**2)**0.5 #Módulo de la velocidad (m/s)
        
        r = RT + yl  # Distancia al centro de la Tierra (m).
        g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
        
        rho = density(yl)  # Densidad (kg/m3).
        T = temperature(yl)  # Temperatura (K).
        Mu_Visc = viscosity(yl) # Viscosidad 
        
        Ml=vl/((GAMMA*R_AIR*T)**0.5) #Mach de vuelo
        
        
        Ratio_areas=(Sref_misil-Sgases)/Sref_misil #Relación de áreas 
                        
              
        Cdl=Cdll(Ml)                     
        D_misil=0.5*rho*Cdl*Sref_misil*vl**2 #Fuerza de resistencia (N)
        Dx=D_misil*cos(thetal) #Componente horizontal de la fuerza de resistencia (N)
        Dy=D_misil*sin(thetal) #Componente vertical de la fuerza de resistencia (N)
        
        dvxl_perdidas=-(Dx/masa_misil)*dtl #Diferencial de la componente horizontal de la velocidad (m/s)
        dvyl_perdidas=-g0*dtl-(Dy/masa_misil)*dtl #Diferencial de la componente vertical de la velocidad (m/s)
            
        dvxl=dvxl_perdidas+Empuje_misil*cos(thetal)*dtl/masa_misil
        dvyl=dvyl_perdidas+Empuje_misil*sin(thetal)*dtl/masa_misil
        masa_misil=masa_misil-gasto2*dtl

        dthetal=dtl*g0*cos(thetal)/(-vl) #Diferencial del ángulo de asiento
        dxl=vxl*dtl #Diferencial de la posición (m)
        dyl=vyl*dtl #Diferencial de la altitud (m)
        dsl=vl*dtl #Diferencial del arco recorrido (m)
        
        Velocidad_perdida = Velocidad_perdida + (dvxl_perdidas**2 + dvyl_perdidas**2)**0.5
    masa_misil = masa_misil - masa_estructura2
    print(gasto)
    print('12tapa velocidad',vl)
    print('2etapa altura',yl)
    print('2etapa',masa_misil)
    print('2etapa',thetalgrados)
    while thetal>0 and yl<500000:
        '''Con este bucle se calculan todas las variables correspondientes al tiro balístico del misil SIN EMPUJE . 
        La condición de parada del bucle es que el ángulo de asiento del misil deje de ser positivo. Esto es cuando
        el misil se encuentre completamente en posición horizontal. La particularidad de este caso
        es que el ángulo de lanzamiento del misil (ángulo de asiento) es el ángulo de asiento de la VELOCIDAD del avión.
        '''

        
        tl=tl+dtl #Evolución temporal (s)
        xl=xl+dxl #Posición horizontal (m)
        yl=yl+dyl #Altitud (m) 
        sl=sl+dsl #Arco recorrido (m)
        
        thetal=thetal+dthetal #Ángulo de asiento
        thetalgrados=thetal*(180/pi) #Conversión de radianes a grados
        
        
        vxl=vxl+dvxl #Componente horizontal de la velocidad (m/s)
        vyl=vyl+dvyl #Componente vertical de la velocidad (m/s)
        vl=(vxl**2+vyl**2)**0.5 #Módulo de la velocidad (m/s)
        
        r = RT + yl  # Distancia al centro de la Tierra (m).
        g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
        
        rho = density(yl)  # Densidad (kg/m3).
        T = temperature(yl)  # Temperatura (K).
        Mu_Visc = viscosity(yl) # Viscosidad 
        
        Ml=vl/((GAMMA*R_AIR*T)**0.5) #Mach de vuelo
        
        
        Ratio_areas=1 #Relación de áreas 
                        
              
        Cdl=Cdll(Ml)                     
        D_misil=0.5*rho*Cdl*Sref_misil*vl**2 #Fuerza de resistencia (N)
        Dx=D_misil*cos(thetal) #Componente horizontal de la fuerza de resistencia (N)
        Dy=D_misil*sin(thetal) #Componente vertical de la fuerza de resistencia (N)
        
        dvxl_perdidas=-(Dx/masa_misil)*dtl #Diferencial de la componente horizontal de la velocidad (m/s)
        dvyl_perdidas=-g0*dtl-(Dy/masa_misil)*dtl #Diferencial de la componente vertical de la velocidad (m/s)
        
        dvxl=dvxl_perdidas
        dvyl=dvyl_perdidas
        
        dthetal=dtl*g0*cos(thetal)/(-vl) #Diferencial del ángulo de asiento
        dxl=vxl*dtl #Diferencial de la posición (m)
        dyl=vyl*dtl #Diferencial de la altitud (m)
        dsl=vl*dtl #Diferencial del arco recorrido (m)
        
        Velocidad_perdida = Velocidad_perdida + (dvxl_perdidas**2 + dvyl_perdidas**2)**0.5 
        
    print(gasto)
    print('2tapa velocidad',vl)
    print('2tapa altura',yl)
    print('2etapa',masa_misil)
    print('2etapa',thetalgrados)
    f.write('%.8f\t' %gasto) #Velocidad (m/s)
    f.write('%.8f\t' %tl) #Tiempo(s))        
    f.write('%.8f\t' %yl) #Altitud (m)
    f.write('%.8f\t' %vl) #Velocidad (m/s)
    f.write('%.8f\t' %thetalgrados) #Ángulo de asiento
    f.write('%.8f\t' %masa_misil) #Ángulo de asiento
    f.write('%.8f\n' %Velocidad_perdida) #Energía mecánica final (J) 

        
        



        
    
f.close()
           
'''
Como resumen:
1) El código ha empezado en una condición de vuelo uniforme.
2) La siguiente maniobra es un giro ascendente, a factor de carga máximo y
constante, y con mínima resistencia.
3) La última maniobra es un ascenso con el ángulo final del giro, con
coeficiente de sustentación óptimo.

Este programa nos exportará un total de 9 archivos (uno para cada incremento de
10º del ángulo final de giro) que, exportados a Excel nos permiten observar
cómo cambian las variables según las condiciones de vuelo.
'''
    

            
        
             