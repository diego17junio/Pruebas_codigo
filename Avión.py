# -*- coding: utf-8 -*-
"""
@author: Team REOS
"""
from math import exp, radians, cos, pi, sin, degrees, atan, log10, log

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
CAM = 4.33
ALT_avion = 2.11
ALT_cuerpo = 1.46
Sbase = BF*ALT_cuerpo 

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
    if M < .4:
        angulo_perdida = radians(16)
    else:
        angulo_perdida = radians(.8262 * mach**3 - .37724 * mach**2 - 6.4264
                                 * mach + 18.05)
    if alfa_posible < angulo_perdida:
        return alfa_posible
    return angulo_perdida

def cd0(mach):
    
    a_ala = 2 * S_W
    a_fuselaje = .75 * pi * LF * BF
    a_gondolas = pi * LM * DM * NE
    a_cola = 2 * (S_H + S_V)
    FRcuerpo = LF/(BF*ALT_avion)
        
    FFcuerpo = 1+ 60/FRcuerpo**3+0.0025*FRcuerpo
        
    FFalas = 1.18
        
    FI = 1.05
    Rec = (rho*v*LF)/Mu_Visc
    Realas = (rho*v*CAM)/Mu_Visc
    
    
    if mach<0.65:
        mach = 0.65
        v065 = mach * (GAMMA * R_AIR * T)**.5
        Rec = (rho*v065*LF)/Mu_Visc
        Realas = (rho*v065*CAM)/Mu_Visc 
        cfriccioncuerpo = a_fuselaje/S_W*FFcuerpo*FI*0.42/((log(0.056*Rec))**2)
        cfriccionalas = a_ala/S_W*FFalas*FI*0.42/((log(0.056*Realas))**2)
        cdbase = (0.1+0.1222*mach**8)*Sbase/S_W
        
        return cfriccioncuerpo + cfriccionalas + cdbase
    
    
    if 0.65<=mach<0.8:
        cfriccioncuerpo = a_fuselaje/S_W*FFcuerpo*FI*0.42/((log(0.056*Rec))**2)
        cfriccionalas = a_ala/S_W*FFalas*FI*0.42/((log(0.056*Realas))**2)
        cdbase = (0.1+0.1222*mach**8)*Sbase/S_W
        
        return cfriccioncuerpo + cfriccionalas + cdbase
    
    if 0.8<mach<1:
        cfriccioncuerpo = a_fuselaje/S_W*FFcuerpo*FI*0.42/((log(0.056*Rec))**2)
        cfriccionalas = a_ala/S_W*FFalas*FI*0.42/((log(0.056*Realas))**2)
        cdbase = (0.1+0.1222*mach**8)*Sbase/S_W
        
        return cfriccioncuerpo + cfriccionalas + cdbase 
    
    if 1<mach<1.2:
        cfriccioncuerpo = a_fuselaje/S_W*0.42/((log(0.056*Rec))**2)
        cfriccionalas = a_ala/S_W*0.42/((log(0.056*Realas))**2)
        cdbase = 0.2222*Sbase/S_W
        cdwavefus = 9/8*pi**2*(1/(5.88**2))*Sbase/S_W
        cdwavealas = 16/3/(mach**-1)**(-1/2)*(0.0601)**2
        
        return cfriccioncuerpo + cfriccionalas + cdbase + cdwavefus + cdwavealas
    
    if 1.2<mach<1.8:
        cfriccioncuerpo = a_fuselaje/S_W*0.42/((log(0.056*Rec))**2)
        cfriccionalas = a_ala/S_W*0.42/((log(0.056*Realas))**2)
        cdbase = 0.2222*Sbase/S_W
        cdwavefus = 9/8*pi**2*(1/(5.88**2))*Sbase/S_W
        cdwavealas = 16/3/(mach**-1)**(-1/2)*(0.0601)**2
        
        return cfriccioncuerpo + cfriccionalas + cdbase + cdwavefus + cdwavealas
        
    if 1.8<=mach:
        cfriccioncuerpo = a_fuselaje/S_W*0.42/((log(0.056*Rec))**2)
        cfriccionalas = a_ala/S_W*0.42/((log(0.056*Realas))**2)
        cdbase = (1.42*Sbase/S_W)/(3.15+mach**2)
        cdwavefus = 9/8*pi**2*(1/(5.88**2))*Sbase/S_W
        cdwavealas = 16/3/(mach**-1)**(-1/2)*(0.0601)**2
        
        return cfriccioncuerpo + cfriccionalas + cdbase + cdwavefus + cdwavealas

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



beta = 89
beta = radians(beta)
f = open('resultados_avion', 'w')  # Fichero de escritura sin extensión.
f.write('ALTURA lanzamiento (m)\tVELOCIDAD (m/s)')
f.write('\tTHETA (deg)'  '\n')
#Cabezas de tabla.
'''-------------------------CONDICIONES INICIALES-------------------------
Ahora, para los próximos cálculos, se definen las variables termodinámicas
obtenidas del modelo ISA.
'''
h = 12000  # Altitud inicial (m).
r = RT + h  # Distancia desde el centro de la Tierra (m).
g0 = MU / r**2  # Aceleración g0itatoria (m/s2).
W = 14273 * g0
rho = density(h)  # Densidad inicial del aire (kg/m3).
p = pressure(h)  # Presión inicial del aire (Pa).
T = temperature(h)  # Temperatura inicial del aire (K).
Mu_Visc = viscosity(h) # Viscosidad 
#A la altura inicial el avión vuela en vuelo estacionario.
M = 1.8  # Número de Mach inicial.
v = M * (GAMMA * R_AIR * T)**.5  # Velocidad inicial (m/s).
CL_alfa1 = cl_alfa(M)  # Pendiente del coeficiente de sustentación.
#Ángulos de asiento, de ataque y de asiento de la velocidad iniciales.
alfa_numerico = 2 * W / (rho * v**2 * S_W * CL_alfa1)
alfa = angulo_ataque(alfa_numerico, M)  # Ángulo de ataque.
alfa_grados = degrees(alfa)  # Ángulo de ataque en grados.
gama = 0  # Ángulo de asiento.
gama_grados = degrees(gama)  # Ángulo de asiento de la velocidad en grados.
theta = gama + alfa  # Ángulo de asiento de la velocidad.
theta_grados = degrees(theta)  # Ángulo de asiento en grados.
#Coeficientes aerodinámicos
CL = alfa*CL_alfa1  # Coeficiente de sustentación inicial.
k1 = k(M)
CD01 = cd0(M)
CD_inducida1 = cd_inducida(k1, CL)
CD = CD01 + CD_inducida1  # Polar del avión.  Coeficiente de resistencia.
'''-------------------------INICIO DE LA MANIOBRA-------------------------
'''
#Se impone un valor constante del radio de giro, es decir, giro ascendente
# a radio constante.  En futuras versiones de este cálculo, esto se
# cambiará para buscar una solución más realista.  Por ahora, con objeto de
# facilitar cálculos, se dejará así.
radius = v**2 / (g0 * (N - 1))  #Radio de giro inicial (m).
#Este radio de giro se obtiene para la velocidad inicial en vuelo
# estacionario y para un factor de carga máximo según los pilones de carga
# n = 3,5.
dt = 0.1  # Diferencial de tiempo (s).
'''------SISTEMA DE ECUACIONES PARA PRIMER TRAMO: VUELO ESTACIONARIO------
'''
#Energías.
ecinetica = .5 * (W/g0) * v**2  # Energía cinética (J).
epotencial = (W/g0) * g0 * h  # Energía potencial (J).
emecanica = ecinetica + epotencial  # Energía mecánica (J).
#Fuerzas.
D = resistencia(v, rho, CD)  # Resistencia aerodinámica (N).
L = sustentacion(v, rho, CL)  # Sustentación aerodinámica (N).
Th = thrust(M, rho)  # Empuje (N).
diferencia_T_D = Th - D
#Esto nos va a permitir calcular en qué momento el empuje se verá superado
# por la resistencia
n = L / W  # Factor de carga.
#Condiciones iniciales para la integración.
t = 0
x = 0
omega = v / radius  # Velocidad angular en la maniobra de giro (rad/s).
#Segunda ley de Newton en el eje horizontal (ejes cuerpo).
dv = dt * (Th * cos(alfa) - D - (W/g0) * g0 * sin(gama)) / (W/g0)
vx = v * cos(gama)  # Velocidad horizontal del avión en ejes tierra (m/s).
vy = v * sin(gama)  # Velocidad vertical del avión en ejes tierra (m/s).
dx = vx * dt  # Variación horizontal de la posición en ejes tierra (m).
dh = vy * dt  # Variación vertical de la posición en ejes tierra (m).
dtheta = v * dt / radius  # Variación del ángulo de empuje.



'''-------SISTEMA DE ECUACIONES PARA SEGUNDO TRAMO: MANIOBRA DE GIRO-------
Ahora comienza el bucle relativo al giro ascendente, que analiza la
trayectoria con nuevas ecuaciones y condiciones de vuelo que se detallan
más afelante.  El significado de theta < beta implica que el bucle realice
el cálculo requerido siempre que el ángulo theta sea menor que beta.  Se ha
obligado a que beta sea el ángulo de final de giro (al inicio del programa
se le ha dado un rango de valores).  Por tanto, una vez que theta sea igual
a beta, se dará por concluida la maniobra de giro y se comenzará con el
siguiente tramo. Todo ello mientras la velocidad sea positiva, condición
que permite ver si el resultado es realista o no.
'''
while gama < beta and v > 0:
    '''A continuación, se escriben en un fichero todas las variables aquí
    detalladas para cada valor de theta < beta.
    '''
    
    f.write('%.8f\t' %h)  # Altitud (m).
    f.write('%.8f\t' %v)  # Velocidad (m/s).   
    f.write('%.8f\n' %gama_grados)  # Asiento de la velocidad (grados).
    
    #Ya que este análisis de maniobra, a diferencia del anterior, lleva un
    # cálculo para distintos valores de tiempo y velocidad, se debe
    # programar su evolución en términos de sus variaciones diferenciales.
    #Aquí se escriben los incrementos diferenciales de las coordenadas
    # espaciales, temporales y de velocidad.  Los diferenciales se obtienen
    # del tramo anterior, y sirven para que estos nuevos valores actúen en
    # las nuevas condiciones para calcular nuevas variables.
    
    '''Inicialización de variables y diferenciales para la maniobra del misil'''
    
    thetal = gama #Inicialización del ángulo de asiento
    thetalgrados=thetal*(180/pi) #Conversión de radianes a grados del ángulo de asiento
    yl = h #Inicialización de la altitud 
    vl = v #Inicialización de la velocidad
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
    thetal = gama #Ángulo de asiento del avión
                  #es igual al ángulo de asiento de la velocidad del misil
                                
    thetalgrados=thetal*(180/pi)
    Ddsl=0
    
    
    
    t = t + dt  # Evolución temporal (s).
    v = v + dv  # Velocidad (m/s).
    x = x + dx  # Posición horizontal (m).
    h = h + dh  # Altitud (m)
    r = RT + h  # Distancia al centro de la Tierra (m).
    g0 = MU / r**2  # Aceleración de la g0edad (m/s2).
    W = 14273 * g0
    #Las variables termodinámicas habrán variado con la altura.
    rho = density(h)  # Densidad (kg/m3).
    T = temperature(h)  # Temperatura (K).
    Mu_Visc = viscosity(h) # Viscosidad 
    M = v / (GAMMA * R_AIR * T)**.5 # Mach de vuelo.
    n = 3.5  # Tomamos la condición de factor de carga máximo y constante.
    radius = v**2 / (g0 * (n - 1))  # Radio de giro varía con la velocidad.
    #Las características aerodinámicas varían con el número de Mach.
    CL_alfa1 = cl_alfa(M)
    k1 = k(M)
    CD01 = cd0(M)
    alfa_numerico = 2 * W * n / (rho * v**2 * S_W * CL_alfa1)
    #El nuevo ángulo de ataque resultará del nuevo factor de carga (antes 1
    # y ahora máximo), la nueva velocidad y las nuevas características
    # aerodinámicas.
    alfa = angulo_ataque(alfa_numerico, M)  # Ángulo de ataque.
    theta = theta + dtheta  # Ángulo de asiento (empuje horizontal).
    gama = theta - alfa  # Ángulo de asiento de la velocidad.
    # Ángulos en grados para la exportación a los ficheros.
    alfa_grados = degrees(alfa)  # Ángulo de ataque (grados).
    theta_grados = degrees(theta)  # Ángulo de asiento (grados).
    gama_grados = degrees(gama)  # Asiento de la velocidad (grados).
    CL = alfa * CL_alfa1  # Coeficiente de sustentación.
    n = .5 * rho * v**2 * S_W * CL / W
    #Nuevas proyecciones de la velocidad (ejes tierra).
    vx = v * cos(gama)  # Proyección horizontal de la velocidad (m/s).
    vy = v * sin(gama)  # Proyección vertical de la velocidad (m/s).
    CD_inducida1 = cd_inducida(k1, CL)
    CD = CD01 + CD_inducida1  # Polar del avión
    D = resistencia(v, rho, CD)  # Fuerza de resistencia (N).
    L = n * W  # Fuerza de sustentación (N).
    #Energías.
    ecinetica = .5 * (W/g0) * v**2  # Energía cinética (J).
    epotencial = (W/g0) * g0 * h  # Energía potencial (J).
    emecanica = ecinetica + epotencial  # Energía mecánica (J).
    #Empuje.
    Th = thrust(M, rho)  # Empuje (N).
    diferencia_T_D = Th - D
    #Velocidad angular.
    omega = v / radius  # Velocidad angular (rad/s).
    #Nuevas variaciones diferenciales.  Segunda ley de Newton.
    dv = dt * (Th * cos(alfa) - D - W * sin(gama)) / (W/g0)
    #Esta última ecuación nos permite obtener un nuevo incvremento de
    # velocidad (m/s) que, al reiniciar el bucle será sumado al último
    # valor de la velocidad, obteniendo la variación de la velocidad.
    dx = vx * dt  # Variación horizontal de la posición (m).
    dh = v * sin(gama) * dt  # Variación vertical de la posición (m).
    dtheta = omega * dt  # Variación del ángulo de asiento.


    

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
    

            
        
             