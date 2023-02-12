#Etapes 1: Importation des modules

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

# Etapes 2: Donner d'entrées (variables)

L=1.5                   # Longueur des 4 cotés de la plaque en m
tf= float(2*3600)       # Temps d'évolution de 2h
alpha=2.10E-5           # Diffusivité thermique en m2/s
T1= 400.                # Température des cotés horizontaux  en °C
T2=300.                 # Temperature des cotés vertiticaux  en °C
T0= 293                 # Température initial sur les 4 cotés en °C
#lamb= 0                 # Conductivité thermique du matériaux en W.m^-1.K^-1
#Rho=0                   # Masse volumique du matériau en kg.^-3
#c= 0                    # Capacité thermique massique du materiau en J.kg^-1.K^-1

Nx= 100               # Nombres de points de calcule sur la barre
Ny= 100              # Nombres de points de calcule sur la barre
dt=1.                 # Fréquence de calcule
Nt= int(tf/dt)        # Nombre de point de calcule dans le temps
dx=float(L/100)       # Pas d'espace
dy=float(L/100)
K=alpha*(dt/(dx*dx*dy*dy))  # Constante K
x=np.linspace(0,L,num=Nx+1)
y=np.linspace(L,0,num=Ny+1)

# Etape 3: Vérification de la condition de stabilité CFL
# K > 0.5 pour être stable

CFL=K
print("CFL = ",CFL)
if CFL<=0.5:
    print("La condition CFL est respectée")
else:
    print("La condition CFL n'est pas respectée")

#Etapes 4: Condition initial

Tx=sp.ones((Nt+1,Nx+1))        # Matrice de température T(x,t)
Tx=293*Tx                       #  Ou 293*np.ones((Nx+1,Nt+1))
for i in range(0,Nx+1):
    Tx[i,1]= 293.

Ty=sp.ones((Ny+1,Nt+1))       # Matrice de température T(y,t)
Ty=293*Ty
for j in range(0,Ny+1):
    Ty[j,1]= 293.

# Etapes 5: Les conditions aux limites

for n in range(0,Nt+1):
    Tx[n,0]=T2
    Tx[n,Nx]=T2

for n in range(0, Nt+1):
    Ty[0,n]=T1
    Ty[Ny,n]=T1

# Profile de température initial dans le barreau
plt.figure(1)
plt.plot(x,Tx[:,1])
plt.plot(y, Ty[:,1])
plt.legend(["Tx pour t=0s","Ty pour t=0s"])

# Etape 6: Calcul de la température à tous les Nt ponts et pour tous les noeuds
# du barreau

for n in range(1,Nt):
    for i in range(1,Nx):
        for j in range(1,Ny):
            Tx[i,n+1]=K*((Tx[i+1,n]+Tx[i-1,n])*dy*dy+(Ty[j+1,n]+Ty[j-1,n])*dx*dx-2*(Tx[i,n]*dx*dx+Ty[j,n]dy*dy)+Tx[i,n])

# Etape 7: Graphiques

#Graphique pour tracer l'évolution de la température au début
#Aumilieur et à la fin du processus

plt.figure(2)
plt.plot(x,Tx[:,1],'-b',x,Tx[:,(int(Nt/2))],'--k',x,Tx[:,Nt],'-r')
# ATTENTION que toutes les valeurs soit bien interpreté comme un nombre entier
plt.legend(['T pour t=0s','Tpour t=1h','T pour t=2h'])
plt.xlabel('Longueur x (m)')
plt.ylabel('Température T(K)')