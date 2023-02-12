import numpy as np
M=100
N=100
T=1
#mesh in space
x=np.linspace(0,1,N+1)
#mesh in time
t=np.linspace(0,T,M+1)
#les coefficients
D=0.0000001 #diffusion coefficient (le schéma est stable lorsque D*N²<<s)
s=6 #intesity

p=D*(T/M)*(N**2)
q=s*(T/M)

u=np.zeros((M+1,N+1)) # matrice pour stocker les données de u(x,t),
# la i'eme ligne de la matrice correspond à u en fonction de x à l'instant ti

#Rq: le premier coefficient de la fonction u est x et le deuxiéme est t ,
# mais dans notre code , le premier coefficient de la matrice u correspond à t et le deuxiéme correspond à x.

B=np.zeros(N+1) # les coefficients B(i) de la diagonale de la matrice A
A=np.zeros((M+1,N+1)) # la matrice tridiagonale A pour la résolution de systéme linéaire

#résolution du système d'équations
from numpy.linalg import solve

# u(x,0) = 1/(1 + exp(x))² ; 0 ≤ x ≤ 1 (les conditions au bord pour le temps)
for i in range(N + 1):
    u[0][i] = 1 / ((1 + np.exp(i / N)) ** 2)

# les coefficients de la matrice A (à linstant 0)
for i in range(1,N):
    B[i]=1+2*p-q*(1-u[0][i]) #les coefficients de la diagonale de A à linstant 0
    A[i,i-1]=-p
    A[i,i+1]=-p
    A[i,i]=B[i]
#les coefficients de A qui ne sont pas calculés dans la derniére boucle
B[0]=1+2*p-q*(1-u[0][0])
B[N]=1+2*p-q*(1-u[0][N])
A[0,1]=-p
A[0,0]=B[0]
A[N,N]=B[N]

#Boucle pour faire la resolution du sytéme linéaire à chaque instant ti
# puis stocker la i'éme résolution dans la i'éme ligne de la matrice u
for n in range(M):
    for i in range(N+1):#boucle pour recalculer les coefficients de la diagonale de la matrice A à l'instant ti
        B[i] = 1+2*p-q*(1-u[n][i])
        A[i,i]=B[i]
    # les conditions au bord pour l'espace à chaque instant ti
    u[n][0]=1/(1+np.exp(5*n*T/M))**2  #u(0,t) = 1/(1 + exp(5t))²  0 ≤ t ≤ T
    u[n][N]=1/(1+np.exp(1-(5*n*T)/M))**2 #u(1,t) = 1/(1 + exp(1-5t))²  0 ≤ t ≤ T
    u[n+1]=solve(A,u[n])

for i in range (M+1):
    u[i][0]=u[i][1]

u=np.transpose(u) # faire la transposé de la matrice u pour representer u(x,t)

#visualisation de la fonction temperature u
import matplotlib.pyplot as plt

fig = plt.figure()

plt.title("Plotting 1-D array")
plt.xlabel("t")
plt.ylabel("u(5,t) ")
plt.plot(t,u[:][5] , color = "red", marker = "o", label = "Array elements")
plt.legend()
plt.show()