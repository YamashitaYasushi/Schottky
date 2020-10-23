import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

# parameters : theta in Project 4.2
tau = np.pi/4

r = 1
theta = np.pi*(1/2)
# generators
st = np.sin(tau)
ct = np.cos(tau)
tt = np.tan(tau)
w  = r*np.exp(theta*(1j))
a = np.matrix([[1, ct*w], [ct/w, 1]])/st
b = np.matrix([[1, ct  ], [ct  , 1]])/st
B = LA.inv(b)
A = LA.inv(a)
gens = [a, B, A, b]
# initial disks: [center, radius] for a, B, A, b
D = [(w/ct,r*tt), (-1/ct,tt), (-w/ct,r*tt), (1/ct,tt)]
     
# cyclic order
co = [[3, 0, 1], [0, 1, 2], [1, 2, 3], [2, 3, 0]]

def plot_circle(c):
    ax.add_artist(plt.Circle((c[0].real, c[0].imag), c[1], alpha=0.2))

# p.84 BOX 10: mobius_on_circle
# M: mobius map = 2x2 matrix, i = 0, 1, 2, 3 --> D[i] is the circle
def mobius(M, i):
    z = D[i][0]-D[i][1]**2/(M[1,1]/M[1,0]+D[i][0]).conjugate()
    Q = (M[0,0]*z+M[0,1])/(M[1,0]*z+M[1,1])
    s = abs(Q-(M[0,0]*(D[i][0]+D[i][1])+M[0,1])/(M[1,0]*(D[i][0]+D[i][1])+M[1,1]))
    return (Q, s)

# depth first search 
# M : mobius map = 2x2 matrix, lastgen = 0:a 1:B, 2:A, 3
def dfs(M, lastgen, depth=0):
    for i in co[lastgen]:
        c = mobius(M, i)
        plot_circle(c)
        if c[1] > 0.01 and depth < 5:
            dfs(M*gens[i], i, depth+1)

# main
fig, ax = plt.subplots()
for i in range(4):
    plot_circle(D[i])
    dfs(gens[i], i)

plt.gca().set_aspect('equal', adjustable='box')
plt.xlim([-2.5, 2.5])
plt.ylim([-2.5, 2.5])
plt.show()
