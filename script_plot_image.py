import math
import copy
import matplotlib.pyplot as plt
import os
from matplotlib.image import imread

def cum_sum(A):
    if A and isinstance(A[0], list):
        aux = [val for row in A for val in row]
    else:
        aux = A
    if not aux:
        return []
    rez = [0] * len(aux)
    rez[0] = aux[0]
    for i in range(1, len(aux)):
        rez[i] = rez[i-1] + aux[i] 
    return rez

def mean_3d(matrix, axis=None):
    Z = len(matrix)
    Y = len(matrix[0]) if Z > 0 else 0
    X = len(matrix[0][0]) if Y > 0 else 0

    if axis == -1:
        axis = 2

    if axis is None:
        toate_elementele = [matrix[z][y][x] for z in range(Z) for y in range(Y) for x in range(X)]
        return sum(toate_elementele) / len(toate_elementele) if toate_elementele else 0

    if axis == 0:
        return [[sum(float(matrix[z][y][x]) for z in range(Z)) / Z for x in range(X)] for y in range(Y)]

    if axis == 1:
        return [[sum(float(matrix[z][y][x]) for y in range(Y)) / Y for x in range(X)] for z in range(Z)]

    if axis == 2:
        return [[sum(float(matrix[z][y][x]) for x in range(X)) / X for y in range(Y)] for z in range(Z)]

    raise ValueError("Axa trebuie sa fie 0, 1, 2, -1 sau None.")

def rang(A):
    n = len(A)
    m = len(A[0])
    rang=0
    for k in range(0,m-1):
        max=0
        index=k
        for a in range(k,n):
            if abs(A[a][k])>max:
                max=abs(A[a][k])
                index=a
        if index !=k:
            for a in range(0,m):
                A[k][a],A[index][a]=A[index][a], A[k][a]
        for i in range(k+1,n):
            v = A[i][k] / A[k][k]
            for j in range(k,m):
                A[i][j] -=v*A[k][j]
    for i in range(0,n):
        if A[i][i] !=0: rang+=1
    return rang

def identitate(n):
    idn = [[0] * i + [1] + [0] * (n - i - 1) for i in range(n)]
    return idn

def block(blocuri):
    m = []
    for rand_blocuri in blocuri:
        for rand_u in zip(*rand_blocuri):
            rand_complet = []
            for rand in rand_u:
                rand_complet.extend(rand)
            m.append(rand_complet)
    return m

def zeros(randuri, coloane):  
    return [[0 for j in range(coloane)] for i in range(randuri)]

def transpusa(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

def scadere(A,B):
    return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def inmultire_rapida_J(M, J, i, j):
    c = J[i][i]
    s = J[i][j]
    
    for k in range(len(M)):
        val_i = M[k][i]
        val_j = M[k][j]
        M[k][i] = val_i * c - val_j * s
        M[k][j] = val_i * s + val_j * c
        
    return M

def inmultire_matrice(A, B):
    n = len(A)
    m = len(A[0])
    p = len(B[0])
    rez = [[0] * p for _ in range(n)]
    for i in range(n):
        for j in range(p):
            for k in range(m):
                rez[i][j] += A[i][k] * B[k][j]
    return rez

def inmultire_scalar(A, k):
    return [[k*A[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def house_holder(a):
    n,m = len(a), len(a[0])
    r = copy.deepcopy(a)
    q = identitate(n)
    limita = min(n - 1, m)
    for i in range(limita):
        x = [[r[k][i]] for k in range(i,n)]
        norma_x = math.sqrt(sum(x[k][0]**2 for k in range(len(x))))
        alpha = -norma_x if x[0][0] < 0 else norma_x
        u = copy.deepcopy(x)
        u[0][0] -= alpha
        norma_u = math.sqrt(sum([u[k][0]**2 for k in range(len(u))]))
        v = inmultire_scalar(u,(1/norma_u))
        Q_sub = scadere((identitate(n-i)),(inmultire_scalar((inmultire_matrice(v, transpusa(v))),2)))
        
        if i > 0:
            Qn = block([
                [identitate(i), zeros(i,n-i)],
                [zeros(n-i,i), Q_sub]
            ])
        else:
            Qn = Q_sub
            
        r = inmultire_matrice(Qn, r)
        q = inmultire_matrice(q, transpusa(Qn))
        return q,r

def calcul_rotatie(B, i, j):
    if abs(B[i][j]) < 1e-15:
        return [[1 if r == c else 0 for c in range(len(B[0]))] for r in range(len(B))]

    tau = (B[j][j] - B[i][i]) / (2.0 * B[i][j])
    
    if abs(tau) > 1e10:
        c = 1.0
        s = 0.0
    else:
        if tau >= 0:
            t = 1.0 / (abs(tau) + math.sqrt(1.0 + tau**2))
        else: 
            t = -1.0 / (abs(tau) + math.sqrt(1.0 + tau**2))
            
        c = 1.0 / math.sqrt(1.0 + tau**2)
        s = t * c
        
    J = [[1 if r == c else 0 for c in range(len(B[0]))] for r in range(len(B))]
    J[i][i] = c
    J[j][j] = c
    J[i][j] = s
    J[j][i] = -s
    
    return J

def svd(a):
    N = len(a[0])    
    M_len = len(a)   
    
    V = [[1.0 if i == j else 0.0 for j in range(N)] for i in range(N)]
    U = [[0.0 for j in range(M_len)] for i in range(M_len)]
    sigma = [[0.0 for j in range(N)] for i in range(M_len)]
    
    done = False
    treceri = 0
    
    while not done and treceri < 20:
        done = True
        treceri += 1
        print(f"-> SVD lucreaza: Treapta {treceri}/20...")
        
        for i in range(N):
            for j in range(N):
                if i < j:
                    b_ii = sum(a[k][i] ** 2 for k in range(M_len))
                    b_jj = sum(a[k][j] ** 2 for k in range(M_len))
                    b_ij = sum(a[k][i] * a[k][j] for k in range(M_len))
                    
                    if abs(b_ij) > 1e-9 * math.sqrt(abs(b_ii * b_jj)) and abs(b_ij) > 1e-15:
                        tau = (b_jj - b_ii) / (2.0 * b_ij)
                        
                        if abs(tau) > 1e10: 
                            c, s = 1.0, 0.0
                        else:
                            t = 1.0 / (abs(tau) + math.sqrt(1.0 + tau**2))
                            if tau < 0: t = -t
                            c = 1.0 / math.sqrt(1.0 + t**2)
                            s = t * c
                        
                        for k in range(M_len):
                            val_i = a[k][i]
                            val_j = a[k][j]
                            a[k][i] = val_i * c - val_j * s
                            a[k][j] = val_i * s + val_j * c
                            
                        for k in range(N):
                            val_i = V[k][i]
                            val_j = V[k][j]
                            V[k][i] = val_i * c - val_j * s
                            V[k][j] = val_i * s + val_j * c
                            
                        done = False 

    print("-> SVD a terminat rotatiile! Calculez U si Sigma...")

    for i in range(N):
        norm_col = math.sqrt(sum(a[k][i] ** 2 for k in range(M_len)))
        if i < M_len:
            sigma[i][i] = norm_col
        
        if norm_col > 1e-10:
            for k in range(M_len):
                U[k][i] = a[k][i] / norm_col
                
    return U, V, sigma

plt.rcParams['figure.figsize'] = [16,8]

A = imread('picture.jpg')
X = mean_3d(A, 2)

print("1")
X_inversat = [[256 - val for val in rand] for rand in X]
img = plt.imshow(X_inversat)
img.set_cmap('gray')
plt.axis('off')
plt.show()

print("2")
pas = 8 
X_mic = [X[i][::pas] for i in range(0, len(X), pas)]

print("3")
U, V, S = svd(X_mic)
V = transpusa(V)

nr_svd = min(len(S), len(S[0]))
sigma = [S[k][k] for k in range(nr_svd) if S[k][k] != 0]
sigma.sort(reverse=True)
print("Calcul svd complet")
print("4")
j=0
for r in (5,20,100):    
    r_efectiv = min(r, len(U[0])) 
    
    Xapprox = inmultire_matrice(
        inmultire_matrice(
            [[U[i][c] for c in range(r_efectiv)] for i in range(len(U))],
            [[S[i][c] for c in range(r_efectiv)] for i in range(r_efectiv)]        
        ), 
        [[V[i][c] for c in range(len(V[0]))] for i in range(r_efectiv)] 
    )
    plt.figure(j+1)
    j +=1
    Xapprox_inversat = [[256 - val for val in rand] for rand in Xapprox]
    img = plt.imshow(Xapprox_inversat)
    img.set_cmap('gray')
    plt.axis('off')
    plt.title('r = ' + str(r_efectiv))
    plt.show()
    
plt.figure(j+1)
plt.semilogy(sigma)
plt.title('Singular Values')
plt.show()

plt.figure(j+2)
suma_totala = sum(sigma)
plt.plot(cum_sum([val / suma_totala for val in sigma]))
plt.title('Singular Values: Cumulative Sum')
print("5")
plt.show()