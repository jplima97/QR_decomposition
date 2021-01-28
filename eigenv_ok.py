import numpy as np
import sys
import math 
sys.stdout = open("eigenvalues.txt", "w")

'''Program that calculates the eigenvalues ​​of an array of order n \nNOTE:
The txt file must have the values ​​separated by single space \nand for lines, a line must be skipped.'''

def UpperTriangularMatrix(newa):  
    newa = np.array(newa, float)
    n = len(newa)
    for i in range(n):
        for j in range(i):
            if newa[i][j]>0 and 1e-10 > newa[i][j] or newa[i][j] == 0:
                return True
            else:
                return False

def QRfactorization(A):
    A = np.array(A, float)
    n = len(A)
    I = np.identity(n)
    c = np.zeros([n,1])
    e = np.zeros([n,1])
    v = np.zeros([n,1])
    newC = np.zeros([n,1])
    newE = np.zeros([n,1])
    newV = np.zeros([n,1])    
    Q = np.zeros([n,n])
    H = np.zeros([n,n])
    R = np.zeros([n,n])
    newa = np.zeros([n,n]) 
    normac = np.zeros([n,1])
    
    def norm(c):
        nor = 0
        for i in range(n):
            nor = nor + c[i]**2
        norma = nor**(1/2)
        return norma
    
    #Step 1 =======================================================
    
    #Matrix c -----------------------------------------------------
    for i in range(n):
        c[i][0] = A[i][0]
    
    #Matrix e -----------------------------------------------------
    if c[0][0] > 0: 
        e[0][0] = 1
    if c[0][0] < 0:
        e[0][0] = -1
    # Norm of C --------------------------------------------------
    normac = norm(c)
    
    # Householder Matrix =========================================
    
    #Matrix v -----------------------------------------------------
    c2 = normac*e    
    v = c+c2
    
    # Householder ------------------------------------------------
    H = I - (2/np.dot(v.T,v))*np.outer(v,v.T)
    Q = H
    R = np.dot(H,A)
        
    # Resetting the matrices --------------------------------------
    c = newC
    e = newE
    v = newV

    # Step 2 ======================================================
    for i in range(1,n-1):
        for j in range(i,n):
            c[j] = R[j][i]
        if c[i] > 0:
            e[i] = 1
        if c[i] < 0:
            e[i] = -1

        normac = norm(c)
        v = c + normac*e
        H = I-(2/np.dot(v.T,v))*np.outer(v,v.T)
        Q = np.dot(Q,H)
        R = np.dot(H,R)   
    A = np.dot(R,Q)
    
    for i in range(n):
        for j in range(n):
            newa[i][j] = A[i][j]
    return newa

def read():
    f = open('qr.txt', 'r')
    temp = f.read().split('\n')
    MatrixA = []
    for i in temp:
        if i == '':
            continue
        MatrixA.append(i.split(" "))
    for i in range(len(MatrixA)):
        for j in range(len(MatrixA[i])):
            MatrixA[i][j] = float(MatrixA[i][j])
    return MatrixA
    return print(len(MatrixA))

def main():
    MatrixA = read()
    Matrix = np.array(MatrixA)
    print('\nTyped Matrix: \n')
    print(Matrix)

    n = len(MatrixA)
    eigenvalues = np.zeros((n,n))
    condition = False
    
    while condition == False: 
        newa = QRfactorization(MatrixA)
        condition = UpperTriangularMatrix(newa)
        
        if condition == False:
            MatrixA = newa
        if condition == True:
            print("\nWOW! We have found the upper triangular matrix!")
 
    if condition == True:
        print("\nUpper triangular matrix: \n")
        print(newa)
        
        # Print the matrix with the eigenvalues on the main diagonal-
        for i in range(n):
            for j in range(n):
                eigenvalues[i][i] = newa[i][i]
        print("\nEigenvalue Matrix: \n")
        print(eigenvalues)
        
        # Method for printing the highest eigenvalue ----------------
        for i in range(n-1):
            for j in range(n-1):
                if eigenvalues[0][0]<eigenvalues[i][j]:
                    value = eigenvalues[i][j]
                else:
                    value = eigenvalues[0][0]   
        print('\n Highest eigenvalue: {}'.format(value))
             
main() 