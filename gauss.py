# Resolución de sistemas de ecuaciones líneales por el método de Gauss  
import numpy as np

# open sirve para abrir ficheros (imagenes, documentos, texto) en nuestra computadora.
# la instruccion with provoca el cierre del fichero cuando se ejecuta el programa
# sintaxis de open: nombre del archivo; r --> read
with open("/home/estebanjrz/Nuevo-Semestre/Python/matriz_coeficientes.inp", "r") as file:
    lines = file.readlines()

# Determina el numero de ecuaciones a partir del numero de renglones n
n = len(lines)

# Convierte las lineas de texto de coefficients.inp en los elementos de A
A = np.array([list(map(float, line.split())) for line in lines], 
             dtype=float)

# Ciclo sobre las columnas
for i in range(n - 1):
    for p in range(i, n):
        if A[p, i] != 0:
            if p == i:
                break
            else:
                A[[i, p]] = A[[p, i]]
                break
        else:
            print("El sistema no tiene solucion unica")
            exit()
    
# Transformacion y reduccion del sistema
    for j in range(i + 1, n):
        m = A[j, i] / A[i, i]
        A[j, i:] -= m * A[i, i:]

# Verifica si el sistema tiene solucion unica
if A[n - 1, n - 1] == 0:
    print("El sistema no tiene solucion unica")
    exit()

# Sustitucion regresiva
x = np.zeros(n)
x[n - 1] = A[n - 1, n] / A[n - 1, n - 1]
for i in range(n - 2, -1, -1):
    sum_ax = np.dot(A[i, i + 1:n], x[i + 1:n])
    x[i] = (A[i, n] - sum_ax) / A[i, i]

# Muestra los resultados en pantalla
print("Calculo terminado satisfactoriamente")
print("Matriz final:")
print(A)
print("solucion:")
for i, xi in enumerate(x, 1):
    print(f"x{i} = {xi:.8f}")