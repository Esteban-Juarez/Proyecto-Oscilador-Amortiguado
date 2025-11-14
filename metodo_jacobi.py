# Programa para resolver sistemas de ecuaciones lineales mediante el Metodo de Jacobi
import numpy as np

# Lectura de datos
with open("/home/estebanjrz/Nuevo-Semestre/Python/Datos_entrada.inp", "r") as file1:
    lineas1 = file1.readlines()

# Extrae los valores de n, m y epsilon
n, m, epsilon = map(float, lineas1[1].strip().split())
n, m = int(n), int(m)

with open("/home/estebanjrz/Nuevo-Semestre/Python/condiciones_iniciales.inp", "r") as file2:
    lineas2 = file2.readlines()

# Convierte la linea en el arreglo de aproximacion inicial
x0 = np.array(list(map(float, lineas2[0].strip().split())), dtype=float)

with open("/home/estebanjrz/Nuevo-Semestre/Python/matriz_coeficientes.inp", "r") as file3:
    lineas3 = file3.readlines()

# Convierte las lineas de texto en los elementos de la matriz extendida
A = np.array([list(map(float, line.split())) for line in lineas3], dtype=float)

# Inicializa la solucion
x = np.zeros(n)

# Contador de iteraciones
k = 1

while k <= m:
    for i in range(n):
        suma = np.dot(A[i, :-1], x0) - A[i, i] * x0[i]
        x[i] = A[i, n]/A[i, i] - suma/A[i, i]
    
    # Error relativo
    delta = np.linalg.norm(x - x0)/np.linalg.norm(x)

    # Criterio de convergencia
    if delta < epsilon:
        print(f"Terminado exitosamente en interaciones {k}")
        print("Con solucion:")
        for i, xi in enumerate(x, 1):
            print(f"x{i} = {xi:.8f}")
        exit()
    else:
        x0[:] = x[:]

    k = k + 1

print("Iteraciones excedidas")