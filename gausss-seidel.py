import numpy as np

# Leer los datos desde el archivo de entrada 
with open("/home/estebanjrz/Nuevo-Semestre/Python/Datos_entrada.inp", "r") as file1:
    lineas1 = file1.readlines()

# Extraer los valores de n, m y epsilon
n, m, epsilon = map(float, lineas1[1].strip().split())
n, m = int(n), int(m)

with open("/home/estebanjrz/Nuevo-Semestre/Python/condiciones_iniciales.inp", "r") as file2:
    lineas2 = file2.readline().strip()

# Convierte la linea en el arreglo de aproximacion inicial
x0 = np.array(list(map(float, lineas2.split())), dtype=float)

with open("/home/estebanjrz/Nuevo-Semestre/Python/matriz_coeficientes.inp", "r") as file:
    lineas = file.readlines()

# Convierte las lineas de texto en los elementos de la matriz extendida
A = np.array([list(map(float, line.split())) for line in lineas], 
             dtype=float)

# Aseguramiento que los argumentos de la matriz tenga n renglones y n+1 columnas
if A.shape[1] != n + 1:
    print("Error: la matriz aumentada debe tener n filas y n + 1 columnas.")
    exit()

x = np.zeros(n)

# Inicializa el contador de iteraciones
k = 1

while k < m:
    for i in range(n):
        suma = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:n], x0[i+1:])
        x[i] = A[i,n]/A[i,i] - suma/A[i,i]

    delta = np.linalg.norm(x - x0) / np.linalg.norm(x)
    
    if delta < epsilon:
        print(f"Terminado exitosamente en iteraciones {k}")
        print("Con solucion")
        for i, xi in enumerate(x, 1):
            print(f"x{i} = {xi:.8f}")
        exit()
    else:
        x0[:] = x[:] 
    k = k + 1

print("Iteraciones excedidas")
