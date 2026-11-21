# OPTIMIZADOR DE RECURSOS - Distribución de horas técnicas
# Trabajo Práctico Integrador AE4 (Programación I - Matemática Discreta - Ing. de Requerimientos)
# Integrantes: Lucas Nahuel Bianchi - Fabricio Sabugo Robles (Comisión A)
#
# Modelo: sistema de n ecuaciones lineales con n incógnitas (x1 = Redes, x2 = Bases de Datos,
# x3 = Desarrollo). Se guarda como matriz aumentada n x (n+1) en una lista de listas,
# se reduce con Gauss-Jordan y se clasifica según las filas nulas que quedan.
#
# Formato del archivo de entrada: una ecuación por línea, coeficientes y término
# independiente separados por comas. Ejemplo: 1,1,1,100  ->  x1 + x2 + x3 = 100

import os

# Cambiamos la carpeta de trabajo a la carpeta donde está este archivo .py,
# así el programa encuentra los .txt de prueba aunque se ejecute desde VS Code.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

EPSILON = 0.0000001   # tolerancia: cualquier |valor| menor a esto se toma como cero


# ------------------------------------------------------------------------------
# LECTURA Y VALIDACIÓN DE DATOS (RF01 / RF02)
# Devuelve (matriz, n) si los datos son válidos, o (None, 0) si hay algún error.
# ------------------------------------------------------------------------------
def cargar_matriz(nombre_archivo):
    matriz = []
    try:
        archivo = open(nombre_archivo, "r", encoding="utf-8")
        lineas = archivo.readlines()
        archivo.close()

        if len(lineas) == 0:
            print("Error: El archivo esta vacio.")
            return None, 0

        columnas_esperadas = -1

        for linea in lineas:
            linea = linea.strip()
            if linea != "":
                partes = linea.split(",")
                fila = []
                for valor in partes:
                    fila.append(float(valor))      # si no es número salta ValueError

                # Todas las filas deben tener la misma cantidad de columnas
                if columnas_esperadas == -1:
                    columnas_esperadas = len(fila)
                elif len(fila) != columnas_esperadas:
                    print("Error: Las filas no tienen la misma cantidad de columnas.")
                    return None, 0

                matriz.append(fila)

        n = len(matriz)

        # Validamos que sea matriz aumentada (n ecuaciones por n+1 columnas)
        if n <= 0 or n + 1 != columnas_esperadas:
            print("Error: La cantidad de ecuaciones no coincide con las incognitas.")
            return None, 0

        return matriz, n

    except FileNotFoundError:
        print("Error: No se encontro el archivo:", nombre_archivo)
        return None, 0
    except ValueError:
        print("Error: El archivo contiene datos que no son numeros.")
        return None, 0


# ------------------------------------------------------------------------------
# OPERACIONES ELEMENTALES DE FILA (una función por operación del método) (RF03)
# ------------------------------------------------------------------------------
def intercambiar_filas(matriz, n, f1, f2):
    # Fi <-> Fj
    for j in range(n + 1):
        aux = matriz[f1][j]
        matriz[f1][j] = matriz[f2][j]
        matriz[f2][j] = aux


def multiplicar_fila(matriz, n, fila, escalar):
    # Fi <- k * Fi
    for j in range(n + 1):
        matriz[fila][j] = matriz[fila][j] * escalar


def sumar_filas(matriz, n, f_origen, f_destino, escalar):
    # Fdestino <- Fdestino + k * Forigen
    for j in range(n + 1):
        matriz[f_destino][j] = matriz[f_destino][j] + (matriz[f_origen][j] * escalar)


# ------------------------------------------------------------------------------
# MÉTODO DE GAUSS-JORDAN (RF03)
# ------------------------------------------------------------------------------
def gauss_jordan(matriz, n):
    fila = 0   # fila donde va el próximo pivote (solo avanza cuando se encuentra uno)
    for col in range(n):
        if fila < n:
            # Si el pivote es cero, buscamos (búsqueda secuencial) una fila inferior
            # que tenga un valor distinto de cero en esa columna y las intercambiamos
            if abs(matriz[fila][col]) < EPSILON:
                encontro = False
                k = fila + 1
                while k < n and not encontro:
                    if abs(matriz[k][col]) > EPSILON:
                        intercambiar_filas(matriz, n, fila, k)
                        encontro = True
                    k = k + 1

            pivote = matriz[fila][col]

            # Si hay pivote: lo llevamos a 1 y hacemos ceros arriba y abajo.
            # Si no hay (columna toda en cero), esa variable queda libre y
            # la misma fila se usa para la columna siguiente.
            if abs(pivote) > EPSILON:
                multiplicar_fila(matriz, n, fila, 1.0 / pivote)

                for k in range(n):
                    if k != fila and abs(matriz[k][col]) > EPSILON:
                        factor = -matriz[k][col]
                        sumar_filas(matriz, n, fila, k, factor)

                fila = fila + 1


# ------------------------------------------------------------------------------
# CLASIFICACIÓN DEL SISTEMA (RF03)
# Una fila 0 0 ... 0 | c con c != 0 es la ecuación imposible 0 = c -> INCOMPATIBLE.
# Una fila 0 0 ... 0 | 0 indica una ecuación que sobra -> COMPATIBLE INDETERMINADO.
# Se revisan TODAS las filas: si aparece una sola fila imposible, el sistema es
# incompatible aunque antes haya aparecido una fila 0 | 0.
# ------------------------------------------------------------------------------
def clasificar_sistema(matriz, n):
    estado = "COMPATIBLE DETERMINADO"
    for i in range(n):
        todos_ceros = True
        for j in range(n):
            if abs(matriz[i][j]) > EPSILON:
                todos_ceros = False

        if todos_ceros:
            if abs(matriz[i][n]) > EPSILON:
                return "INCOMPATIBLE"
            estado = "COMPATIBLE INDETERMINADO"   # seguimos revisando las demás filas

    return estado


# ------------------------------------------------------------------------------
# SALIDA POR PANTALLA Y A ARCHIVO (RF04 / RF05)
# ------------------------------------------------------------------------------
def mostrar_matriz(matriz, n):
    for i in range(n):
        linea = ""
        for j in range(n + 1):
            num = round(matriz[i][j], 2)
            if num == 0:
                num = 0.0          # evita que se muestre -0.0

            if j == n:
                linea = linea + "| " + str(num)
            else:
                linea = linea + str(num) + "\t"
        print(linea)


def guardar_reporte(nombre_salida, matriz, n, estado):
    try:
        archivo = open(nombre_salida, "w", encoding="utf-8")
        archivo.write("REPORTE DE OPTIMIZACION DE RECURSOS\n")
        archivo.write("===================================\n\n")
        archivo.write("Estado del sistema: " + estado + "\n\n")
        archivo.write("Matriz escalonada reducida:\n")

        for i in range(n):
            linea = ""
            for j in range(n + 1):
                num = round(matriz[i][j], 2)
                if num == 0:
                    num = 0.0
                linea = linea + str(num) + "\t"
            archivo.write(linea + "\n")

        archivo.write("\nAsignacion optima:\n")
        if estado == "COMPATIBLE DETERMINADO":
            for i in range(n):
                valor = round(matriz[i][n], 2)
                archivo.write(" -> Variable x" + str(i + 1) + ": " + str(valor) + " horas\n")
        elif estado == "COMPATIBLE INDETERMINADO":
            archivo.write("El sistema tiene infinitas soluciones posibles.\n")
        else:
            archivo.write("El sistema no tiene solucion (restricciones contradictorias).\n")

        archivo.close()
        print("\nReporte guardado en:", nombre_salida)
    except OSError:
        print("\nError: No se pudo guardar el reporte en", nombre_salida)


# ------------------------------------------------------------------------------
# ALGORITMO PRINCIPAL
# ------------------------------------------------------------------------------
def main():
    print("=== OPTIMIZADOR DE HORAS DE SOPORTE ===")
    nombre_archivo = input("Ingrese el nombre del archivo de restricciones (ej. prueba1.txt): ")

    matriz, n = cargar_matriz(nombre_archivo)

    if matriz is not None:
        print("\nMatriz inicial:")
        mostrar_matriz(matriz, n)

        gauss_jordan(matriz, n)

        print("\nMatriz escalonada reducida:")
        mostrar_matriz(matriz, n)

        estado = clasificar_sistema(matriz, n)
        print("\nEstado del sistema:", estado)

        if estado == "COMPATIBLE DETERMINADO":
            print("Asignacion optima:")
            for i in range(n):
                valor = round(matriz[i][n], 2)
                print(" -> Variable x" + str(i + 1) + ":", valor, "horas")
        elif estado == "COMPATIBLE INDETERMINADO":
            print("El sistema tiene infinitas soluciones: falta una restriccion independiente.")
        else:
            print("No existe asignacion posible: las restricciones se contradicen.")

        guardar_reporte("resultado_optimizacion.txt", matriz, n, estado)


if __name__ == "__main__":
    main()
