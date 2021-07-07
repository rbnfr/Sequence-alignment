# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 17:46:10 2015

@author: Rubén
"""

import datetime
import time
import pprint
import pandas as pd
from os import path
import csv
from Bio.SeqIO.FastaIO import SimpleFastaParser

def obtener_matriz_similitud(input_file:str) -> dict:
    """
    Obtiene los datos de similitud de un fichero y devuelve una matriz de similitud entre bases nitrogenadas.
    :Fichero_entrada: Archivo que contiene la matriz de similitud entre bases.
    """

    _, file_type = path.splitext(input_file)

    if file_type == '.xlsx':
        score_excel = pd.read_excel(input_file, index_col=0)
        return score_excel.to_dict('dict')

    elif file_type == '.csv':
        with open(input_file, 'r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.readline())
            csv_separator = dialect.delimiter        
        
        score_csv = pd.read_csv(input_file, index_col=0, sep = csv_separator)
        return score_csv.to_dict('dict')

    elif file_type == '.txt':
        #Obtención de la lista de bases.
        with open(input_file, "r") as f:
            documento = f.readlines()
            bases = documento[0].split()[1:]
            matriz = {}

            # Incializar la matriz a cero.
            for base in bases:
                matriz[base] = {base : 0}

            for i in range(len(bases)):
                for j in range(len(bases)):
                    linea_punt = documento[j+1]
                    linea_punt_lista = linea_punt.split(" ")
                    puntuacion = linea_punt_lista[i+1].rstrip()
                    matriz[bases[i]][linea_punt_lista[0]] = puntuacion

            # print("Score matrix")
            # pp = pprint.PrettyPrinter(indent=4)
            # pp.pprint(matriz)
            f.close()
            return matriz
    else:
        raise Exception('File format not supported')

def analizar_fasta(file) -> dict:
    """
    Obtiene las secuencias y los identificadores de un archivo en formato FASTA y los devuelve en una tabla.
    """
    sequences_dict = {}

    try:
        with open(file, 'r') as fasta_file:
            fasta_parsed = SimpleFastaParser(fasta_file)
            for values in fasta_parsed:
                sequences_dict[values[0]] = sequences_dict.get(values[0], values[1])
            fasta_file.close()
    except IOError:
        """
        Si se detecta un error de entrada por no poner correctamente el nombre del archivo, su ruta o por no estar
        este en el mismo directorio que el script, lo informará y dará la opción de introducir los datos manualmente
        """
        print ("No se puede abrir el archivo, compruebe que el nombre o la ruta sean correctos.")
        identificador1 = input("Identificador de la secuencia 1: ")
        secuencia1     = input("Secuencia 1: ")
        identificador2 = input("Identificador de la secuencia 2: ")
        secuencia2     = input("Secuencia 2: ")        
        
        sequences_dict = {identificador1:secuencia1,identificador2:secuencia2}

    return sequences_dict
    # Trata de abrir el archivo con las secuencias
    """
    try:
        archivo_fasta = open(archivo)
        for linea in archivo_fasta:
            if linea.startswith(">"):
                secuencia = ""
                cabecera = linea.rstrip()
                # Dividimos la cadena por por el símbolo ">" y
                # tomamos la primera palabra a partir del primer carácter para eliminar el ">" inicial.
                identificador = cabecera.split(">")[1]
                tabla[identificador] = secuencia
            else:
                secuencia+=linea.rstrip()
                tabla[identificador]=secuencia

    # Si el archivo no puede ser abierto, se informa del error.
    # Se da la opción de introducir manualmente las secuencias para los casos de prueba
    except IOError:
        
        #Si se detecta un error de entrada por no poner correctamente el nombre del archivo, su ruta o por no estar
        #este en el mismo directorio que el script, lo informará y dará la opción de introducir los datos manualmente
        
        print ("No se puede abrir el archivo, compruebe que el nombre o la ruta sean correctos.")
        # secuencia1 = raw_input("Secuencia 1: ")
        # secuencia2 = raw_input("Secuencia 2: ")
        # identificador1 = raw_input("Identificador de la secuencia 1: ")
        # identificador2 = raw_input("Identificador de la secuencia 2: ")
        # tabla = {identificador1:secuencia1,identificador2:secuencia2}

    return tabla
    """

def generar_salida(tabla_secuencias,gap,matriz_similitud,archivo_salida) -> None:
    """
    Subrutina a la que se da una tabla con dos secuencias y sus identificadores, un archivo de salida y las puntuaciones para un alineamiento.
    Genera un fichero txt con los identificadores, las secuencias alineadas y la puntuación del alineamiento.
    Por defecto, sobreescribirá un archivo de salida si ya existe otro con el mismo nombre.
    """
    lista_secuencias = list(tabla_secuencias.values())
    lista_identificadores = list(tabla_secuencias.keys())
    alineamiento1,alineamiento2,puntuacion = alineamiento(lista_secuencias[0],lista_secuencias[1],gap,matriz_similitud)

    # Abrimos un fichero en el que escribiremos los alineamientos junto a los identificadores de las secuencias, así como la puntuación de ese alineamiento.
    # Cada vez que se ejecute el algoritmo, si se especica un nombre de fichero idéntico al anterior o a otro que ya exista en el mismo directorio, este será sobreescrito.
    # Para ir añadiendo alineamientos a un mismo fichero en caso de querer ejecutarlo de forma recursiva, cambiar la "w" de la siguiente línea por "a".
    with open(archivo_salida,"w") as output_file:
        output_file.write(str(lista_identificadores[0]))
        output_file.write(str("\n"))
        output_file.write(str(alineamiento1))
        output_file.write(str("\n"))
        output_file.write(str(alineamiento2))
        output_file.write(str("\n"))
        output_file.write(str(lista_identificadores[1]))
        output_file.write(str("\n"))
        output_file.write(str("Puntuación: "))
        output_file.write(str(puntuacion))
        output_file.write(str("\n\n"))
        output_file.close()

def alineamiento(secuencia1,secuencia2,gap,matriz_similitud):
    """
    Subrutina de alineamiento global de secuencias (algoritmo de Needleman-Wunsch).
    Alinea dos secuencias de la misma longitud y devuelve el alineamiento junto a la puntuación obtenida.
    La puntuación dependerá del valor que le hayamos definido a los eventos de coincidencia, error o gap.
    Devuelve cada una de las secuencias con los gaps (si los hubiera) y la puntuación de ese alineamiento.
    """

    # Tabla para asociar cada base a una posición dentro de la matriz de similitud
    bases = {'A':0,'T':1,'C':2,'G':3}

    # La matriz de similitud almacena los valores de coincidencia y error. Se pueden editar los valores para asignar puntuaciones concretas según el caso.

    # Variables para recorrer las secuencias
    I = range(len(secuencia1)+1)
    J = range(len(secuencia2)+1)

    # Generación de la matriz vacía
    F = [[0 for i in J] for j in I]

    # Generación de casos base
    for i in I:
        F[i][0] = gap*i
    for j in J:
        F[0][j] = gap*j

    # Cálculo de la matriz de puntuación
    for i in range(1, len(secuencia1)+1):
        for j in range(1, len(secuencia2)+1):
            # Bases
            base1 = secuencia1[i-1]
            base2 = secuencia2[j-1]

            # Posición en la matriz de similitud
            # posicion1 = bases[base1]
            # posicion2 = bases[base2]
            # prueba = matriz_similitud[base1][base2]
            # print "prueba:",prueba

            diagonal  = F[i-1][j-1] + int(matriz_similitud[base1][base2])
            izquierda = F[i-1][j]   + gap
            arriba    = F[i][j-1]   + gap
            F[i][j] = max(diagonal, izquierda, arriba)

    # Reconstrucción del alineamiento
    alineamiento1 = ""
    alineamiento2 = ""

    i = len(secuencia1)
    j = len(secuencia2)

    while (i > 0 and j > 0):
        # Puntuaciones
        puntuacion = F[i][j]
        punt_diagonal = F[i-1][j-1]
        punt_arriba = F[i][j-1]
        punt_izquierda = F[i-1][j]

        # Bases
        base1 = secuencia1[i-1]
        base2 = secuencia2[j-1]

        # Posición en la matriz de similitud
        # posicion1 = bases[base1]
        # posicion2 = bases[base2]

        if puntuacion == punt_diagonal + int(matriz_similitud[base1][base2]):
            alineamiento1 = secuencia1[i-1] + alineamiento1
            alineamiento2 = secuencia2[j-1] + alineamiento2
            i = i-1
            j = j-1

        elif puntuacion == punt_izquierda + gap:
            alineamiento1 = secuencia1[i-1] + alineamiento1
            alineamiento2 = "-" + alineamiento2
            i = i-1

        elif puntuacion == punt_arriba + gap:
            alineamiento1 = "-" + alineamiento1
            alineamiento2 = secuencia2[j-1] + alineamiento2
            j = j-1

    while i > 0:
        alineamiento1 = secuencia1[i-1] + alineamiento1
        alineamiento2 = "-" + alineamiento2
        i = i-1

    while j > 0:
        alineamiento1 = "-" + alineamiento1
        alineamiento2 = secuencia2[j-1] + alineamiento2
        j = j-1

    return alineamiento1, alineamiento2, F[len(secuencia1)][len(secuencia2)]

def main():

    # Definir archivos de entrada, salida y generar la tabla con cada identificador ligado a su secuencia
    archivo_fasta = input("Archivo FASTA de entrada:")
    while not archivo_fasta:
        archivo_fasta = input("Archivo FASTA de entrada:")
    
    archivo_similitud = input("Archivo con la matriz de similitud:")
    while not archivo_similitud:
        archivo_similitud = input("Archivo con la matriz de similitud:")

    archivo_salida_nombre = input("Archivo de salida:")
    while not archivo_salida_nombre:
        archivo_salida_nombre = input("Archivo de salida:")

    archivo_salida = archivo_salida_nombre + ".txt"
    tabla_secuencias = analizar_fasta(archivo_fasta)

    # Comprobación de las secuencias y los identificadores obtenidos
    #    for identificador in tabla_secuencias.keys():
    #        print "Identificador:",identificador,"\n","Secuencia:",tabla_secuencias[identificador],"\n"

    # Definición de las puntuaciones de coincidencia, error y gap
    # coincidencia = int(raw_input("Puntuación de coincidencia:"))
    # error = int(raw_input("Puntuación de error:"))

    gap = int(input("Puntuación de gap: "))
    matriz_similitud = obtener_matriz_similitud(archivo_similitud)

    ti=time.perf_counter() # Tiempo inicial

    # Se escribe el alineamiento, junto a los identificadores de cada secuencia y la puntuación obtenida en el archivo de salida
    generar_salida(tabla_secuencias,gap,matriz_similitud,archivo_salida)

    tf=time.perf_counter() # Tiempo final
    print ("ETA -->",tf-ti)


if __name__ == "__main__":
    main()
