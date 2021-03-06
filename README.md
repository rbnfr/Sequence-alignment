# Sequence-alignment
**Input:** Sequences, similitude matrix.

**Output:** Similitude file in txt format with aligned sequences, gaps and total punctuation according to our similitude matrix.

## Methods used
### Obtener matriz de similitud
**Input:** Txt file with similitude matrix.

**Output:** similitude table.

```python
def obtener_matriz_similitud(fichero_entrada):    

    #Obtención de la lista de bases.
    with open(fichero_entrada, "r") as f:
        documento = f.readlines()
        bases = documento[0].split()[1:]
        matriz = {}

        # Incializar la matriz a cero.
        for base in bases:
            matriz[base] = {base : 0}

        for i in xrange(len(bases)):
            for j in xrange(len(bases)):
                linea_punt = documento[j+1]
                linea_punt_lista = linea_punt.split(" ")
                puntuacion = linea_punt_lista[i+1].rstrip()
                matriz[bases[i]][linea_punt_lista[0]] = puntuacion

        print matriz
        return matriz
```

### Analizar FASTA
**Input:** FASTA file with DNA sequences and headers.

**Output:** Table with headers as keys and sequences as values.

```python
def analizar_fasta(archivo):
    """
    Obtiene las secuencias y los identificadores de un archivo en formato FASTA y los devuelve en una tabla.
    """
    tabla={}

    # Trata de abrir el archivo con las secuencias
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
        """
        Si se detecta un error de entrada por no poner correctamente el nombre del archivo, su ruta o por no estar
        este en el mismo directorio que el script, lo informará y dará la opción de introducir los datos manualmente
        """
        print "No se puede abrir el archivo, compruebe que el nombre o la ruta sean correctos."
    #    secuencia1 = raw_input("Secuencia 1: ")
    #    secuencia2 = raw_input("Secuencia 2: ")
    #    identificador1 = raw_input("Identificador de la secuencia 1: ")
    #    identificador2 = raw_input("Identificador de la secuencia 2: ")
    #    tabla = {identificador1:secuencia1,identificador2:secuencia2}

    return tabla
```
### Generar salida
**Input:**
* Sequences table.
* Similitude table.
* GAP penalization (User input)
* Output filename (User input)

**Output:** Txt file with headers, alignment and punctuation.

```python
def generar_salida(tabla_secuencias,gap,matriz_similitud,archivo_salida):
    """
    Subrutina a la que se da una tabla con dos secuencias y sus identificadores, un archivo de salida y las puntuaciones para un alineamiento.
    Genera un fichero txt con los identificadores, las secuencias alineadas y la puntuación del alineamiento.
    Por defecto, sobreescribirá un archivo de salida si ya existe otro con el mismo nombre.
    """
    lista_secuencias = tabla_secuencias.values()
    lista_identificadores = tabla_secuencias.keys()
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
```

### Alineamiento
**Input:**
* Sequences
* GAP punctuation
* Similitude table

**Output:** Alignments and alignment punctuation.

```python
def alineamiento(secuencia1,secuencia2,gap,matriz_similitud):
    """
    Subrutina de alineamiento global de secuencias (algoritmo de Needleman-Wunsch).
    Alinea dos secuencias de la misma longitud y devuelve el alineamiento junto a la puntuación obtenida.
    La puntuación dependerá del valor que le hayamos definido a los eventos de copincidencia, error o gap.
    Devuelve cada una de las secuencias con los gaps (si los hubiera) y la puntuación de ese alineamiento.
    """

    # Tabla para asociar cada base a una posición dentro de la matriz de similitud
    bases = {'A':0,'T':1,'C':2,'G':3}

    # La matriz de similitud almacena los valores de coincidencia y error. Se pueden editar los valores para asignar puntuaciones concretas según el caso.

    # Variables para recorrer las secuencias
    I = xrange(len(secuencia1)+1)
    J = xrange(len(secuencia2)+1)

    # Generación de la matriz vacía
    F = [[0 for i in J] for j in I]

    # Generación de casos base
    for i in I:
        F[i][0] = gap*i
    for j in J:
        F[0][j] = gap*j

    # Cálculo de la matriz de puntuación
    for i in xrange(1, len(secuencia1)+1):
        for j in xrange(1, len(secuencia2)+1):
            # Bases
            base1 = secuencia1[i-1]
            base2 = secuencia2[j-1]

            # Posición en la matriz de similitud
#            posicion1 = bases[base1]
#            posicion2 = bases[base2]
#            prueba = matriz_similitud[base1][base2]
#            print "prueba:",prueba

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
#        posicion1 = bases[base1]
#        posicion2 = bases[base2]

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
```
