from Bio import SeqIO

archivo_entrada = r"C:\Users\ACER\Desktop\taller1\pregunta2\align_pregunta2.fasta"
archivo_salida = r"C:\Users\ACER\Desktop\taller1\pregunta2\alineamiento2_limpio.fasta"

# 1. Diagnóstico: Ver cuántas 'X' tiene CADA secuencia
print("--- DIAGNÓSTICO DE 'X' POR SECUENCIA ---")
registros = list(SeqIO.parse(archivo_entrada, "fasta"))

for rec in registros:
    count_x = str(rec.seq).upper().count('X')
    print(f"{rec.id.ljust(20)}: {count_x} 'X'")

# 2. Filtrado: Eliminar las secuencias de Plutella (Pxyl) o con cualquier 'X'
# Opción A: Eliminar directamente las secuencias que contengan "Pxyl" en el ID
# Opción B: Eliminar si tienen más de 0 'X' (ajusta el número si deseas)

UMBRAL_MAX_X = 0  # Permite máximo 0 'X' (elimina cualquier secuencia con 'X')

secuencias_limpias = []
for rec in registros:
    count_x = str(rec.seq).upper().count('X')
    
    # Filtro doble: ni que contenga 'Pxyl' ni que tenga 'X' desmedidas
    if "Pxyl" not in rec.id and count_x <= UMBRAL_MAX_X:
        secuencias_limpias.append(rec)
    else:
        print(f"\n--> Eliminando: {rec.id} ({count_x} 'X')")

# 3. Guardar el nuevo FASTA
SeqIO.write(secuencias_limpias, archivo_salida, "fasta")
print(f"\n¡Listo! Guardado '{archivo_salida}' con {len(secuencias_limpias)} secuencias.")