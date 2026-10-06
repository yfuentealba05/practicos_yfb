from Bio import AlignIO

# 1. Cargar el archivo de alineamiento en formato FASTA
# Asegúrate de poner el nombre correcto de tu archivo
NOMBRE_ARCHIVO = r"C:\Users\ACER\Desktop\taller1\pregunta1\align_pregunta1.fasta"

alignment = AlignIO.read(NOMBRE_ARCHIVO, "fasta")
total_sitios = alignment.get_alignment_length()

# Confirmación de longitud del alineamiento (Debe dar 1299)
print(f"--- INFORMACIÓN DEL ALINEAMIENTO ---")
print(f"Secuencias analizadas: {len(alignment)}")
print(f"Longitud total detectada: {total_sitios} sitios\n")

# 2. Identificar columnas 100% conservadas (estricta identidad)
conserved_positions = []
for i in range(total_sitios):
    # Convertir a mayúsculas e ignorar guiones de gaps
    col = [res.upper() for res in alignment[:, i]]
    
    # Comprobar si todas las secuencias coinciden en el mismo aminoácido (sin vacíos)
    if len(set(col)) == 1 and col[0] not in ['-', 'X', '.']:
        conserved_positions.append(i + 1)  # Posición basada en 1 (de 1 a 1299)

# 3. Agrupar posiciones consecutivas para formar bloques/regiones
regions = []
if conserved_positions:
    current_region = [conserved_positions[0]]
    for pos in conserved_positions[1:]:
        if pos == current_region[-1] + 1:
            current_region.append(pos)
        else:
            regions.append(current_region)
            current_region = [pos]
    regions.append(current_region)

# 4. Filtrar regiones con MÁS DE 3 aminoácidos (> 3 AA, es decir, de 4 AA o más)
valid_regions = [r for r in regions if len(r) > 3]

# 5. RESULTADO FINAL
x = len(valid_regions)

print("=" * 40)
print(f"Hay {x} regiones conservadas")
print("=" * 40 + "\n")

# Detalle de cada región encontrada dentro de los 1299 sitios
if x > 0:
    print("Detalle de las regiones (coordenadas entre 1 y 1299):")
    for idx, reg in enumerate(valid_regions, 1):
        inicio = reg[0]
        fin = reg[-1]
        longitud = len(reg)
        secuencia = str(alignment[0].seq[inicio - 1:fin]).upper()
        print(f"Región {idx}: Posiciones {inicio} a {fin} ({longitud} AA) -> Secuencia: {secuencia}")