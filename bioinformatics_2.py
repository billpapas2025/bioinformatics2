import streamlit as st
from Bio import SeqIO, pairwise2
from Bio.SeqUtils import molecular_weight, gc_fraction
import pandas as pd
import re
import io  # Importar el módulo io

# Función para contar las bases nitrogenadas
def count_bases(sequence):
    base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base in sequence:
        if base in base_counts:
            base_counts[base] += 1
    return base_counts

# Función para realizar alineamiento de secuencias
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments

# Función para buscar motivos en una secuencia
def search_motif(sequence, motif):
    return [match.start() for match in re.finditer(motif, sequence)]

# Función principal para la interfaz de usuario
def main():
    st.title('Análisis de Secuencias de ADN')
    st.sidebar.title('Opciones')

    # Cargar archivo FASTA
    uploaded_file = st.sidebar.file_uploader("Cargar archivo FASTA", type=["fasta"])

    if uploaded_file is not None:
        st.sidebar.text('Archivo cargado exitosamente.')
        # Asegurarse de leer el archivo como texto
        uploaded_file_text = uploaded_file.read().decode("utf-8")
        records = list(SeqIO.parse(io.StringIO(uploaded_file_text), "fasta"))

        if len(records) == 0:
            st.sidebar.text('El archivo no contiene secuencias válidas.')
            return

        sequence = str(records[0].seq)

        # Mostrar información básica de la secuencia
        st.header('Información de la Secuencia')
        st.write(f'**ID:** {records[0].id}')
        st.write(f'**Descripción:** {records[0].description}')
        st.write(f'**Longitud:** {len(sequence)} bases')

        # Contar bases nitrogenadas
        base_counts = count_bases(sequence)
        df_counts = pd.DataFrame([base_counts], columns=base_counts.keys())
        
        st.header('Conteo de Bases')
        st.write(df_counts)

        # Calcular porcentaje de bases
        total_bases = len(sequence)
        percentage_counts = {base: (count / total_bases) * 100 for base, count in base_counts.items()}
        df_percentage = pd.DataFrame([percentage_counts])

        st.header('Porcentaje de Bases')
        st.write(df_percentage)

        # Calcular contenido GC
        gc_content = gc_fraction(sequence) * 100
        st.header('Contenido GC')
        st.write(f'**Contenido GC:** {gc_content:.2f}%')

        # Calcular índice de masa molecular
        mol_weight = molecular_weight(sequence, seq_type='DNA')
        st.header('Índice de Masa Molecular')
        st.write(f'**Índice de Masa Molecular:** {mol_weight:.2f} Da')

        # Búsqueda de motivos
        st.header('Búsqueda de Motivos')
        motif = st.text_input('Introduce el motivo a buscar (ejemplo: ATG):')
        if motif:
            motif_positions = search_motif(sequence, motif)
            st.write(f'Motivo encontrado en las posiciones: {motif_positions}')
        
        # Alineamiento de secuencias
        if len(records) > 1:
            st.header('Alineamiento de Secuencias')
            seq1 = str(records[0].seq)
            seq2 = str(records[1].seq)
            alignments = align_sequences(seq1, seq2)
            st.write(f'Alineamientos entre la primera y segunda secuencia:')
            for alignment in alignments:
                st.text(pairwise2.format_alignment(*alignment))

    else:
        st.sidebar.text('Carga un archivo FASTA para comenzar.')

if __name__ == '__main__':
    main()
