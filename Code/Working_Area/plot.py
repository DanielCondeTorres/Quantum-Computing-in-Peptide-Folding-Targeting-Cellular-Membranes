# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 11:11:50 2023

@author: DANIEL
"""
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

#def get_directory(file_path):
#    return os.path.dirname(os.path.abspath(file_path))

# Función para parsear los argumentos de línea de comandos
def parse_arguments():
    parser = argparse.ArgumentParser(description='Animación de puntos en 3D con plano fijo')
    parser.add_argument('-p', '--iteractions', required=True, help='Archivo de texto con la ecuación del plano')
    parser.add_argument('-x', '--frequencies', required=True, help='Archivo XYZ con las coordenadas de los puntos')
    return parser.parse_args()

# Parsea los argumentos de línea de comandos
args = parse_arguments()

#xyz_directory = get_directory(args.iteractions)

# Extrae los coeficientes del plano desde el archivo de texto

# Cargar el archivo
data = np.loadtxt(args.iteractions)

# Separar los datos en columnas
counts = data[:, 0]  # Primera columna
values = data[:, 1]  # Segunda columna

fig = plt.figure()
plt.plot(counts, values)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel("Conformation Energy",fontsize=20)
plt.xlabel("VQE Iterations",fontsize=20)
fig.add_axes([0.44, 0.51, 0.44, 0.32])
plt.plot(counts[40:], values[40:])

plt.ylabel("Conformation Energy", fontsize=18)
plt.xlabel("VQE Iterations",fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()
#BITSTRING


data = np.loadtxt(args.frequencies,dtype=str)
# Separar los datos en columnas
representamos = 50
energies = data[:, 0].astype(float)[-representamos:]#[-representamos:]  # Primera columna
bitstrings = data[:, 1][-representamos:] # Segunda columna
# Contar la frecuencia de cada bitstring
bitstring_counts = {}
# Contar la frecuencia de cada bitstring
total_bitstrings = len(bitstrings)
for bitstring in bitstrings:
    if bitstring in bitstring_counts:
        bitstring_counts[bitstring] += 1
    else:
        bitstring_counts[bitstring] = 1
        
        
# Calcular el total de frecuencias
total_frequencies = sum(bitstring_counts.values())
# Calcular los porcentajes
bitstring_percentages = {key: (count / total_frequencies) * 100 for key, count in bitstring_counts.items()}
# Crear el histograma
keys_as_strings = [str(int(key)) for key in bitstring_counts.keys()]
fig, ax = plt.subplots()
bars = ax.bar(bitstring_percentages.keys(),bitstring_percentages.values(), color='b')
ax.set_xlabel('Bitstrings',fontsize=18)
ax.set_ylabel('Population %',fontsize=18)

sorted_indices = np.argsort(energies)
# Obtener los 5 valores más pequeños
smallest_energies = np.copy(energies[sorted_indices[:5]])

# Agregar los valores de energía encima de las barras
for bar, energy in zip(bars, energies):
       if energy in smallest_energies:
         height = bar.get_height()
         ax.annotate(f'{energy:.2f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3), textcoords='offset points',
                ha='center', va='bottom',fontsize = 16)
# Rotar las etiquetas del eje x para que sean legibles
plt.xticks(rotation=90)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
# Mostrar el gráfico
plt.tight_layout()
plt.show()










#
# Crear un diccionario para almacenar las bitstrings más comunes para cada energía
representamos = 100
energies = data[:, 0].astype(float)[-representamos:]#[-representamos:]  # Primera columna
bitstrings = data[:, 1][-representamos:] # Segunda columna
energy_bitstrings = {}

# Calcular la bitstring más común para cada energía
# Calcular la bitstring más común para cada energía
for energy, bitstring in zip(energies, bitstrings):
    if energy in energy_bitstrings:
        current_most_common_bitstring = energy_bitstrings[energy]
        current_bitstring_count = bitstrings.tolist().count(current_most_common_bitstring)
        new_bitstring_count = bitstrings.tolist().count(bitstring)
        if new_bitstring_count > current_bitstring_count:
            energy_bitstrings[energy] = bitstring
    else:
        energy_bitstrings[energy] = bitstring

# Calcular el total de frecuencias
total_frequencies = len(energies)

# Calcular los porcentajes
energy_percentages = {key: (bitstrings.tolist().count(value) / total_frequencies) * 100 for key, value in energy_bitstrings.items()}

# Obtener la energía más baja y su bitstring asociado
lowest_energy = min(energy_bitstrings.keys())
lowest_energy_bitstring = energy_bitstrings[lowest_energy]

# Obtener la energía más alta y su bitstring asociado
highest_energy = max(energy_bitstrings.keys())
highest_energy_bitstring = energy_bitstrings[highest_energy]

# Crear el histograma
fig, ax = plt.subplots()
bars = ax.stem(energy_percentages.keys(), energy_percentages.values())

# Calcular la energía más baja y su bitstring asociado
lowest_energy_index = np.argmin(energies)
lowest_energy = energies[lowest_energy_index]
lowest_energy_bitstring = bitstrings[lowest_energy_index]
# Calcular la energía más alta y su bitstring asociado
highest_energy_index = np.argmax(energies)
highest_energy = energies[highest_energy_index]
highest_energy_bitstring = bitstrings[highest_energy_index]


maxima= max(energy_percentages.values())+0.5
minima =  min(energy_percentages.values())
# Anotar la energía más baja y su bitstring asociado
plt.annotate(f' Energy: {lowest_energy:.2f}\nBitstring: {lowest_energy_bitstring}', 
             (lowest_energy, maxima), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=12)

# Anotar la energía más alta y su bitstring asociado
plt.annotate(f' Energy: {highest_energy:.2f}\nBitstring: {highest_energy_bitstring}', 
             (highest_energy, minima), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=12)




ax.set_xlabel('Energy', fontsize=18)
ax.set_ylabel('Population %', fontsize=18)

plt.xticks(rotation=90)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

plt.show()

#########
#Energias: 
representamos = 100
energies = data[:, 0].astype(float)[-representamos:]#[-representamos:]  # Primera columna
bitstrings = data[:, 1][-representamos:] # Segunda columna
energy_counts={}
for energy in energies:
    if energy in energy_counts:
        energy_counts[energy] += 1
    else:
        energy_counts[energy] = 1
        
        
# Calcular el total de frecuencias
total_frequencies = sum(energy_counts.values())
# Calcular los porcentajes
energy_percentages = {key: (count / total_frequencies) * 100 for key, count in energy_counts.items()}
# Crear el histograma
fig, ax = plt.subplots()
bars = ax.stem(energy_percentages.keys(),energy_percentages.values())
ax.set_xlabel('Energy',fontsize=18)
ax.set_ylabel('Population %',fontsize=18)


plt.xticks(rotation=90)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
# Mostrar el gráfico
plt.tight_layout()
plt.show()
print('E')
plt.close()
