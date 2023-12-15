import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import imageio
import os
import argparse
from matplotlib.lines import Line2D
from mayavi import mlab
# Función para obtener el directorio del archivo
def get_directory(file_path):
    return os.path.dirname(os.path.abspath(file_path))

# Función para parsear los argumentos de línea de comandos
def parse_arguments():
    parser = argparse.ArgumentParser(description='Animación de puntos en 3D con plano fijo')
    parser.add_argument('-p', '--plano', required=True, help='Archivo de texto con la ecuación del plano')
    parser.add_argument('-x', '--xyz', required=True, help='Archivo XYZ con las coordenadas de los puntos')
    return parser.parse_args()

# Parsea los argumentos de línea de comandos
args = parse_arguments()

# Obtiene el directorio del archivo XYZ
xyz_directory = get_directory(args.xyz)

# Extrae los coeficientes del plano desde el archivo de texto
with open(args.plano, 'r') as plano_file:
    lines = plano_file.readlines()

if len(lines) < 5:
    print("El archivo no contiene suficientes líneas para calcular la ecuación del plano.")
    exit(1)
A = float(lines[1].strip())
B = float(lines[2].strip())
C = float(lines[3].strip())
D = float(lines[4].strip())

print('C',C)
# Lee el archivo XYZ (asegúrate de que el archivo tenga el formato adecuado)
with open(args.xyz, 'r') as xyz_file:
    lines = xyz_file.readlines()
# Límites de los ejes y puntos del plano (ajusta estos valores según tus preferencias)
xlim = (-10, 10)
ylim = (-10, 10)
zlim = (-10, 10)
x_plane_points, y_plane_points = np.meshgrid(
    np.linspace(xlim[0], xlim[1], 100),
    np.linspace(ylim[0], ylim[1], 100)
)

if C!=0:
  z_plane_points = (-A * x_plane_points - B * y_plane_points - D) / C
# Lee el archivo .xyz (asegúrate de que el archivo tenga el formato adecuado)



# Listas para almacenar las coordenadas x, y, z y las etiquetas de cada frame
frames = []
frame = {'coords': [], 'labels': []}

# Variable para el análisis del archivo
leyendo_coordenadas = False
valores_x=[]
valores_y=[]
valores_z=[]
aminoacido=[]
for line in lines:
    # Ignora líneas en blanco
    if line.strip() == '':
        continue

    # Si encontramos un número, indica el inicio de un nuevo conjunto de coordenadas
    if line.strip().isdigit():
        if frame['coords']:
            frames.append(frame)
        frame = {'coords': [], 'labels': []}
        leyendo_coordenadas = True
        continue

    # Si estamos leyendo coordenadas, agrega las coordenadas y etiquetas a este frame
    if leyendo_coordenadas:
        parts = line.strip().split()
        x, y, z, label = float(parts[1]), float(parts[2]), float(parts[3]), parts[0]
        valores_x.append(x)
        valores_y.append(y)
        valores_z.append(z)
        aminoacido.append(label)
        frame['coords'].append([x, y, z])
        frame['labels'].append(label)

# Agrega el último frame a la lista
if frame['coords']:
    frames.append(frame)

# Configura la figura 3D de Matplotlib
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Evolución de puntos')

# Función para inicializar el gráfico
def init():
    return []
# Diccionario  para almacenar los colores del diccionario de aminoacidos
aminoacid_info = {
    'R': {'code': 'Arg', 'color': 'blue','color_mayavi': (0, 0, 1)},
    'C': {'code': 'Cys', 'color': 'orange', 'color_mayavi': (1, 0.5, 0)},
    'G': {'code': 'Gly', 'color': 'cyan', 'color_mayavi': (0, 1, 1)},
    'E': {'code': 'Glu', 'color': 'magenta', 'color_mayavi': (1, 0, 1)},
    'F': {'code': 'Phe', 'color': 'yellow', 'color_mayavi': (1, 1, 0)},
    'D': {'code': 'Asp', 'color': 'red', 'color_mayavi': (1, 0, 0)},
    'H': {'code': 'His', 'color': 'pink', 'color_mayavi': (1, 0.5, 0.5)},
    'I': {'code': 'Ile', 'color': 'brown', 'color_mayavi': (0.5, 0.25, 0)},
    'K': {'code': 'Lys', 'color': 'teal', 'color_mayavi': (0, 0.5, 0.5)},
    'L': {'code': 'Leu', 'color': 'olive', 'color_mayavi': (0, 0.5, 0)},
    'M': {'code': 'Met', 'color': 'gold','color_mayavi': (1, 0.8, 0)},
    'N': {'code': 'Asn', 'color': 'lime', 'color_mayavi': (0, 1, 0)},
    'P': {'code': 'Pro', 'color': 'cyan', 'color_mayavi': (0, 1, 1)},
    'Q': {'code': 'Gln', 'color': 'purple', 'color_mayavi': (0.5, 0, 0.5)},
    'A': {'code': 'Ala', 'color': 'pink', 'color_mayavi': (1, 0.5, 0.5)},
    'W': {'code': 'Trp', 'color': 'brown', 'color_mayavi': (0.5, 0.25, 0)},
    'T': {'code': 'Thr', 'color': 'teal', 'color_mayavi': (0, 0.5, 0.5)},
    'V': {'code': 'Val', 'color': 'olive', 'color_mayavi': (0, 0.5, 0)},
    'S': {'code': 'Ser', 'color': 'gold', 'color_mayavi': (1, 0.8, 0)},
    'X': {'code': 'Unknown', 'color': 'brown', 'color_mayavi': (0.5, 0.25, 0)},
    'Y': {'code': 'Tyr', 'color': 'navy', 'color_mayavi': (0, 0, 0.5)}
    # Agrega más aminoácidos y colores según sea necesario
}





# Convert amino acid labels to RGB colors
colors = [aminoacid_info[aa]['color_mayavi'] for aa in aminoacido]

# Convert the list of RGB tuples to a Nx3 numpy array
colors_array = np.array(colors)

# Create the visualization with colored points
mlab.figure(bgcolor=(1, 1, 1))

# Display the colored points, connect them with lines, and label them
for i in range(len(valores_x)):
    mlab.points3d(valores_x[i], valores_y[i], valores_z[i], mode='sphere', color=colors[i], scale_factor=1., opacity=1,resolution=10)#opacity 0.7
   # mlab.text3d(valores_x[i], valores_y[i], valores_z[i], aminoacido[i], scale=0.3, color=(0,0,0))
    if i < len(valores_x) - 1:
        mlab.plot3d(valores_x[i:i+2], valores_y[i:i+2], valores_z[i:i+2], tube_radius=0.03, color=(0., 0., 0),opacity=0.6)


# Create the visualization
if C!=0:
  mlab.mesh(x_plane_points/2, y_plane_points/2, z_plane_points/2, color=(0.5, 0., 0.5), opacity=0.7)  # You can change the color if needed


# Función para actualizar el gráfico en cada frame
def update(frame_num):
    ax.clear()
    used_aminoacids = {} ## Crea un diccionario para almacenar los aminoácidos utilizados
    coords = frames[frame_num]['coords']
    labels = frames[frame_num]['labels']
    x, y, z = zip(*coords)

    #ax.scatter(x, y, z, c='r', marker='o', s=100)

    for i in range(len(coords) ):
        #color = aminoacid_colors.get(labels[i], 'gray')  # Usa 'gray' como color predeterminado si no se encuentra la letra
        letter = labels[i]  # Letra de un solo carácter
        info = aminoacid_info.get(letter, {'code': 'Unknown', 'color': 'gray'})  # Si la letra no se encuentra, usa 'Unknown' y 'gray'
        code = info['code']  # Código de tres letras
        color = info['color']  # Color
        ax.scatter(x[i], y[i], z[i], c=color, marker='o', s=100,alpha= 0.5)
        used_aminoacids[letter] = (code, color)
    for i in range(len(coords) - 1):
        ax.plot([x[i], x[i + 1]], [y[i], y[i + 1]], [z[i], z[i + 1]], color='r')
        #used_aminoacids[labels[i]] = color
    for i, label in enumerate(labels):
        ax.text(x[i], y[i], z[i], label, fontsize=12)

    # Dibuja el plano fijo
    x_plane_points, y_plane_points = np.meshgrid(
        np.linspace(-2.5, 3, 100),
        np.linspace(-2.5, 3, 100)
        )
    if C!=0:
      z_plane_points = (-A * x_plane_points - B * y_plane_points - D) / C

      ax.plot_surface(
        x_plane_points, y_plane_points, z_plane_points,
        alpha=0.5, cmap='viridis', rstride=100, cstride=100
        )
    # Crea la leyenda utilizando los datos almacenados
    #legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label) for label, color in used_aminoacids.items()]
    # Crea la leyenda utilizando solo las tres letras y el color correspondiente
    legend_elements = [Line2D([0], [0], marker='o', color=color,linestyle='' ,markerfacecolor=color, markersize=10, label=f'{letter} ({code})') for letter, (code, color) in used_aminoacids.items()]

    ax.legend(handles=legend_elements, loc='upper right')

    # Agrega la leyenda al gráfico
    ax.add_artist(ax.legend(handles=legend_elements, loc='upper right'))


    ax.set_xlabel('Axis X')
    ax.set_ylabel('Axis Y')
    ax.set_zlabel('Axis Z')
    ax.set_title(f'Conformation {frame_num + 1}/{len(frames)}')
    ax.grid(False) # New one
    ax.axis('off')
    return []

from matplotlib.animation import FuncAnimation
ani = FuncAnimation(fig, update, frames=len(frames), init_func=init, blit=True)


plt.show()

