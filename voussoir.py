"""
Este módulo implementa la solución de M.S. Diederichs y P.K. Kaiser* para el
análisis de la estabilidad de techo de excavaciones en macizos resistentes
estratificados con el modelo de dovelas.

*Stability of large excavations in laminated hard rock masses: the voussoir
analogue revisited.
International Journal of Rock Mechanics and Mining Sciences 36 (1999) 97-117
https://www.sciencedirect.com/science/article/pii/S0148906298001806?via%3Dihub

Realiza los cálculos para un rango dado del valor del espesor del estrato (t):
desde t hasta t x f_t

Tiene en cuenta una carga parabólica por encima del estrato.
Pueden implementar ustedes otras variantes fácilmente.

@author: Fernando García Bastante
Universidad de Vigo
"""
#  Comentar si no se desea cambiar el directorio de trabajo
import os
# Introduzca la ruta a su directorio de trabajo
W_D = r'C:\Users\usuario\SynologyDrive\Proyecto_colaborativo\rock_mechanics'
os.chdir(W_D)

# Se importan los módulos requeridos
import numpy as np
import pandas as pd

# #  Datos de entrada
s = 14               # Span (m)
t = 0.5              # Thickness (m)
e = 19900            # Young Module (MPa)
gamma = 28.9           # Density (kN/m3)
sigma_c = 95         # RCS (MPa)
alpha = 18           # Inclination (º)
fhi = 36             # Friction (º)
upp_density = 25.9    # Upper Layer Density (kN/m3)
upp_thicker = 1      # Thickness (m)

# Nombre de fichero de salida y espesor máximo a analizar (t  x  f_t)
out_file = "vouss_stab.xlsx"
f_t = 1.5

# Conversión previa a los cálculos de grados a radianes
alpha = alpha * np.pi/180
fhi = fhi * np.pi/180


def calc_charges(upp_density, upp_thicker, t):
    '''
    Función para calcular la carga (parabólica) debida a estratos superiores.
    Puede incorporar usted otras cargas: v.g. presión de sostenimiento.
'''
    # Carga de los estratos superiores (parabólica)
    subcharge_load_horiz = 7/9 * upp_density * upp_thicker / t

    # Carga total (según inclinación)
    gamma_e = (gamma + subcharge_load_horiz) * np.cos(alpha) / 1000

    return gamma_e


def z(thick):
    '''
Función que calcula fuerzas y factores de seguridad en el modelo de dovelas
'''
    # Lista para guardar resultados: [n, zo, z_chk, z, fm, fav]
    save_results = []
    # Nomino el espesor como t
    t = thick

    # Carga
    gamma_e = calc_charges(upp_density, upp_thicker, t)

    # Inicio del paso de n y del contador de casos de vuelco inestables
    step = 0.01
    buck_limit = 1

    #  Algoritmo que itera sobre n
    for n in np.arange(0.01, 1 + step, step):
        zo = t * (1 - 2/3 * n)
        L = s + 8 / (3 * s) * zo**2

        #  Iniciación de dl y de su variación entre iteraciones
        dl = 0
        inc_dl = 1  # Un número grande en el comienzo para que itere

        #  Iteración en dl
        while inc_dl > 0.0000001:
            z_chk = (8 / (3 * s) * zo**2 - dl)

            # Caso inestable: pasa al siguiente valor de n
            if z_chk < 0:
                buck_limit += 1
                # fm es igual a infinito para imposibilitar su selección
                fm = float('inf')
                break

            #  Iteraciones en incrementos de dl
            else:
                z = np.sqrt((3 * s / 8) * z_chk)
                fm = gamma_e * s**2 / (4*n*z)
                fav = fm / 3 * (2/3 + n)
                dlp = dl
                dl = fav/e * L
                inc_dl = dl - dlp

        #  Guardado de datos al salir de iterar en incremento de dl
        save_results.append(([n, zo, z_chk, z, fm, fav]))

    #  Resultados a array para obtener la solución con fm mínima
    results_array = np.array(save_results)
    # Determinación de las posiciones correspondientes a valores mínimos
    posicion = np.argmin(results_array, axis=0)
    # Selección de la fila con menor fm [índice 4 en la array]
    min_fm_row = np.append(results_array[posicion[4], :], buck_limit)
    min_fm = min_fm_row[4]

    #  Se calculan los coeficientes de seguridad y la deflexión
    fs_crush_by2 = (sigma_c/2) / min_fm
    fs_sliding = min_fm * min_fm_row[0] * np.tan(fhi)/(gamma_e * s)
    deflection = min_fm_row[1] - min_fm_row[3]
    fs_array = np.array([fs_crush_by2, fs_sliding, deflection])

    return np.append(min_fm_row, fs_array)


def sensib(thick=t, ft=f_t):
    '''
Función que varía el espesor t hasta un máximo de (t x f_t) y lanza para cada
uno de ellos a la función z.
'''
    #  Iniciación del almacén de resultados y creación del grid de valores de t
    resultados = np.empty((0, 10))
    t_grid = np.linspace(thick, ft * thick, 50)

    #  Se lanza cada t a la función z
    for i in range(len(t_grid)):
        sol = (z(t_grid[i]).reshape(1, -1))
        resultados = np.append(resultados, sol, axis=0)

    #  Guardado de resultados en un dataframe
    df = pd.DataFrame(resultados, index=t_grid,
                      columns=['n', 'zo', 'z_chk', 'z', 'fm', 'fav',
                               'buck_limit', 'fs_crush', 'fs_sliding',
                               'deflection (mm)'])

    #  Guardado en fichero excel
    df.to_excel(out_file)

    return df


# Se corren las funciones con los inputs definidos arriba
df = sensib()

#  Se pasan los espesores a una columna más y se grafican los resultados
df['t (m)'] = df.index
df.plot(x='t (m)', figsize=(9, 6), y=['fs_sliding', 'fs_crush'],
        color='Green', logx=False, subplots=True)
df.plot(x='t (m)', figsize=(9, 6), y=['deflection (mm)', 'buck_limit'],
        color='Green', logx=False, subplots=True)
