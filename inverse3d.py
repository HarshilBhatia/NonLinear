import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

from helper.Inverse import inverse

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def colorFader(c1, c2, mix=0):
    # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)



P = 0.9
eta = 0.88
xx, yy = np.meshgrid(np.linspace(0.01, 2, 5), np.linspace(0.01, 2, 5))

fig = plt.figure()
ax = plt.axes(projection="3d")

c1 = "#1f77b4"  
c2 = "red"
n = 15


# Plots the stepwise evolution of points on the bloch sphere
def plot(Uo, Vo, Wo):
    ul, vl, wl = [], [], []
    for itr in range(15):
        if np.iscomplex(Wo) or np.iscomplex(Uo) or np.iscomplex(Vo):
            return

        if Uo ** 2 + Vo ** 2 + Wo ** 2 > 1:
            ax.scatter3D(Uo, Vo, Wo, color="red")
            ul.append(Uo)
            vl.append(Vo)
            wl.append(Wo)
            ax.plot3D(ul, vl, wl, color="gray")
            return

        elif Uo ** 2 + Vo ** 2 + Wo ** 2 <= 1:
            ax.scatter3D(Uo, Vo, Wo, color=colorFader(c1, c2, itr / n))
            ul.append(Uo)
            vl.append(Vo)
            wl.append(Wo)

        u, v, w = inverse(Uo, Vo, Wo, eta, positive= False)
        Uo, Vo, Wo = u, v, w

    ax.plot3D(ul, vl, wl, color="gray")


for row in range(len(xx)):
    Z_row = []
    for col in range(len(xx[row])):
        u = (
            2
            * ((2 * P - 1) ** 0.5)
            * xx[row][col]
            / (xx[row][col] ** 2 + yy[row][col] ** 2 + 1)
        )
        v = (
            2
            * ((2 * P - 1) ** 0.5)
            * yy[row][col]
            / (xx[row][col] ** 2 + yy[row][col] ** 2 + 1)
        )
        w = (
            2 * ((2 * P - 1) ** 0.5) / (xx[row][col] ** 2 + yy[row][col] ** 2 + 1)
            - (2 * P - 1) ** 0.5
        )

        if u != 0:
            plot(u, v, w)


#Plot sphere
u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color="#d3d3d3")


plt.show()
