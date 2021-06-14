import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

from NonLinear.chaos import inverse

import warnings


def colorFader(c1, c2, mix=0):
    # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)


warnings.filterwarnings("ignore", category=FutureWarning)

fig = plt.figure()
ax = plt.axes(projection="3d")

c1 = "#1f77b4"  # blue
c2 = "black"


n = 8 #max number of iterations 
eta = 0.88
P = 0.56 
xx, yy = np.meshgrid(np.linspace(0.01, 2, 3), np.linspace(0.01, 2, 3)) 
# meshgrid with evenly distributed points


def plot(Uo, Vo, Wo):
    ul, vl, wl = [], [], []

    for itr in range(n):
        if np.iscomplex(Wo) or np.iscomplex(Uo) or np.iscomplex(Vo):
            return

        if Uo ** 2 + Vo ** 2 + Wo ** 2 > 1:
            ax.scatter3D(Uo, Vo, Wo, color="red")
            # print(Uo, Vo, Wo, (Uo ** 2 + Vo ** 2 + Wo ** 2 + 1) / 2)
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

ax.set_xlabel('u')
ax.set_ylabel('v')
ax.set_zlabel('w')
plt.show()

