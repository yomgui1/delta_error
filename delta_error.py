# Delta Printer Error 3D Visualizer
#
# Copyright (C) 2018 Guillaume Roguez
#
# Author: Guillaume Roguez <yomgui1@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.
#

from math import *

class Delta:
    def __init__(self, L, r):
        self.L = float(L)
        self.r = float(r)
        self.ae = [0, 0, 0]
        self.reset()

    def reset(self):
        At = radians(90 + self.ae[0])
        Bt = radians(210 + self.ae[1])
        Ct = radians(330 + self.ae[2])
        self.ColPos = [
            [self.r*cos(At), self.r*sin(At), 0],
            [self.r*cos(Bt), self.r*sin(Bt), 0],
            [self.r*cos(Ct), self.r*sin(Ct), 0]
        ]

    def angular_error(self, *err):
        self.ae = err
        self.reset()

    def world_to_delta(self, *world):
        L2 = self.L**2

        XA = world[0] - self.ColPos[0][0]
        XB = world[0] - self.ColPos[1][0]
        XC = world[0] - self.ColPos[2][0]

        YA = world[1] - self.ColPos[0][1]
        YB = world[1] - self.ColPos[1][1]
        YC = world[1] - self.ColPos[2][1]

        return sqrt(L2 - XA**2 - YA**2), \
            sqrt(L2 - XB**2 - YB**2), \
            sqrt(L2 - XC**2 - YC**2)

    def delta_to_world(self, *delta):
        p1 = [self.ColPos[0][0], self.ColPos[0][1], delta[0]]
        p2 = [self.ColPos[1][0], self.ColPos[1][1], delta[1]]
        p3 = [self.ColPos[2][0], self.ColPos[2][1], delta[2]]

        p12 = [a-b for a,b in zip(p2, p1)]
        p13 = [a-b for a,b in zip(p3, p1)]

        d = sqrt(sum(a**2 for a in p12))

        ex = [a/d for a in p12]

        i = sum(a*b for a,b in zip(ex, p13))
        iex = [i*a for a in ex]

        eyv = [a-b for a,b in zip(p13, iex)]
        eyd = sqrt(sum(a**2 for a in eyv))
        ey = [a/eyd for a in eyv]

        ez = [ex[1]*ey[2] - ex[2]*ey[1],
              ex[2]*ey[0] - ex[0]*ey[2],
              ex[0]*ey[1] - ex[1]*ey[0]]

        j = sum(a*b for a,b in zip(ey, p13))

        x = d/2
        y = ((i**2 + j**2)/2 - i*x)/j
        z = sqrt(self.L**2 - x**2 - y**2)

        return tuple(round(a + x*b + y*c - z*d, 2) for a,b,c,d in zip(p1, ex, ey, ez))


def compute_surface(axe, delta_virtual, delta_real):
    # Create the mesh in polar coordinates and compute corresponding Z
    r = np.linspace(0, 1, 50)
    p = np.linspace(0, 2*np.pi, 50)
    R, P = np.meshgrid(r, p)
    Z = np.zeros(R.shape)

    # Express the mesh in the cartesian system.
    X, Y = R*np.cos(P), R*np.sin(P)

    # Make data.
    s = Z.shape
    for i in range(s[0]):
        for j in range(s[1]):
            x = delta_real.r * X[i,j]
            y = delta_real.r * Y[i,j]
            res = delta_virtual.world_to_delta(x, y, 0)
            z = delta_real.delta_to_world(*res)[2]
            Z[i,j] = z

    # Plot.
    return axe.plot_surface(X, Y, Z, cmap=plt.cm.coolwarm_r,
                            linewidth=0, antialiased=False)


if __name__ == "__main__":
    ## 3D visualization based on matplotlib.org examples

    # My delta settings, change them for your version
    rod_length = 368.
    col_radius = 163.

    # Delta configuration as seen by the computer
    my_virtual_delta = Delta(rod_length, col_radius)

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import pyplot as plt
    import numpy as np

    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.suptitle("Rods length: %umm, Columns radius: %umm" % (rod_length, col_radius))

    #===================
    # +1mm error on rods
    #===================
    my_real_delta = Delta(rod_length + 1, col_radius)
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.set_title("+1mm error on rods")

    surf = compute_surface(ax, my_virtual_delta, my_real_delta)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    ax.set_zlim(-2.00, 2.00)
    ax.zaxis.set_major_locator(plt.LinearLocator(10))
    ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

    #============================
    # -1mm error on column radius
    #============================
    my_real_delta = Delta(rod_length, col_radius - 1)

    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ax.set_title("-1mm error on column radius")

    surf = compute_surface(ax, my_virtual_delta, my_real_delta)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    # Customize the z axis.
    ax.set_zlim(-2.00, 2.00)
    ax.zaxis.set_major_locator(plt.LinearLocator(10))
    ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

    #=====================================
    # 1 degree error on first column angle
    #=====================================
    my_real_delta = Delta(rod_length, col_radius)
    my_real_delta.angular_error(1, 0, 0)

    ax = fig.add_subplot(2, 2, 3, projection='3d')
    ax.set_title("1 degree error on first column angle")

    surf = compute_surface(ax, my_virtual_delta, my_real_delta)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    # Customize the z axis.
    ax.set_zlim(-2.00, 2.00)
    ax.zaxis.set_major_locator(plt.LinearLocator(10))
    ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

    plt.tight_layout()
    plt.show()
