from core import *

import matplotlib.pyplot as plt
import matplotlib

WAVELENGTH = 1550e-9
N = 1
nstep = 300


grid = Grid((160, 120, N))
grid[80,60,0] = PointSource( period=WAVELENGTH / SPEED_LIGHT, name="source")

for i in range(nstep):
    grid.step()  # running simulation 1 timestep a time and animating
    if i % 5 == 0:
        plt.figure()
        norm = matplotlib.colors.Normalize(vmax=0.05, vmin=-0.05)
        plt.imshow(grid.E[:,:,0,2], interpolation="sinc", cmap=plt.cm.jet, norm=norm)
        plt.colorbar()
        plt.show()
