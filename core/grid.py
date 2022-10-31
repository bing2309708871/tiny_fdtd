import numpy as np
import core.constants as const

from tqdm import tqdm
import os
from os import path, makedirs, chdir, remove
from glob import glob
from datetime import datetime

def curl_E(E) :
    curl = np.zeros(E.shape)

    curl[:, :-1, :, 0] += E[:, 1:, :, 2] - E[:, :-1, :, 2]
    curl[:, :, :-1, 0] -= E[:, :, 1:, 1] - E[:, :, :-1, 1]

    curl[:, :, :-1, 1] += E[:, :, 1:, 0] - E[:, :, :-1, 0]
    curl[:-1, :, :, 1] -= E[1:, :, :, 2] - E[:-1, :, :, 2]

    curl[:-1, :, :, 2] += E[1:, :, :, 1] - E[:-1, :, :, 1]
    curl[:, :-1, :, 2] -= E[:, 1:, :, 0] - E[:, :-1, :, 0]
    return curl


def curl_H(H) :
    curl = np.zeros(H.shape)

    curl[:, 1:, :, 0] += H[:, 1:, :, 2] - H[:, :-1, :, 2]
    curl[:, :, 1:, 0] -= H[:, :, 1:, 1] - H[:, :, :-1, 1]

    curl[:, :, 1:, 1] += H[:, :, 1:, 0] - H[:, :, :-1, 0]
    curl[1:, :, :, 1] -= H[1:, :, :, 2] - H[:-1, :, :, 2]

    curl[1:, :, :, 2] += H[1:, :, :, 1] - H[:-1, :, :, 1]
    curl[:, 1:, :, 2] -= H[:, 1:, :, 0] - H[:, :-1, :, 0]

    return curl

## FDTD Grid Class
class Grid:


    def __init__( self, shape, grid_spacing=0.05e-6, permittivity=1.0, permeability=1.0):
        self.grid_spacing = float(grid_spacing)
        self.Nx, self.Ny, self.Nz = shape
        self.D = int(self.Nx > 1) + int(self.Ny > 1) + int(self.Nz > 1)
        self.courant_number = 0.99 * float(self.D) ** (-0.5)
        self.time_step = self.courant_number * self.grid_spacing / const.SPEED_LIGHT
        self.E = np.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.H = np.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.inverse_permittivity = np.ones((self.Nx, self.Ny, self.Nz, 3)) / float(permittivity)
        self.inverse_permeability = np.ones((self.Nx, self.Ny, self.Nz, 3)) / float(permeability)

        self.time_steps_passed = 0
        self.sources = []

    def run(self, total_time, progress_bar = True):
        time = range(0, int(total_time), 1)
        if progress_bar:
            time = tqdm(time)
        for _ in time:
            self.step()

    def step(self):
        self.update_E()
        self.update_H()
        self.time_steps_passed += 1

    def update_E(self):
        curl = curl_H(self.H)
        self.E += self.courant_number * self.inverse_permittivity * curl
        for src in self.sources:
            src.update_E()

    def update_H(self):
        curl = curl_E(self.E)
        self.H -= self.courant_number * self.inverse_permeability * curl
        for src in self.sources:
            src.update_H()

    def add_source(self, name, source):
        source._register_grid(self)
        self.sources[name] = source


    def __setitem__(self, key, attr):
        x, y, z = key
        attr._register_grid(grid=self,x=x,y=y,z=z)