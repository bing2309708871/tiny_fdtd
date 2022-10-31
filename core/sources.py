import numpy as np

class PointSource():
    def __init__(self, period=15, amplitude=1.0, phase_shift=0.0, name=None):
        self.grid = None
        self.period = period
        self.amplitude = amplitude
        self.phase_shift = phase_shift
        self.name = name

    def _register_grid(self, grid, x, y, z):
        self.grid = grid
        self.grid.sources.append(self)
        setattr(grid, self.name, self)
        self.x, self.y, self.z = (x, y, z)
        self.period = int(float(self.period) / self.grid.time_step + 0.5)

    def update_E(self):
        q = self.grid.time_steps_passed
        src = self.amplitude * np.sin(2 * np.pi * q / self.period + self.phase_shift)
        self.grid.E[self.x, self.y, 0, 2] += src

    def update_H(self):
        pass