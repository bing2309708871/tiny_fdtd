""" often used constants go here """

from math import pi

SPEED_LIGHT: float = 299_792_458.0  # [m/s] speed of light
""" speed of light """

eps0: float = 4e-7 * pi
""" vacuum permittivity """

mu0: float = 1.0 / (eps0 * SPEED_LIGHT ** 2)
""" vacuum permeability """
