### 2、添加点光源

#### 2.1 光源波形

​		本处点光源设置为正弦波，具有幅度、周期和相位三个参数。其中q表示经过的时间

```python
src = amplitude * sin(2 * pi * q / period + phase_shift)
E[x, y, z, 2] += src
```

#### 2.2 网格设置

​		仿真任意器件都需要设置网格大小。新建网格类如下，其中shape表示输入求解域三维的长度数据，为

```python
class Grid:
    def __init__(self, shape, grid_spacing=0.05e-6, permittivity=1.0, permeability=1.0):
        self.Nx, self.Ny, self.Nz = shape
		self.grid_spacing = float(grid_spacing)
        self.D = int(self.Nx > 1) + int(self.Ny > 1) + int(self.Nz > 1)
        self.courant_number = 0.99 * float(self.D) ** (-0.5)
        self.time_step = self.courant_number * self.grid_spacing / const.c
        self.E = np.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.H = np.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.inverse_permittivity = np.ones((self.Nx, self.Ny, self.Nz, 3)) / float(permittivity)
        self.inverse_permeability = np.ones((self.Nx, self.Ny, self.Nz, 3)) / float(permeability)

        self.time_steps_passed = 0
        self.sources = []
```

#### 2.3 更新电磁场

```python
    def run(self, total_time, progress_bar=True):
        time = range(0, int(total_time), 1)
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
```

2.4 设置属性

​		在python中__setitem__(self,key,value)方法时python魔法方法的一种，这个方法会让类按照一定的方法存储和key映射的value。该值可以使用另一种魔法方法__getitem__(self,key)来获取。

使用场景：当期望定义的类具备按照键存储值时，即类能够执行data['key']=value

目的：如果给类定义了__setitem__方法，则可以方便的给类进行赋值。


```python
    def __setitem__(self, key, attr):
        x, y, z = key
        attr._register_grid(grid=self,x=x,y=y,z=z)
```



#### 2.4 建立点光源类

```python
class PointSource:
    def __init__( self, period = 15, amplitude = 1.0, phase_shift = 0.0,name= None):
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

    def update_E(self):
        q = self.grid.time_steps_passed
        src = self.amplitude * sin(2 * pi * q / self.period + self.phase_shift)
        self.grid.E[self.x, self.y, 0, 2] += src

    def update_H(self):
        pass
```

#### 2.5 建立简单的test文件

```python
WAVELENGTH = 1550e-9
SPEED_LIGHT = 299_792_458.0  # [m/s] speed of light

N = 1
nstep = 200
grid = Grid((160, 120, N))

grid[80, 60, 1] = PointSource(period=WAVELENGTH / SPEED_LIGHT, name="source")

import matplotlib.pyplot as plt
import matplotlib

for i in range(nstep):
    grid.step()  # running simulation 1 timestep a time and animating
    if i % 5 == 0:
        plt.figure()
        norm = Norm(vmax=-0.1,vmin=0.1)
       plt.imshow(grid.E[:,:,0,2],interpolation='bilinear',cmap=plt.cm.jet,norm=norm)
        plt.colorbar()
        plt.show()
```

