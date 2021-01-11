import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation 
from matplotlib.colors import LogNorm
from matplotlib import rc

def read_parameters(file):
	f = open(file, "r")
	l = f.readlines()
	x = {}
	for i in range(len(l)):
		it = iter(l[i].split())
		x.update(dict(zip(it,it)))
	return x

file = "DATA/Method.dat"

parameters = read_parameters(file)

method = parameters["Method"]
neqs = int(parameters["neqs"])
Nx   = int(parameters["Nx"])
Ny   = int(parameters["Ny"])
n_im = int(parameters["it"])
factor = int(parameters["factor"])

# Domain simulation
x = np.fromfile("DATA/x.dat", dtype="d", count=Nx)
y = np.fromfile("DATA/y.dat", dtype="d", count=Ny)

# scale of the image
x_m = round(np.max(x),1)
y_m = round(np.max(y),1)

# definition variable to animate
density    = 0
x_velocity = 1
y_velocity = 2
presion    = 3

# choose variable
variable = density

# array definition of Hidrodynamics variables
Prim = np.zeros((n_im, neqs, Nx, Ny))

# empty array to save images to simulate
ims = []

print("Method = ", method)
print("number equations = ", neqs)
print("system size = ", Nx)
print("Number if iterations = ", n_im)

print("")
if variable == 0:
	name = "density"
	print("Begins 2D animation of density")
if variable == 1:
	name = "x-velocity"
	print("Begins 2D animation of x-velocity")
if variable == 2:
	name = "y_velocity"
	print("Begins 2D animation of y-velocity")
if variable == 3:
	name = "presion"
	print("Begins 2D animation of presion")



fig, ax = plt.subplots()
fig.suptitle(variable, fontsize=10)
for i in range(int(n_im/factor)):
	# Lecture of each data at simulation time 
	Prim[i] = np.fromfile("DATA/Primitive"+str(i*factor)+".dat",
						  dtype="d",
						  count=neqs*Nx*Ny).reshape((neqs, Nx, Ny))
	 
	im = ax.imshow(Prim[i, variable, :, :].T, 
			           cmap="inferno", 
			           norm=LogNorm(), 
			           extent=[0, x_m, 0, y_m],
			           interpolation="gaussian"
			           )
	ax.set_xlabel("$x[ cm ]$")
	ax.set_ylabel("$y[ cm ]$")
	ax.set_aspect('equal')
	ims.append([im])
clb = fig.colorbar(im, extend="both")
clb.ax.set_title('$\\rho$')
clb.ax.set_title('$\\rho[\\frac{g}{cm^{3}}]$')
ani = animation.ArtistAnimation(fig, ims, 
								interval=1, 
								blit=False,
								repeat_delay=1000)

ani.save(name+"_"+method+".gif",fps=1000 ,writer="pillow")
# plt.show()

