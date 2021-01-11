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
density       = 0
x_velocity    = 1
y_velocity    = 2
presion       = 3
# derived variables
norm_velocity = 5
Temperature   = 6

# choose variable
# when using a job file to compile all plots
# variable = sys.argv[0]
# manual plot
variable = Temperature

# array definition of Hidrodynamics variables
Prim = np.zeros((n_im, neqs, Nx, Ny))

# empty array to save images to simulate
ims = []

print("Method = ", method)
print("number equations = ", neqs)
print("system size = ", Nx)
print("Number if iterations = ", n_im)

print("")

if variable == density:
	name       = "density" 
	title_plot = '$\\rho$'
	units_plot = '$\\rho[\\frac{g}{cm^{3}}]$'
	print("Begins 2D animation of density")
if variable == x_velocity:
	name     = "x-velocity" 
	title_plot = '$v_{x}$'
	units_plot = '$v_{x}[\\frac{m}{s}]$'
	print("Begins 2D animation of x-velocity")
if variable == y_velocity:
	name      = "y_velocity"
	title_plot = '$v_{y}$'
	units_plot = '$v_{y}[\\frac{m}{s}]$'
	print("Begins 2D animation of y-velocity")
if variable == presion:
	name     = "presion"
	title_plot = '$P$'
	units_plot = '$P[Pa]$'
	print("Begins 2D animation of presion")
if variable == norm_velocity:
	name     = "Norm Velocity"
	title_plot = '$|v|$'
	units_plot = '$v[\\frac{m}{s}]$'
	print("Begins 2D animation of norm velocity")
if variable == Temperature:
	name     = "Temperature"
	title_plot = '$T$'
	units_plot = '$T[^{o}K]$'
	print("Begins 2D animation of Temperature")

fig, ax = plt.subplots()
fig.suptitle(name, fontsize=10)
interpolation = "gaussian"

for i in range(int(n_im/factor)):
	# Lecture of each data at simulation time 
	Prim[i] = np.fromfile("DATA/Primitive"+str(i*factor)+".dat",
						  dtype="d",
						  count=neqs*Nx*Ny).reshape((neqs, Nx, Ny))

	if variable == density:
		im = ax.imshow(Prim[i, variable, :, :].T, 
				           cmap="Blues",
				           norm=LogNorm(), 
				           extent=[0, x_m, 0, y_m],
				           interpolation=interpolation
				           )
	if variable == x_velocity:
		im = ax.imshow(Prim[i, variable, :, :].T, 
			           cmap="seismic",
			           norm=None, 
			           extent=[0, x_m, 0, y_m],
			           interpolation=interpolation
			           )
		
	if variable == y_velocity:
		im = ax.imshow(Prim[i, variable, :, :].T, 
			           cmap="seismic",
			           norm=None, 
			           extent=[0, x_m, 0, y_m],
			           interpolation=interpolation
			           )
		
	if variable == presion:
		im = ax.imshow(Prim[i, variable, :, :].T, 
			           cmap="inferno",
			           norm=LogNorm(), 
			           extent=[0, x_m, 0, y_m],
			           interpolation=interpolation
			           )
		
	if variable == norm_velocity:
		im = ax.imshow(np.sqrt(Prim[i, x_velocity, :, :].T**2 + Prim[i, y_velocity, :, :].T**2), 
			           cmap="seismic",
			           norm=None, 
			           extent=[0, x_m, 0, y_m],
			           interpolation=interpolation
			           )
		
	if variable == Temperature:
		im = ax.imshow(Prim[i, presion, :, :].T/Prim[i, density, :, :].T, 
			           cmap="jet",
			           norm=None, 
			           extent=[0, x_m, 0, y_m],
			           interpolation=interpolation
			           )
		
	ax.set_xlabel("$x[ cm ]$")
	ax.set_ylabel("$y[ cm ]$")
	ax.set_aspect('equal')
	ims.append([im])
clb = fig.colorbar(im, extend="both")
clb.ax.set_title(title_plot)
clb.ax.set_title(units_plot)
ani = animation.ArtistAnimation(fig, ims, 
								interval=1, 
								blit=False,
								repeat_delay=1000)

ani.save(name+"_"+method+".gif",fps=1000 ,writer="pillow")
# plt.show()

