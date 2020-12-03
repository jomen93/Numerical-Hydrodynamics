import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


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
n_im = int(parameters["it"])
factor = int(parameters["factor"])
# n_im = int(n_im)

print("Method = ", method)
print("number equations = ", neqs)
print("system size = ", Nx)
print("Number if iterations = ", n_im)


x = np.fromfile("DATA/x.dat",
				 dtype="d",
				 count=Nx).reshape(Nx)

print(np.shape(x))

im = np.fromfile("DATA/Primitive"+str(n_im)+".dat",
				 dtype="d",
				 count=neqs*Nx).reshape((neqs, Nx))

def an(var, var_name):
	# Set up the plotting area
	fig, ax = plt.subplots()
	# ax.set(xlim=(-0.1, 1.1), ylim=(0.2,2.3))

	# plot the first line
	im0 = np.fromfile("DATA/Primitive0.dat",
				 dtype="d",
				 count=neqs*Nx).reshape((neqs, Nx))
	line = ax.plot(x, im0[var, :], "b-", lw=2)[0]
	ax.grid(color='k',
	         alpha=0.5,
	         linestyle='dashed',
	         linewidth=0.9)
	ax.set_ylim(0,1.1)



	# create function to update the line
	def animate(i):
		im = np.fromfile("DATA/Primitive"+str(i*factor)+".dat",
				 dtype="d",
				 count=neqs*Nx).reshape((neqs, Nx))
		print("Read the Primitive"+str(i*factor)+".dat file in the animation")
		line.set_data(x, im[var, :])
		# ax.set_title(method+" $t$ = {:.2f} s".format(t[i]))
		ax.set_title(method)
		ax.grid(color='k',
	         alpha=0.5,
	         linestyle='dashed',
	         linewidth=0.9)
		ax.set_xlabel("$x$")
		ax.set_ylabel("$u$")
		ax.set_ylim(0,1.1)


	# Call FuncAnimation and show
	n_frames = int(n_im/factor)
	anim = FuncAnimation(fig, animate, interval=100, frames=n_frames)
	# anim.save(var_name+".gif",fps=1000 ,writer="pillow")
	plt.show()

an(0, "Density")
an(1, "Velocity")
an(2, "Presion")