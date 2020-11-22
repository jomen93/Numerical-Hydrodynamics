import numpy as np
import matplotlib.pyplot as plt

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

neqs = int(parameters["neqs"])
Nx   = int(parameters["Nx"])

im = np.fromfile("DATA/Primitive1.dat",
				 dtype="d",
				 count=neqs*Nx).reshape((neqs, Nx))


