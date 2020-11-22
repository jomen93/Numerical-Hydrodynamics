#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <string>

using namespace std;

// Selection of method
// Rusanov method           ----------> 1
// HLL method               ----------> 2
const int Method = 1;

// numerical difusivity
// const double beta = 0.1;
// system size
const int Nx = 128;
// number equations
const int neqs = 3;
// Courant constant
const double C = 0.5;
// adiabatic constant 
const double gam = 1.4;
// initial point x-axis
const double xo = 0;
// final point x-axis
const double xf = 1;
// space resolution x-axis
const double dx = (xf-xo)/(Nx-1);

// time variable
// time variable 
double t = 0.0;
// final time of simulation 
double tf = 1.0;
// int iteration variable
int it = 0.0;
// time step variable 
double dt;


// x position
double x[Nx+2];

// conserved variables 
double U[neqs][Nx+2];
double U_t[neqs][Nx+2];
double UP[neqs][Nx+2];

// primitive variables
double P[neqs][Nx+2];
double P_t[neqs][Nx+2];

// Flux variables
double F[neqs][Nx+2];
double F_t[neqs][Nx+2];

// type of boundary conditions

// Free boundaries  ----------------> 1
const int boundary_type = 1;


// set up system 
void Init(void)
{
	const double rho_L  = 1.0;
	const double u_L    = 0.0;
	const double P_L    = 1.0;

	const double rho_R  = 0.125;
	const double u_R    = 0.0;
	const double P_R    = 0.1;

	const double E_L    = 0.5*rho_L*u_L*u_L + P_L/(gam-1);
	const double E_R    = 0.5*rho_R*u_R*u_R + P_R/(gam-1);

	// diaphragm position
	double xd = 0.5;

	for (int i = 1; i <= Nx; ++i)
	{
		// space 
		x[i]    = xo + (i-1)*dx;
		// Conserved
		if (x[i] <= xd)
		{
			U[0][i] = rho_L;
			U[1][i] = rho_L*u_L;
			U[2][i] = E_L;
		}

		if (x[i] >= xd)
		{
			U[0][i] = rho_R;
			U[1][i] = rho_R*u_R;
			U[2][i] = E_R;
		}
	}
}

// trnsformation to primitive variables
void UtoP(double U[][Nx+2], double P[][Nx+2])
{
	for(int i = 0; i <= Nx+1; ++i)
	{
		P[0][i] = U[0][i]; 
		P[1][i] = U[1][i]/U[0][i];
		P[2][i] = (gam-1)*(U[2][i] - 0.5*pow(U[1][i], 2)/U[0][i]); 
	}
}

// Bonudary implementation 
void boundary()
{
	switch(boundary_type)
	{
		case 1:
			for (int n = 0; n < neqs; ++n)
			{
				U[n][0]      = U[n][1];
				U[n][Nx+1]   = U[n][Nx];
				U_t[n][0]    = U_t[n][1];
				U_t[n][Nx+1] = U_t[n][Nx];
			}
			break;
	}
}


void Output(int number)
{
	string primitive = "Primitive"+to_string(number)+".dat";
	char name_primitive[30];
	char iteration[14] = "iteration.dat";
	char name_method[11] = "Method.dat";
	strcpy(name_primitive, primitive.c_str());

	FILE* Primitive_dat = fopen(name_primitive, "wb");
	FILE* x_dat = fopen("x.dat", "wb");

	ofstream name_dat;
	ofstream iteration_dat;

	name_dat.open(name_method);
	iteration_dat.open(iteration);

	switch(Method)
	{
		case 1:
			name_dat << "Method" << " Lax_Method" << endl;
			break;
		case 2:
			name_dat << "Method" << " Marcormack_Method" << endl;
			break;
	}
	name_dat << "Nx " << Nx << endl;
	name_dat << "neqs " << neqs << endl;


	iteration_dat << number << endl;

	fwrite(x, sizeof(double), Nx+2, x_dat);

	UtoP(U, P);

	for (int n = 0; n < neqs; ++n)
	{
		for (int i = 1; i <= Nx; ++i)
		{
			
			fwrite(&(P[n][i]), sizeof(double),1 ,Primitive_dat);
		}
	}
	name_dat.close();
	iteration_dat.close();
}

double CFL(void)
{
	double u = 0.0;
	double umax = 1.0;

	for (int i = 1; i <= Nx; ++i)
	{
		u = fabs(P[2][i])+sqrt(gam*P[3][i]/P[0][i]);
		if (umax < u)
		{
			umax = u;
		}
	}
	return (C*dx/umax);
}


int main()
{
	Init();
	UtoP(U, P);
	printf("%f\n", CFL());
	Output(1);
	return 0;
}