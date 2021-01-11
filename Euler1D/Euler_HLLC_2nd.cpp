#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <cmath>

using namespace std;

// Method definitions
const int HLLC_METHOD = 1;

// MEthod Selection
const int Method = HLLC_METHOD;

// system size 
const int Nx = 128;
//  number of equations 
const int neqs = 3;
// Courant constant
const double C = 0.5;
// adiabatic constant 
const double gam = 1.4;
// initial  point x-axis
const double xo = 0;
// final point x-axis
const double xf = 1;
// space resolution x-axis
const double dx = (xf - xo)/(Nx - 1);

// time variables
double t;
// final time of simulation
double tf = 0.5;
// time step variable
double dt;
// iterator variable
int it;

// factor output to print less files
const int factor = 1;

// x position
double x[Nx+4];
// conserverd variables
double U[neqs][Nx+4];
double UP[neqs][Nx+4];
// Prmimitivs variables
double P[neqs][Nx+4];
// Flux variables
double F[neqs][Nx+4];
double F_t[neqs][Nx+4];
double FF[neqs];
double Fl[neqs];
double Fr[neqs];



// set up system 
void Initial_conditions()
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

	for (int i = 2; i <= Nx+1; ++i)
	{
		// space 
		x[i] = xo + (i-1)*dx;
		// Conserved
		if (x[i] <= xd)
		{
			U[0][i] = rho_L;
			U[1][i] = rho_L*u_L;
			U[2][i] = E_L;
		}
		else
		{
			U[0][i] = rho_R;
			U[1][i] = rho_R*u_R;
			U[2][i] = E_R;
		}
	}
	// reset the counters an time to zero
	t  = 0;
	it = 0;
}

// computes the primitives as a function of U vector
void UtoP(double U[][Nx+4], double P[][Nx+4])
{
	for(int i = 0; i <= Nx+3; ++i)
	{
		P[0][i] = U[0][i]; 
		P[1][i] = U[1][i]/U[0][i];
		P[2][i] = (gam-1)*(U[2][i] - 0.5*pow(U[1][i], 2)/U[0][i]); 
	}
}

// Computes the primitives as a function of the U's, only one cell
void PtoU(double pp[], double uu[])
{
	uu[0] = pp[0];
	uu[1] = pp[0]*pp[1];
	uu[2] = 0.5*pp[0]*pp[1]*pp[1] + pp[2]/(gam - 1.0);
}

void FluxCell(double pp[], double ff[])
{
	ff[0] = pp[0]*pp[1];
	ff[1] = pp[0]*pp[1]*pp[1] + pp[2];
	ff[2] = pp[1]*( 0.5*pp[0]*pp[1]*pp[1] + gam*pp[2]/(gam-1.0) );
}

void Output(int number)
{
	if (number % factor == 0)
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
			case HLLC_METHOD:
				name_dat << "Method" << " HLLC " << endl;
				break;
		}

		name_dat << "Nx " << Nx << endl;
		name_dat << "neqs " << neqs << endl;
		name_dat << "it " << it << endl;
		name_dat << "factor " << factor << endl;

		iteration_dat << number << endl;

		fwrite(x, sizeof(double), Nx+4, x_dat);

		UtoP(U, P);

		for (int n = 0; n < neqs; ++n)
		{
			for (int i = 2; i <= Nx+1; ++i)
			{
				
				fwrite(&(P[n][i]), sizeof(double),1 ,Primitive_dat);
			}
		}
		name_dat.close();
		iteration_dat.close();
	}
}

double time_step(void)
{
	double cs;
	double del = 1e30;

	for (int i = 0; i <= Nx+3; ++i)
	{
		// cs = fabs(P[1][i])+sqrt(gam*P[2][i]/P[0][i]);
		cs = sqrt(gam*P[2][i]/P[0][i]);
		del = min(del, dx/(fabs(P[1][i]) + cs));
	}
	return (C*del);
}


double average(double a,double b)
{
	double s;
	s = copysign(1.0, a);
	return s*max(0.0, min(fabs(a), s*b));
}

// Obtain the HLLC Fluxes
void P2HLLC(double Pl[], double Pr[])
{
	double Csl;
	double Csr;

	double Sl;
	double Sr;
	double Ss;

	double Ek;

	double ust[neqs];
	
	// double FF[neqs];
	double uu[neqs];

	Csl = sqrt(gam*Pl[2]/Pl[0]);
	Csr = sqrt(gam*Pr[2]/Pr[0]);

	Sl = min(Pl[1] - Csl, Pr[1] - Csr);
	Sr = max(Pl[1] + Csl, Pr[1] + Csr);

	Ss = (Pr[2] - Pl[2] + Pl[0]*Pl[1]*(Sl-Pl[1]) - Pr[0]*Pr[1]*(Sr-Pr[1]))/(Pl[0]*(Sl - Pl[1]) - Pr[0]*(Sr - Pr[1]));

	if ( Sl > 0.0 )
	{
		FluxCell(Pl, FF);
	}
	else if (Ss >= 0.0)
	{
		ust[0] = Pl[0]*(Sl - Pl[1])/(Sl-Ss);
		ust[1] = ust[0]*Ss;

		Ek     = (0.5*Pl[0]*Pl[1]*Pl[1]) + Pl[2]/(gam - 1.0);

		ust[2] = ust[0]*((Ek/Pl[0] + (Ss - Pl[1]))*(Ss + Pl[2]/(Pl[0]*(Sl - Pl[1]))));

		FluxCell(Pl, Fl);
		PtoU(Pl, uu);

		for (int n = 0; n < neqs; ++n)
		{
			FF[n] = Fl[n] + Sl * (ust[n] - uu[n]);
		}
	}
	else if (Ss <= 0.0)
	{
		ust[0] = Pr[0]*(Sr - Pr[1])/(Sr-Ss);
		ust[1] = ust[0]*Ss;

		Ek     = (0.5*Pr[0]*Pr[1]*Pr[1]) + Pr[2]/(gam - 1.0);

		ust[2] = ust[0]*((Ek/Pr[0] + (Ss - Pr[1]))*(Ss + Pr[2]/(Pr[0]*(Sr - Pr[1]))));

		FluxCell(Pr, Fr);
		PtoU(Pr, uu);

		for (int n = 0; n < neqs; ++n)
		{
			FF[n] = Fr[n] + Sr * (ust[n] - uu[n]);
		}
	}
	else if (Sr < 0.0)
	{
		FluxCell(Pr, FF);
	}
}


// Computes the HLLC fluxes in the entires domain 
void HLLC_Flux(double P[][Nx+4],double F[][Nx+4], int order)
{
	double Pl[neqs];
	double Pr[neqs];

	double Pll[neqs];
	double Prr[neqs];

	double dl;
	double dm;
	double dr;
	double al;
	double ar;

	double FF[neqs];

	for (int n = 0; n < neqs; ++n)
	{
		for (int i = 2; i <= Nx+1; ++i)
		{
			Pl[n] = P[n][i];
			Pr[n] = P[n][i+1];

			if (order == 2)
			{
				Pll[n] = P[n][i-1];
				Prr[n] = P[n][i+2];

				dl = Pl[n]-Pll[n];
				dm = Pl[n]-Pll[n];
				dr = Prr[n]-Pr[n];
				
				al = average(dl, dm);
				ar = average(dm, dr);

				Pl[n] = Pl[n] + 0.5*al;
				Pr[n] = Pr[n] + 0.5*ar;
			}
			else
			{
				P2HLLC(Pl, Pr);
				F[n][i] = FF[n];			
			}
			
		}
	}
}

void Boundaries(double U[][Nx+4])
{
	for (int n = 0; n < neqs; ++n)
	{
		U[n][0] = U[n][2];
		U[n][1] = U[n][2];
		U[n][Nx+2] = U[n][Nx+1];
		U[n][Nx+3] = U[n][Nx+1]; 
	}
}

void integration (double dt, double t)
{
	// Obtain the fluxes 
	HLLC_Flux(P, F, 1);

	// first step
	for (int n = 0; n < neqs; ++n)
	{
		for (int i = 2; i <= Nx+1; ++i)
		{
			UP[n][i] = U[n][i] - 0.5*dt*(F[n][i] - F[n][i-1])/dx;
		}
	}

	// Boundary conditions to the U^n+1
	Boundaries(UP);

	// second step
	// updates primitives
	UtoP(U, P);
	// Obtain Fluxes
	HLLC_Flux(P, F, 1);

	// second step
	for (int n = 0; n < neqs; ++n)
	{
		for (int i = 2; i <= Nx+1; ++i)
		{
			UP[n][i] = U[n][i] - 0.5*dt*(F[n][i] - F[n][i-1])/dx;
		}
	}

	// Boundary conditions to the U^n+1
	Boundaries(UP);

    // copy the UP to the U
    for (int n = 0; n < neqs; ++n)
     {
     	for (int i = 2; i <= Nx+1; ++i)
     	{
     		U[n][i] = UP[n][i];
     	}
    } 
}

int main()
{
	// generate the initial conditions
	Initial_conditions();
	// main loop, iterate until maximun time is reached 
	while (it <= 100)
	{
		// updates the primitives 
		UtoP(U, P);
		// Obtain time step by CFL criterium
		dt = time_step();
		// Integrate U from t ---> t + dt
		integration(dt, t);
		// time increment
		t += dt;
		// iterator increment
		it += 1;
		printf("%f\n", t);
		Output(it);
	}
}