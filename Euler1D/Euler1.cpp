#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <string>

using namespace std;

// Method definitions
const int LAX_METHOD = 1;
const int RUSANOV_METHOD = 2;
const int HLL_METHOD = 3;
const int HLLC_METHOD = 4;
const int MUSCL_METHOD = 5;

// Method selection
const int Method = MUSCL_METHOD;

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
double tf = 0.5;
// int iteration variable
int it = 0.0;
// time step variable 
double dt;

// factor output
const int factor = 3;

// x position
double x[Nx+2];

// definitions or the variables

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

// Auxiliary variables to MUSCL method
double deltaL;
double deltaR;
double s;
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
	}

	for (int i = 1; i <= Nx; ++i)
	{
		
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
void boundary(double U[][Nx+2])
{
	switch(boundary_type)
	{
		case 1:
			for (int n = 0; n < neqs; ++n)
			{
				U[n][0]      = U[n][1];
				U[n][Nx+1]   = U[n][Nx];
			}
			break;
	}
}

void Flux(double F[][Nx+2], double P[][Nx+2])
{
	for (int i = 0; i <= Nx+1; ++i)
	{
		F[0][i] = P[0][i]*P[1][i];
		F[1][i] = P[0][i]*pow(P[1][i],2) + P[2][i];
		F[2][i] = P[1][i]*(0.5*P[0][i]*pow(P[1][i],2) + (gam/(gam-1))*P[2][i]);	
	}
}

void Rusanov()
{
	double ur;
	double ul;

	double Lambda;

	for (int i = 0; i <= Nx; ++i)
	{
		ul = fabs(P[1][i])+sqrt(gam*P[2][i]/P[0][i]);
		ur = fabs(P[1][i+1])+sqrt(gam*P[2][i+1]/P[0][i+1]);
		// printf("%f %f \n",ul, ur);
		if (ur >= ul) 
			Lambda = ur;	
		else 
			Lambda = ul;

		for (int n = 0; n < neqs; ++n)
		{	
			// F_t[n][i] = 0.5*(F[n][i] - F[n][i+1]);
			F_t[n][i] = 0.5*(F[n][i] + F[n][i+1]) - 0.5*Lambda*(U[n][i+1] - U[n][i]);		
		}
	}
}

void Godunov()
{
	for (int i = 1; i <= Nx; ++i)
	{
	 	for (int n = 0; n < neqs; ++n)
	 	{
	 		UP[n][i] = U[n][i] - (dt/dx)*(F_t[n][i] - F_t[n][i-1]);
	 	}
	} 

	for (int i = 0; i <= Nx+1; ++i)
	{
		for (int n = 0; n < neqs; ++n)
		{
			// stepping
			U[n][i]  = UP[n][i];
		}
	}
}

void Godunov2()
{
	for (int i = 1; i <= Nx; ++i)
	{
	 	for (int n = 0; n < neqs; ++n)
	 	{
	 		UP[n][i] = U[n][i] - 0.5*(dt/dx)*(F_t[n][i] - F_t[n][i-1]);
	 	}
	} 

	for (int i = 0; i <= Nx+1; ++i)
	{
		for (int n = 0; n < neqs; ++n)
		{
			// stepping
			U[n][i]  = UP[n][i];
		}
	}
}

void HLL()
{
	double Sl;
	double Sr;

	double ul;
	double ur;

	double Csl;
	double Csr;

	for (int i = 0; i <= Nx; ++i)
	{
		Csl = sqrt(gam*P[2][i]/P[0][i]);
		Csr = sqrt(gam*P[2][i+1]/P[0][i+1]);
	
		ul = P[1][i];
		ur = P[1][i+1];

		Sl = min(ul - Csl, ur - Csr);
		Sr = max(ul + Csl, ur + Csr);

		for (int n = 0; n < neqs; ++n)
		{
			if( Sl > 0 )
			{
				F_t[n][i] = F[n][i];
			}
			else if( Sl <= 0 || 0 <= Sr )
			{
				F_t[n][i] = (Sr*F[n][i] - Sl*F[n][i+1] + Sl*Sr*(U[n][i+1] - U[n][i]))/(Sr - Sl);
			}
			else
			{
				F_t[n][i] = F[n][i+1];
			}
		}

	}
}

void HLLC()
{
	double Sl;
	double Sr;
	double Ss;

	double rhol;
	double rhor;

	double ul;
	double ur;
	double uol;
	double uor;
	double Ur;
	double Ul;
	
	double Pl;
	double Pr;

	double Csl;
	double Csr;

	double El;
	double Er;

	for (int i = 0; i <= Nx; ++i)
	{
		Csl = sqrt(gam*P[2][i]/P[0][i]);
		Csr = sqrt(gam*P[2][i+1]/P[0][i+1]);
	
		rhol = P[0][i];
		rhor = P[0][i+1];
		
		ul = P[1][i];
		ur = P[1][i+1];
		
		Pl = P[2][i];
		Pr = P[2][i+1];

		Sl = min(ul - Csl, ur - Csr);
		Sr = max(ul + Csl, ur + Csr);
		Ss = (Pr-Pl + rhol*ul*(Sl-ul) - rhor*ur*(Sr-ur))/(rhol*(Sl-ul) - rhor*(Sr-ur));

		uol = rhol*(Sl - ul)/(Sl - Ss);
		uor = rhor*(Sr - ur)/(Sr - Ss);

		El = 0.5*rhol*pow(ul,2) + Pl/(gam-1);
		Er = 0.5*rhor*pow(ur,2) + Pr/(gam-1);

		Ul = uol*((El/rhol) + (Ss - ul)*(Ss + Pl/(rhol*(Sl - ul))));
		Ur = uor*((Er/rhor) + (Ss - ur)*(Ss + Pr/(rhor*(Sr - ur))));

		if( 0 <= Sl )
		{	
			for (int n = 0; n < neqs; ++n)
			{
				F_t[n][i] = F[n][i];
			}
		}
		else if( Sl <= 0 || 0 <= Ss )
		{	
			F_t[0][i] = F[0][i] + Sl*(uol - U[0][i]);
			F_t[1][i] = F[1][i] + Sl*(uol*Ss - U[1][i]);
			F_t[2][i] = F[2][i] + Sl*(Ul - U[2][i]);
		}
		else if( Ss <= 0 || 0 <= Sr )
		{ 
			F_t[0][i] = F[0][i+1] + Sr*(uor - U[0][i+1]);
			F_t[1][i] = F[1][i+1] + Sr*(uor*Ss - U[1][i+1]);
			F_t[2][i] = F[2][i+1] + Sr*(Ur - U[2][i+1]);
		}
		else
		{
			for (int n = 0; n < neqs; ++n)
			{
				F_t[n][i] = F[n][i+1];
			}
		}
	}
}


void MUSCL()
{
	double m; 
	double phi;
	for (int i = 1; i <=Nx+1 ; ++i)
	{
		for (int n = 0; n < neqs; ++n)
		{
			deltaL = P[n][i] - P[n][i-1];
			deltaR = P[n][i+1] - P[n][i];

			s = copysign(1.0, deltaL);
			m = min(fabs(deltaL),s * deltaR);
			phi = max(0.0, m);
			
			P[n][i] = P[n][i] + 0.5*phi;

			deltaL = P[n][i+1] - P[n][i];
			deltaR = P[n][i+2] - P[n][i+1];			

			s = copysign(1.0, deltaL);
			m = min(fabs(deltaL),s * deltaR);
			phi = max(0.0, m);

			P[n][i+1] = P[n][i+1] - 0.5*phi;
		}
	}
}

void Lax()
{
	Flux(F, P);

	for (int i = 1; i <= Nx; ++i)
	{
		for (int n = 0; n < neqs; ++n)
		{
			// integration
			UP[n][i] = 0.5*((U[n][i+1] + U[n][i-1]) - dt*(F[n][i+1]-F[n][i-1])/dx);
		}
	}

	for (int i = 0; i <= Nx+1; ++i)
	{
		for (int n = 0; n < neqs; ++n)
		{
			// stepping
			U[n][i]  = UP[n][i];
		}
	}
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
			case LAX_METHOD:
				name_dat << "Method" << " Lax " << endl;
				break;
			case RUSANOV_METHOD:
				name_dat << "Method" << " Rusanov " << endl;
				break;
			case HLL_METHOD:
				name_dat << "Method" << " HLL " << endl;
				break;
			case MUSCL_METHOD:
				name_dat << "Method" << " MUSCL " << endl;
				break;
			case HLLC_METHOD:
				name_dat << "Method" << " HLLC " << endl;
				break;
		}

		name_dat << "Nx " << Nx << endl;
		name_dat << "neqs " << neqs << endl;
		name_dat << "it " << it << endl;
		name_dat << "factor " << factor << endl;


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
		}
		

double CFL(void)
{
	double u = 0.0;
	double umax = 1.0;

	for (int i = 1; i <= Nx; ++i)
	{
		u = fabs(P[1][i])+sqrt(gam*P[2][i]/P[0][i]);
		if (umax < u)
		{
			umax = u;
		}
	}
	return (C*dx/umax);
}

void report(void)
{
	cout << "=============================================================" << endl;
	switch(Method)
	{
		case LAX_METHOD:
			cout << "Integration Method: Lax" << endl;
			cout << "Iteration = " << it << endl;
			cout << "Time = " << t << endl;
			break;
		case RUSANOV_METHOD:
			cout << "Integration Method: Rusanov" << endl;
			cout << "Iteration = " << it << endl;
			cout << "Time = " << t << endl;
			break;
		case HLL_METHOD:
			cout << "Integration Method: HLL" << endl;
			cout << "Iteration = " << it << endl;
			cout << "Time = " << t << endl;
			break;
		case MUSCL_METHOD:
			cout << "Integration Method: MUSCL" << endl;
			cout << "Iteration = " << it << endl;
			cout << "Time = " << t << endl;
			break;
		case HLLC_METHOD:
			cout << "Integration Method: HLLC" << endl;
			cout << "Iteration = " << it << endl;
			cout << "Time = " << t << endl;
			break;
	}
}


int main()
{
	Init();
	boundary(U);
	Output(0);
	while(t <= tf)
	// while(it <= 100)
	{
		switch(Method)
		{
			case LAX_METHOD: 
				dt = CFL();
				boundary(U);
				UtoP(U, P);
				Lax();
				t  += dt;
				it += 1;
				break;

			case RUSANOV_METHOD: 
				UtoP(U, P);
				dt = CFL();
				Flux(F, P);
				Rusanov();
				Godunov();
				printf("%f\n", dt);
				t  += dt;
				it += 1;
				boundary(U);
				break;

			case HLL_METHOD: 
				UtoP(U, P);
				dt = CFL();
				Flux(F, P);
				HLL();
				Godunov();
				printf("%f\n", dt);
				t  += dt;
				it += 1;
				boundary(U);
				break;

			case MUSCL_METHOD:
				UtoP(U, P);
				dt = CFL();
				printf("%f\n", dt);
				t  += dt;
				it += 1;
				Flux(F, P);
				HLL();
				Godunov2();
				boundary(U);
				UtoP(U, P);
				MUSCL();
				HLL();
				Godunov();
				boundary(U);
				break;

			case HLLC_METHOD:
				UtoP(U, P);
				dt = CFL();
				Flux(F, P);
				HLLC();
				Godunov();
				printf("%f\n", dt);
				t  += dt;
				it += 1;
				boundary(U);
				break;

			default:
				printf("choose an integration method \n");
		}
		Output(it);			
		report();
	}	
}