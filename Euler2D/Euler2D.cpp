#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <string>

using namespace std;

// Method definitions
const int LAX_METHOD = 1;
const int RUSANOV_METHOD = 2;
const int MACORMACK_METHOD =3;
const int HLL_METHOD = 4;
const int HLLC_METHOD = 5;
const int MUSCL_METHOD = 6;

// choose integration method
const int Method = RUSANOV_METHOD;

// numerical difusivity
const double beta = 0.1;
// system size
const int Nx = 128;
const int Ny = 128;
// number equations
const int neqs = 4;
// Courant constant
const double C = 0.5;
// adiabatic constant 
const double gam = 1.4;
// initial points 
const double xo = 0;
const double yo = 0;
// final points
const double xf = 2;
const double yf = 2;
// space resolutions
const double dx = (xf-xo)/(Nx-1);
const double dy = (yf-yo)/(Ny-1);

// time variable 
double t = 0.0;
// final time of simulation 
double tf = 2.0;
// iteration variable
int it = 0.0;
// time step variable 
double dt;
// time metrics
time_t ltime;


// factor output (this is to plot in python, please see fils plotxx.py)
const int factor = 3;

// arrays definitions

// x position
double x[Nx+2];
// y position
double y[Nx+2];

// conserved variables 
double U[neqs][Nx+2][Ny+2];
double U_t[neqs][Nx+2][Ny+2];
double UP[neqs][Nx+2][Ny+2];
double UP2[neqs][Nx+2][Ny+2];

// primitive variables
double P[neqs][Nx+2][Ny+2];
double P_t[neqs][Nx+2][Ny+2];

// Flux variables
double F[neqs][Nx+2][Ny+2];
double F_t[neqs][Nx+2][Ny+2];
double G[neqs][Nx+2][Ny+2];
double G_t[neqs][Nx+2][Ny+2];

// Auxiliary variables to MUSCL method
double deltaL;
double deltaR;
double s;

// type of boundary conditions

// right wall
const int RIGHT_FREE       = 11;
const int RIGHT_REFLECTIVE = 12;

// left wall
const int LEFT_FREE       = 21;
const int LEFT_REFLECTIVE = 22;

// top wall
const int TOP_FREE       = 31;
const int TOP_REFLECTIVE = 32;

// bottom wall
const int BOTTOM_FREE       = 41;
const int BOTTOM_REFLECTIVE = 42;

// walls selection
const int LEFT   = LEFT_REFLECTIVE;
const int RIGHT  = RIGHT_REFLECTIVE;
const int TOP    = TOP_REFLECTIVE;
const int BOTTOM = BOTTOM_REFLECTIVE;

// set up the system
void Init_explotion(void)
{
	// set-up initial conditions of physical variables
	
	const double rho_in  = 1.0;
	const double u_in    = 0.0;
	const double v_in    = 0.0;
	const double P_in    = 1.0;

	const double rho_out  = 0.125;
	const double u_out    = 0.0;
	const double v_out    = 0.0;
	const double P_out    = 0.1;

	const double E_in  = 0.5*rho_in*(pow(u_in,2)+pow(v_in,2)) + P_in/(gam-1);
	const double E_out = 0.5*rho_out*(pow(u_out,2)+pow(v_out,2)) + P_out/(gam-1);

	
	
	for (int i = 1; i <= Nx; ++i)
	{
		// space x-axis
		x[i] = xo + (i-1)*dx;
	}
	for (int j = 1; j <= Ny; ++j)
	{
		// space y-axis
		y[j] = yo + (j-1)*dy;
	}
	// definition of radius of explosion
	double r = 0.4;
	double ro;
	// center of explosion
	double xc = 0.5*xf;
	double yc = 0.5*yf;
	// Matrix filling
	for (int i = 1; i <= Nx; ++i)
	{
		for (int j = 1; j <= Ny; ++j)
		{
			ro = sqrt(pow(x[i]-xc,2)+pow(y[j]-yc,2));
			if (ro <= r)
			{
				U[0][i][j] = rho_in;
				U[1][i][j] = rho_in*u_in;
				U[2][i][j] = rho_in*v_in;
				U[3][i][j] = E_in;
			}
			else
			{
				U[0][i][j] = rho_out;
				U[1][i][j] = rho_out*u_out;
				U[2][i][j] = rho_out*v_out;
				U[3][i][j] = E_out;
			}
		}
	}
}

// Function to transform conserved variables to primitives variables
void UtoP(double U[][Nx+2][Ny+2], double P[][Nx+2][Ny+2])
{
	for(int i = 1; i <= Nx; ++i)
	{
		for(int j = 1; j <= Ny; ++j)
		{
			// density
			P[0][i][j] = U[0][i][j]; 
			// x - velocity
			P[1][i][j] = U[1][i][j]/U[0][i][j];
			// y - velocity 
			P[2][i][j] = U[2][i][j]/U[0][i][j];
			// Pressure
			P[3][i][j] = (gam-1)*(U[3][i][j] - 0.5*(pow(U[1][i][j],2)+pow(U[2][i][j],2))/U[0][i][j]); 
		}
	}
}

//  here we'll define all boundary conditions
void boundary(double U[][Nx+2][Nx+2])
{
	switch(LEFT)
	{
		case LEFT_FREE:
		{
			for (int n = 0; n < neqs; ++n)
			{
				for (int j = 0; j <= Ny+1; ++j)
				{
					U[n][0][j] = U[n][1][j];
				}
			}
		}
		break;

		case LEFT_REFLECTIVE:
		{
				for (int j = 0; j <= Ny+1; ++j)
				{
					U[0][0][j] = U[0][1][j];
					U[1][0][j] = -U[1][1][j];
					U[2][0][j] = U[2][1][j];
					U[3][0][j] = U[3][1][j];
				}
		}
		break;
	}

	switch(RIGHT)
	{
		case RIGHT_FREE:
		{
			for (int n = 0; n < neqs; ++n)
			{
				for (int j = 0; j <= Ny+1; ++j)
				{
					U[n][Nx+1][j] = U[n][Nx][j];
				}
			}
		}
		break;
		
		case RIGHT_REFLECTIVE:
		{
			for (int j = 0; j <= Ny+1; ++j)
			{
				U[0][Nx+1][j] = U[0][Nx][j];
				U[1][Nx+1][j] = -U[1][Nx][j];
				U[2][Nx+1][j] = U[2][Nx][j];
				U[3][Nx+1][j] = U[3][Nx][j];
			}
		}
		break;

	}

	switch(TOP)
	{
		case TOP_FREE:
		{
			for (int n = 0; n < neqs; ++n)
			{
				for (int i = 0; i <= Nx+1; ++i)
				{
					U[n][i][Ny+1] = U[n][i][Ny];
				}
			}
		}
		break;

		case TOP_REFLECTIVE:
		{
			for (int i = 0; i <= Nx+1; ++i)
			{
				U[0][i][Ny+1] = U[0][i][Ny];
				U[1][i][Ny+1] = U[1][i][Ny];
				U[2][i][Ny+1] = -U[2][i][Ny];
				U[3][i][Ny+1] = U[3][i][Ny];
			}
		}
		break;
	}

	switch(BOTTOM)
	{
		case BOTTOM_FREE:
		{
			for (int n = 0; n < neqs; ++n)
			{
				for (int i = 0; i <= Nx+1; ++i)
				{
					U[n][i][0] = U[n][i][1];
				}
			}
		}
		break;

		case BOTTOM_REFLECTIVE:
		{
			for (int i = 0; i <= Nx+1; ++i)
			{
				U[0][i][0] = U[0][i][1];
				U[1][i][0] = U[1][i][1];
				U[2][i][0] = -U[2][i][1];
				U[3][i][0] = U[3][i][1];
			}
		}
		break;
	}
}

void Flux(double F[][Nx+2][Ny+2], double G[][Nx+2][Nx+2], double P[][Nx+2][Ny+2])
{
	for (int i = 0; i <= Nx+1; ++i)
		{
			for (int j = 0; j <= Ny+1; ++j)
			{
				F[0][i][j] = P[0][i][j]*P[1][i][j];
				F[1][i][j] = P[0][i][j]*pow(P[1][i][j],2) + P[3][i][j];
				F[2][i][j] = P[0][i][j]*P[1][i][j]*P[2][i][j];
				F[3][i][j] = P[1][i][j]*(0.5*P[0][i][j]*(pow(P[1][i][j],2)+pow(P[2][i][j],2)) + (1/(gam-1))*P[3][i][j]);
				
				G[0][i][j] = P[0][i][j]*P[2][i][j];
				G[1][i][j] = P[0][i][j]*P[1][i][j]*P[2][i][j];
				G[2][i][j] = P[0][i][j]*pow(P[2][i][j],2) + P[3][i][j];
				G[3][i][j] = P[2][i][j]*(0.5*P[0][i][j]*(pow(P[1][i][j],2)+pow(P[2][i][j],2)) + (1/(gam-1))*P[3][i][j]);
			}
	}
}

void Integration(void)
{
	switch(Method)
	{
		// Lax - Wendroff Method implementation
		case LAX_METHOD:
		{
			boundary(U);
			UtoP(U, P);
			Flux(F, G, P);

			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					for (int n = 0; n < neqs; ++n)
					{
						// integration
						UP[n][i][j] = 0.25*(U[n][i+1][j] + U[n][i-1][j]+U[n][i][j+1] + U[n][i][j-1])
						-0.5*dt*(F[n][i+1][j]-F[n][i-1][j])/dx
						-0.5*dt*(G[n][i][j+1]-G[n][i][j-1])/dy;
					}
				}
			}

			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					for (int n = 0; n < neqs; ++n)
					{
						// stepping
						U[n][i][j]  = UP[n][i][j];
					}
				}
			}	
		}
		break;

		// MacorMack Method Implementation
		case MACORMACK_METHOD:
		{
			boundary(U);
			UtoP(U, P);
			Flux(F, G, P);

			// Predictor
			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{	
					for (int n = 0; n < neqs; ++n)
					{
						U_t[n][i][j] = U[n][i][j] - dt*(F[n][i+1][j] - F[n][i][j])/dx - dt*(G[n][i][j+1] - G[n][i][j])/dy;
					}
				}
			}

			// Boundary update
			boundary(U_t);
			// Primitive variables update
			UtoP(U_t, P_t);
			// Flux calculations
			Flux(F_t, G_t, P_t);
				
		  	// corrector
			for (int i = 1; i <= Nx; ++i)
			 {
			 	for (int j = 1; j <= Ny; ++j)
			 	{
				 	for (int n = 0; n < neqs; ++n)
				 	{
				 		UP[n][i][j] = 0.5*(U[n][i][j] + U_t[n][i][j])
				 					- 0.5*dt*(F_t[n][i][j] - F_t[n][i-1][j])/dx
				 					- 0.5*dt*(G_t[n][i][j] - G_t[n][i][j-1])/dy;
				 	}
			 	}
			 } 

			// Numerical difusivity
			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					for(int n = 0; n < neqs; ++n)
					{
						if((UP[n][i+1][j]-UP[n][i][j])*(UP[n][i][j]-UP[n][i-1][j]) < 0 || (UP[n][i][j+1]-UP[n][i][j])*(UP[n][i][j]-UP[n][i][j-1]) < 0)
						{
							U[n][i][j] = UP[n][i][j] + beta*(UP[n][i+1][j] + UP[n][i-1][j] -4*UP[n][i][j] + UP[n][i][j+1] + UP[n][i][j-1]);
						}
						else
						{
							U[n][i][j] = UP[n][i][j];
						}
					}
				}	
			}	
		}
		break;

		case RUSANOV_METHOD:
		{
			double ur;
			double ul;
			double up;
			double ud;

			double Lambda_x;
			double Lambda_y;

			UtoP(U, P);
			Flux(F, G, P);

			for (int i = 0; i <= Nx+1; ++i)
			{
				for (int j = 0; j <= Ny+1; ++j)
				{
					// sound waves x-direction
					ul = fabs(P[1][i][j])+sqrt(gam*P[2][i][j]/P[0][i][j]);
					ur = fabs(P[1][i+1][j])+sqrt(gam*P[2][i+1][j]/P[0][i+1][j]);
					// Lambda_x = max(ul, ur);
					
					if (ur >= ul)
					{
						Lambda_x = ur;
					}	
					else
					{
						Lambda_x = ul;
					}

					// sound waves in y-direction
					ud = fabs(P[1][i][j])+sqrt(gam*P[2][i][j]/P[0][i][j]);
					up = fabs(P[1][i][j+1])+sqrt(gam*P[2][i][j+1]/P[0][i][j+1]);
					// Lambda_y = max(ud, up);

					if (ud >= up)
					{
						Lambda_y = ur;
					}	
					else
					{
						Lambda_y = ul;
					}

					// intercell flux in both directions
					for (int n = 0; n < neqs; ++n)
					{
						F_t[n][i][j] = 0.5*(F[n][i][j] + F[n][i+1][j]) - 0.5*Lambda_x*(U[n][i+1][j] - U[n][i][j]);
						G_t[n][i][j] = 0.5*(G[n][i][j] + G[n][i][j+1]) - 0.5*Lambda_y*(U[n][i][j+1] - U[n][i][j]);
					}
				}
			}
			// Gonudov scheme of first order in x direction  
			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					for (int n = 0; n < neqs; ++n)
					{
						UP2[n][i][j] = U[n][i][j] - 0.5*(dt/dx)*(F_t[n][i][j] - F_t[n][i-1][j]); 
					}
				}
			}

			UtoP(UP2, P);
			Flux(F, G, P);

			// Godunov scheme of first order in y direction
			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					for (int n = 0; n < neqs; ++n)
					{
						UP[n][i][j] = UP2[n][i][j] - 0.5*(dt/dy)*(G_t[n][i][j] - G_t[n][i][j-1]); 
					}
				}
			}
			
			// Stepping
			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					for (int n = 0; n < neqs; ++n)
					{
						U[n][i][j] = UP[n][i][j];
					}
				}
			}
		boundary(U);
		}
		break;
	}
}

void report(void)
{
	cout << "=============================================================" << endl;
	switch(Method)
	{
		case LAX_METHOD:
			cout << "Method: Lax" << endl;
			break;
		case RUSANOV_METHOD:
			cout << "Method: Rusanov" << endl;
			break;
		case MACORMACK_METHOD:
			cout << "Method: Macormack" << endl;
			break;
		case HLL_METHOD:
			cout << "Method: HLL" << endl;
			break;
		case HLLC_METHOD:
			cout << "Method: HLLC" << endl;
			break;
		case MUSCL_METHOD:
			cout << "Method: MUSCL" << endl;
			break;
	}

	cout << "=============================================================" << endl;
	cout << "Iteration       " << it << "           Complete" << endl;
	ltime = time(NULL);
	cout << asctime(localtime(&ltime)) << endl;
	cout << "Time                                      " << t << endl; 
	cout << "step time                                 " << dt << endl; 
	cout << endl; 
}

double CFL(void)
{
	double u = 0.0;
	double v = 0.0;
	double umax = 1.0;
	double vmax = 1.0;

	for (int i = 1; i <= Nx; ++i)
	{
		for (int j = 1; j <= Ny; ++j)
		{
			v = fabs(P[1][i][j])+sqrt(gam*P[3][i][j]/P[0][i][j]);
			u = fabs(P[2][i][j])+sqrt(gam*P[3][i][j]/P[0][i][j]);
			if (vmax < v)
			{
				vmax = v;
			}
			if (umax < u)
			{
				umax = u;
			}
		}
	}
	// choose the min 
	if (umax >= vmax)
	{
		return (C*dx/vmax)*(1/sqrt(2));
	}
	else
	{
		return (C*dx/umax)*(1/sqrt(2));
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
		FILE* y_dat = fopen("y.dat", "wb");

		ofstream name_dat;
		ofstream iteration_dat;
		
		name_dat.open(name_method);
		iteration_dat.open(iteration);
		
		switch(Method)
		{
			case LAX_METHOD:
				name_dat << "Method Lax" << endl;
				break;
			case RUSANOV_METHOD:
				name_dat << "Method Rusanov" << endl;
				break;
			case MACORMACK_METHOD:
				name_dat << "Method Macormack" << endl;
				break;
			case HLL_METHOD:
				name_dat << "Method HLL" << endl;
				break;
			case HLLC_METHOD:
				name_dat << "Method HLLC" << endl;
				break;
			case MUSCL_METHOD:
				name_dat << "Method MUSCL" << endl;
				break;
		}

		name_dat << "Nx " << Nx << endl;
		name_dat << "Ny " << Nx << endl;
		name_dat << "neqs " << neqs << endl;
		name_dat << "it " << it << endl;
		name_dat << "factor " << factor << endl;
		
		iteration_dat << number << endl;

		fwrite(x, sizeof(double), Nx+2, x_dat);
		fwrite(y, sizeof(double), Ny+2, y_dat);

		UtoP(U, P);

		for (int n = 0; n < neqs; ++n)
		{
			for (int i = 1; i <= Nx; ++i)
			{
				for (int j = 1; j <= Ny; ++j)
				{
					fwrite(&(P[n][i][j]), sizeof(double),1 ,Primitive_dat);
				}
			}
		}
		name_dat.close();
		iteration_dat.close();	
	}
	
}

int main()
{
	ltime = time(NULL);
	printf("\n================ Starting Simulation ========================\n");
	printf("Started: %s", asctime(localtime(&ltime)));
	Init_explotion();
	boundary(U);
	Output(0);
	// Init_KH();
	while (t <= tf)
	// while (it <= 100)
	{
		dt = CFL();
		Integration();
		t  += dt;
		it += 1; 
		report();
		Output(it);			
	}
}




