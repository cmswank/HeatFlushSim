//This class is supposed to trasport a helium-3 particle in superfuid helium-4. It will be a function of temperature.
#include <assert.h>
#include <iostream>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>


class cmsTransport{
private:
		//the randomly selected amount of probability required to scatter, current prob,  etc. 
		
		
	
		static constexpr double k=1.38064852e-23; //boltzmann constant 
 		static constexpr double m=2.25*3.0160293*1.660539040e-27; //he-3 effective mass in superfuild
 		static constexpr double pi=3.141592653589793;
		
		//boost random engine. 
		typedef boost::mt19937 RNG;    // Mersenne Twister
 		typedef boost::random::uniform_real_distribution<double> DIST_flat;
		typedef boost::variate_generator<RNG, DIST_flat> RAND_flat;

		typedef boost::normal_distribution<double> DIST_velocity;   // Normal Distribution
		typedef boost::variate_generator<RNG,DIST_velocity> RAND_velocity;
		RNG rng;
	 	DIST_flat dist_flat;
	 	RAND_flat rand_flat;
	 	DIST_velocity dist_vel;
	 	RAND_velocity rand_vel;
	 	

 		void getInitial(int elem);

 		

public:
	double pscatter,lam,pscatpast,pscatnow;
	double t,rkdt,twall,pcurrent;
	int nnn;
	double* tempT;
	double tempTime;
	double temptwall;
	double temptime;
	double temppcurrent;
	double tempppast;
	int tempelement;
	int tempElement;

	 	//cmsBoundary class. 
 		cmsBoundary* bd;

	//pointer to temperature function. 
	void (*Temperature)(double*,double*);

 	//pointer to B field function (can be anything, or use the spline interpolation with a dataset in class cmsInterp.)
 	void (*Bfield)(double*,double*);
 	//pointer to E field function same as B field function but for E field. 
 	void (*Efield)(double*,double*,double*);
	//pointer to convection velocity field function, maybe use the spline, maybe not.
	void (*Convection)(double*,double*);

	cmsTransport(int flat_seed,int elem,std::vector<int> belement,\
							double* size,std::vector<double> hole,\
							void (*Temperature)(double*,double*),\
							void (*Convection)(double*,double*),\
							void (*Bfield)(double*,double*),\
							void (*Efield)(double*,double*,double*));
	void propagate(double dt);

	//get the Kutta from Runge
	void solve();

	//initializing memory
	//position will be initialized from getInitial in boundary.  
	double* position;
	double* tempSpin;
	double* velocity;
	double* spin;
	double* bfield;
	double* efield;
	double* bmotional;
	double* btotal;
	double* tempU;
	//state saving required for error checking
	double* cash;
	double* karp;
	double* temppos;
	double* tempvel;
	double* tempbouncepos;
	int tempwhatwall;


	//state saving required in the RK integration
	double* tempPos;
	double* tempVel;
	double* tempBouncePos;
	double tempWallt;
	double tempScatp;
	int tempWhatWall;
	double tempScatold;


	//RK constants. 
	double* k1;
	double* k2;
	double* k3;
	double* k4;
	double* k5;
	double* k6;


	double* Position;
	double* Velocity;
	double* Spins;
	double* BField;
	double* EField;
static constexpr double gamma=-203794709.3;
	//Runge-Kutta Cash Karp.
  static constexpr double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;
  static constexpr double b21 = 0.2;
  static constexpr double b31 = 0.075, b32 = 0.225;
  static constexpr double b41 = 0.3, b42 = -0.9, b43 = 1.2;
  static constexpr double b51 = -2.03703703703703692e-1, b52 = 2.5, b53 = -2.59259259259259256, b54 = 1.29629629629629628;
  static constexpr double b61 = 2.94958043981481469e-2,b62 = 3.41796875e-1,b63 = 4.15943287037037063e-2,b64 = 4.00345413773148140e-1,b65 = 6.1767578125e-2;
  static constexpr double c1 = 9.78835978835978782e-2, c2 = 0.,c3 = 4.02576489533011284e-1,c4 = 2.10437710437710451e-1,c5 = 0.,c6 = 2.89102202145680387e-1;
  static constexpr double d1 = 1.02177372685185189e-1,d2 = 0.,d3 = 3.83907903439153431e-1,d4 = 2.44592737268518517e-1,d5 = 1.93219866071428562e-2,d6 = 0.25;

};