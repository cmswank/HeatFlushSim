#include <boost/thread.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

class cmsRK {
private: 
	int N,timesteps,Nactual; //number of spins/time-steps simulated. Nacutal is the real spins number. 
	int startseed; //random seed. 
	int cores; //number of parallel processes to be used. 
	double T; //total time in seconds. 
	double dt; //data save time steps, (does not have to be integer time steps to T.)
	double error; //error per step;
	static constexpr double pi=3.141592653589793;
	double rkdt; //initial guess at the time step of the Runge-Kutta integrator. Should be ~5 times smaller than the collision time of he-3
	std::ofstream  spinfile;
public:
	std::vector<boost::thread *> bthreads;
	std::vector<cmsTransport*> tr;
	

	cmsRK(int N, double T, double dt, double error, double rkdt, int startseed, int cores, std::string filename);

	
	//conduct all spin runs DH = Disk Heavy, more writing to the disk, more memory friendly
	void conductorDH(int elem,std::vector<int> belement,\
							double* size,std::vector<double> hole,\
							void (*Temperature)(double*,double*),\
							void (*Convection)(double*,double*),\
							void (*Bfield)(double*,double*),\
							void (*Efield)(double*,double*,double*));

	//conduct all spins, MH = memory heavy. write to disk once at the end, memory heavy. 
	void conductorMH(bool spinsolve, int elem,std::vector<int> belement,\
							double* size,std::vector<double> hole,\
							void (*Temperature)(double*,double*),\
							void (*Convection)(double*,double*),\
							void (*Bfield)(double*,double*),\
							void (*Efield)(double*,double*,double*));
	//run a single spin
	void run(cmsTransport * tr);
	void runNS(cmsTransport * tr);
	//this is the classes work horse,
	//it is used in run to do integration and rk-error checking/step sizing. 
	void integrate(cmsTransport* tr, double t); 

	void saveSpins();
	void saveSpinsDH();
	double* Position;
	double* Velocity;
	double* Spins;
	double* BField;
	double* EField;

};