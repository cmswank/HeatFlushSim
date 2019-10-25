
#include "sgn.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <boost/thread.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>



class cmsBoundary{
private:
		/***
		OK Folks, the first boundary in z always starts at zero. the last element always ends an L max.  
	    boundaries and holes are always centered in the xy plane. 
		This is meant to RK5 integrate in a heat flush. Heatflush and Fields must come from somewhere else
		I have implemented a cubic spline that can be initialized with a dataset of rectangular grid. 
		OR you can use your own functions and assign the pointers in main(),
		feel free to add tricky stuff.
		***/
	static constexpr double pi=3.141592653589793;
	
	
	
	double dtwall; //time requred from current position to get to the wall. 
	
	double* tempv;
	

public:
		double* sizes; //size of thingy. if its a box, Lx,Ly,Lz, if its a cylinder R,Lz, whatever (not nothing), 
		std::vector<double> holes;     //circular hole radius, RADIUS MUST BE AS SMALL AS THE SMALLEST X and Y for that volume. I feel like that should be obvious, but who knows with you people.  
		std::vector<int> direction; //0 is z , 1 is x, 2 is y (verticle)
		int current_element;
		double* bounce_pos;
		int whatwall;
		bool benter; //did we enter a new volume yes=true. 
		bool bleave;
		std::vector<int> element;       //3D cylinder 0, 3D rectangle 1. Thats all I'm doing. so foff.  
	
	 
		std::vector<double> bstart;  //starting point of elements in z. 

		cmsBoundary(std::vector<int> belement,double *bsize, std::vector<double> hole);
		//randomly select in first volume.
		double* getInitPos(int elem,double randx, double randy, double randz);

		double nextBounce(double* position,double* velocity);
		double getBounce(double rand0, double rand1,double speed,double *position, double *velocity);

};