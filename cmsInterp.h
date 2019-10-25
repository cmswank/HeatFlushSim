#include <iostream>
#include <assert.h>

class cmsInterp {
private:
	double *data;
	double *step;
	double *start;
	int dsize;
	int *Num;

	double cubicInterpolate (double p[4], double x);
	double nCubicInterpolate (int n, double* p, double coordinates[]);
	//double bicubicInterpolate (double p[4][4], double x, double y);
	//double tricubicInterpolate (double p[4][4][4], double x, double y, double z);
	

public:
	//this defines the data in flattened C type array, (last column moves fastest) and the step size for all arrays. 
	cmsInterp(double *dat, double *star,double *ste, int N[3]){
		Num=N;
		dsize=N[0]*N[1]*N[2];
		data=dat;
		start=star;
		step=ste;
	};


	
	//this is the function you want to interp. must be from a rectangular grid. 
	//must be an interpolation (no checking done to increase speed, it will extrapolate)
	//this is a tricubic spline interpolation. 
	//if its too slow I guess it could be made into a linear with ease. 
	double interp3D(double* pos);




};