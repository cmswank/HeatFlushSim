#include <iostream>
#include <assert.h>
#include "cmsInterp.h"

double cmsInterp::cubicInterpolate (double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}
/*
double cmsInterp::bicubicInterpolate (double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}
*/
double cmsInterp::interp3D(double *pos)
{
	double tempx=((pos[0]-start[0])/step[0]);
	double tempy=((pos[1]-start[1])/step[1]);
	double tempz=((pos[2]-start[2])/step[2]);
	
	double tempp[4][4][4];
	int indx=(int)(tempx);
	int indy=(int)(tempy);
	int indz=(int)(tempz);
	int offx=0,offy=0,offz=0;
	if(indx>1){
		offx=indx-1;
		if(indx>Num[0]-3) offx=(Num[0]-4);	
		tempx=tempx-(double)offx;
	}

	if(indy>1){
		offy=indy-1;
		if(indy>Num[1]-3) offy=(Num[1]-4);	
		tempy=tempy-(double)offy;
	
	}
	
	if(indz>1){
		offz=indz-1;
		if(indz>Num[2]-3) offz=(Num[2]-4);	
		tempz=tempz-(double)offz;
	}

	tempx=tempx-1.0;
	tempy=tempy-1.0;
	tempz=tempz-1.0;
    double tempos[3]={tempx,tempy,tempz};

	for(int i = 0; i<4; i++){for(int ii = 0; ii<4; ii++){for(int iii = 0; iii<4; iii++){
	
		tempp[i][ii][iii]=data[(i+offx)*Num[1]*Num[2]+(ii+offy)*Num[2]+(offz+iii)];

	}}}

	
	return nCubicInterpolate(3,(double*) tempp,tempos);// tricubicInterpolate(tempp,tempx,tempy,tempz);

}

/*
double cmsInterp::tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
	double arr[4];
	arr[0] = bicubicInterpolate(p[0], y, z);
	arr[1] = bicubicInterpolate(p[1], y, z);
	arr[2] = bicubicInterpolate(p[2], y, z);
	arr[3] = bicubicInterpolate(p[3], y, z);
	return cubicInterpolate(arr, x);
}
*/

double cmsInterp::nCubicInterpolate (int n, double* p, double coordinates[]) {
	assert(n > 0);
	if (n == 1) {
		return cubicInterpolate(p, *coordinates);
	}
	else {
		double arr[4];
		int skip = 1 << (n - 1) * 2;
		arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
		arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
		arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
		arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
		return cubicInterpolate(arr, *coordinates);
	}
}


int main () {
	
	double dat[1000];
	int iiii=0;
	for(int i = 0; i<10; i++){
		for(int ii = 0; ii<10; ii++){
			for(int iii = 0; iii<10; iii++){
			
				dat[iiii]=(double)(i+ii+iii);
				iiii++;
			}
		}
	}

	double step[3];
	for(int i = 0; i<3; i++) step[i]=1.0;

	double start[3];
	for(int i = 0; i<3; i++) start[i]=0.0;
	int N[3]={10,10,10};
	cmsInterp* Bxi=new cmsInterp(dat,start,step,N);
	double position[3]={5.5,4.4,8.2};
	std::cout<<"interpolation gives "<<Bxi->interp3D(position)<<"\n";
	
}
