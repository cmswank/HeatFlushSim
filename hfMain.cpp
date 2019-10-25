//g++ hfMain.cpp cmsRK.cpp cmsTransport.cpp cmsBoundary.cpp -I/data1/cmswank/BoostCodeSwank/boost_1_64_0 -o Trajtest -lboost_thread -lboost_system
#include "cmsBoundary.h"
#include "cmsTransport.h"
#include "cmsRK.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <boost/thread.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>



//Temperature function... put in what you want!
void Temperature(double* temperature,double* pos ){
		temperature[0]=0.45;//0.45-0.01*pos[2];
		return;
}

//Bfield function
void Bfield(double* bfield,double* pos){
	bfield[0]=-1.5e-7*pos[0];
	bfield[1]=0.0;
	bfield[2]=3e-6+1.5e-7*(pos[2]-0.5);

	return;  //3e-6 T is 30 mG in Tesla 3.3E-10 is 3ppm/cm
}

//Efield function
void Efield(double* efield,double* pos, double* vel){
	efield[0]=0.0;
	efield[1]=0.0;
	efield[2]=0.0;
				 //changed to small value temporarily for debuggin puroposes.

	return; //75 kV/cm in kV/m 75E7
}

//Convection function
void Convection(double* conv,double* pos){
	conv[0]=0.0;
	conv[1]=0.0;
	// 80 cm/s max speed in a r^2 profile. until it gets to the cell...
	conv[2]=pos[2]>1.0 ? 0.1 : 0.8*(1.0-(pos[0]*pos[0]+pos[1]*pos[1])/4.1209e-04);
	return; 
}


int main()
{
	//boundary element (cylinder or rectangle)
	std::vector<int> belem;
	belem.push_back(1);  //cylinder
	belem.push_back(1); //rectangle (thats all there is (coded) folks). 
	
	// exit hole radius in meters
	std::vector<double> hole;
	hole.push_back(0.45*1.6*2.54/100.0/2.0); //this must not let them out, I did not code them leaving!!!
     hole.push_back(-1.0); // a negative number will certianly never exit the hole. 
	
	//					belem=0 		 belem=1
	//size of object (R, Lz ,Lz) or (Lx, Ly, Lz)
	double * tempsize=(double*) new double[2][3];
	tempsize[0]=1.6*2.54/100.0/2.0;
	tempsize[1]=1.6*2.54/100.0/2.0;
	tempsize[2]=1.0;
	tempsize[3]=.4;  //nEDM@SNS cell size.  
	tempsize[4]=0.102;
	tempsize[5]=0.076;
	std::string filename="/data1/cmswank/SpinSimSwank2/data/spindataTest.dat";
	//Runge-Kutta integrator construction
	//double num;
//	for (int i =0; i<100;i++)
//	{
//	std::cout<<"what is the double?\n";
//	std::cin>>num;
//	std::cout<<sgn(num)<<"\n";
//	}


	   				  //Number   Time    dt     error    rkdt       seed      cores
	cmsRK *rk=new cmsRK(100,    50.0,   0.2,    1E-9,   2.9001E-3 ,    540,     20,   filename);
	
	//this 'conducts' the spins (and saves a file)
	//disk heavy, memory light (not that heavy and not that light, not worth writing actually...)
	//don't use this one it doesn't seem to save much memory adds a bunch of waiting. 		
	//rk->conductorDH(0,belem,(double*)tempsize,hole,Temperature,Convection,Bfield,Efield);	

	//this also 'conducts' the spins (and saves a file)
	//memory heavy, disk light
	//first entry is if you want to calculate spins or just trajectories. 
	//true will calculate spins.(bfield efield spins)
				    //sove for spins? starting element, 				
	rk->conductorMH(true,0,belem,tempsize,hole,Temperature,Convection,Bfield,Efield);
	
	//discription of input variables. 
	//conductorMH(solve for spins?, starting element number, element tags, element sizes,  hole sizes, Temperature function... etc)
 	
	



	//Transport Class testing... (pretty good for boundary as well.)
	/*
	cmsTransport * tr= new cmsTransport( 82+96,0,belem,(double*)tempsize,hole,Temperature,Convection,Bfield,Efield);
	


	for(int i = 0 ; i < 1E7; ++i)
    {	

    	tr->propagate(1.1e-4);
 		if (i%(int)1E6 ==0)
 			{ 
 				std::cout<<"hello we made it to "<<i/((int)1E6)<<" million time steps\n";
 				std::cout<<"position: "<<tr->position[0]<<" "<<tr->position[1]<<" "<<tr->position[2]<<"\n";
 				std::cout<<"veloccity: "<<tr->velocity[0]<<" "<<tr->velocity[1]<<" "<<tr->velocity[2]<<"\n";
 			}

    }
    
    */

return 0;

}