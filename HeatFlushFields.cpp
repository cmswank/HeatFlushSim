#include "HeatFlushFields.h"
#include "BFieldInterp.h"


HeatFlushFields::HeatFlushFields(double start[3],double step[3],int N[3],std::string filenamex, std::string filenamey, std::string filenamez,\
						 std::string filenamee, std::string filenamecx, std::string filenamecy, std::string filenamecz, std::string filenamet){

//std::cout<<"howdy"<<std::endl;	
						 
double *star=start;
double *ste=step;
/*
for(int i =0; i<3; i++){
	star[i]=start[i];
	ste[i]=step[i]; 
}*/

							 
bxon=fileExists(filenamex);
byon=fileExists(filenamey);
bzon=fileExists(filenamez);
eon=fileExists(filenamee);
uxon=fileExists(filenamecx);
uyon=fileExists(filenamecy);
uzon=fileExists(filenamecz);
ton=fileExists(filenamet);

if(bxon)
	Bx= new BFieldInterp(filenamex,star,ste,N);
if(byon)
	By= new BFieldInterp(filenamey,star,ste,N);
if(bzon)
	Bz= new BFieldInterp(filenamez,star,ste,N);
if(eon)
	Efield= new BFieldInterp(filenamee,star,ste,N);
if(uxon)
	convectx= new BFieldInterp(filenamecx,star,ste,N);
if(uyon)
	convecty= new BFieldInterp(filenamecy,star,ste,N);
if(uzon)
	convectz= new BFieldInterp(filenamecz,star,ste,N);
if(ton)
	temperature= new BFieldInterp(filenamet,star,ste,N);
	
	
	
}