#ifndef __HeatFlushFields_H__
#define __HeatFlushFields_H__
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "BFieldInterp.h"
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
//#include <stdio.h>
#include <string>

inline bool fileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


class HeatFlushFields {
private:

public:

	bool bxon;
	bool byon;
	bool bzon;
	bool eon;
	bool uxon;
	bool uzon;
	bool uyon;
	bool ton;

	BFieldInterp* Bx;
	BFieldInterp* By;
	BFieldInterp* Bz;
	BFieldInterp* Efield;
	BFieldInterp* convectx;
	BFieldInterp* convecty;
	BFieldInterp* convectz;
	BFieldInterp* temperature;

	HeatFlushFields(double start[3],double step[3],int N[3],\
				     std::string filenamex, std::string filenamey, std::string filenamez,\
					 std::string filenamee,\
					 std::string filenamecx, std::string filenamecy,std::string filenamecz,\
					 std::string filenamet);
	



};
#endif