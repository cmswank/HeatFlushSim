#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "sgn.h"
//#include <boost/thread.hpp>
//#include <boost/math/special_functions.hpp>
//#include <boost/random.hpp>
//#include <boost/random/uniform_real_distribution.hpp>





//easy testing program. 
int main()
{
	//does this work for doubles like it does for classes??
	double* poop=(double*)new double[3];

	poop[0]=1.1234; 

	double* turd= new double(*poop);
	std::cout<<poop[0]<<"\n";
	std::cout<<turd[0]<<"\n";

	turd[0]=2.345;
	
	std::cout<<poop[0]<<"\n";
	std::cout<<turd[0]<<"\n";

	poop[0]=turd[0];


	turd[0]=3.456;

	std::cout<<poop[0]<<"\n";
	std::cout<<turd[0]<<"\n";


	float * answer=(float*) new float[1];
	float * x=(float*) new float[1]; 
	std::cout<<"hello\n";
	x[0]=1.0f;

	float xx = 1.0;
	//FastInvSqrt(answer, x);
	float ans2=FISqrt(x[0]);
	std::cout<<answer[0]<<"\n";
	std::cout<<ans2<<"\n";
	return 0;




}	