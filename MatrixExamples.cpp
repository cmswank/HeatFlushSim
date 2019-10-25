#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec.hpp>
#include <boost/qvm/quat.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/mat_traits_array.hpp>
#include <boost/qvm/vec_traits_array.hpp>
//#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/map_vec_mat.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/mat_access.hpp>
#include <fftw3.h>
#include "T2Decay.h"
#include "cmsBField.h"
#include <complex.h>
#include <cerf.h>


//testing out stuff. 
int main()
{	
	//double complex test= w_of_z(1.0+I*2.0);
	//cout<<test<<"\n";
    //using namespace boost::lambda;
	using namespace boost::qvm;
    quat<float> rx=rotx_quat(3.14159f);
	vec<float,3> v={0,0,7};
	///boost can invert a matrix fast as hell. 
	//might be usefull for minimization procedures 
	mat<float,4,4> tr=translation_mat(v);
	mat<float,4,4> itr=inverse(tr);
	//float v[3]={0,1,2};
	float out1=A<1,1>(itr);
	float out15=A<1,2>(itr);
	float out2=A<2,2>(itr);
	std::cout<<out1<<"\n";
	std::cout<<out15<<"\n";
	std::cout<<out2<<"\n";
	// END OF BOOST MATRIX FUN STUFF. 

	//T2 Decay / FID spectrum stuff. 
	int N=1000;
	T2Decay* signal=new T2Decay(1000.0, 100, 1, .001,N);	  

	std::cout<<signal->SFID[22][0]<<"\n";
	std::cout<<signal->SFID[22][1]<<"\n";    
	
	signal->update(500.0,97.4,10);

	std::cout<<signal->SFID[22][0]<<"\n";
	std::cout<<signal->SFID[22][1]<<"\n";  
	//these commands are needed for destructor... TODO: write destructor. 
	//fftw_destroy_plan(p);
    //fftw_free(signal->FID); fftw_free(out);
    //End of T2 Decay FID spectrum stuff

	//start of Bfield 3D fourier transform stuff. 
    std::vector<float> param;
	param.push_back(0.5e-7);
	param.push_back(0.0);
	param.push_back(0.0);
	param.push_back(0.5e-7);
	param.push_back(0.0);
	std::vector<float> pos;
	pos.push_back(-7.6/2.0);
	pos.push_back(7.6/2.0);
	pos.push_back(-10.2/2.0);
	pos.push_back(10.2/2.0);
	pos.push_back(-40/2.0);
	pos.push_back(40/2.0);
	int num=33;
	BField* field= new BField(0.03, param, pos,num);
	int ix,iy,iz,ireal;
	ix=17;
	iy=17;
	iz=17;
	ireal=0; //0 for real part, 1 for imaginary, 
	
	std::cout<<"T2 = "<<field->T2<<" s \n";
	return 0;
}
