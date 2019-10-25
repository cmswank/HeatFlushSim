#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <boost/random/uniform_real_distribution.hpp>
#include <fftw3.h>

class BField{		
	private:
		std::vector<float> parameters;
		double *Bx,*By,*Bz,*Btemp;
		int N;
		std::vector<float> Bpos(float x, float y, float z);
		void Barray(std::vector<float> pos, int num);
		fftw_plan fieldp;

	public:
		BField(float B0, std::vector<float> params, std::vector<float> pos,int num);
		float T2,df,f;
		fftw_complex *Sx,*Sy,*Sz,*out;
		void bxx();

};