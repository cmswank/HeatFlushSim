#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <boost/random/uniform_real_distribution.hpp>
#include <fftw3.h>

class T2Decay{
	private:
		typedef boost::mt19937 RNG;    // Mersenne Twister
		typedef boost::normal_distribution<float> DIST_norm;   // Normal Distribution
		typedef boost::variate_generator<RNG,DIST_norm> NOISE;    // Noise generator. 

     RNG rng;
     DIST_norm dist_noise;
     NOISE noise;

		static const double pi = 3.141592653589793238462643383;
		void getSpec();
		int len;
		float dt;
		fftw_plan p;
	public:		
		double* FID;
		fftw_complex* SFID;
		T2Decay(float T2, float f, float SNR, float tstep, int N);
		void update(float T2, float f, float SNR);
			
};
