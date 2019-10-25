#include "T2Decay.h"

T2Decay::T2Decay(float T2, float f, float SNR, float tstep,int N):rng(82),dist_noise(0.0f,1.0f),noise(rng,dist_noise){
	this->len=N;
	this->dt=tstep;
	this->FID = (double*) fftw_malloc(sizeof(double)*this->len);
	this->SFID = (fftw_complex*) fftw_malloc(sizeof(double)*this->len);
	this->p = fftw_plan_dft_r2c_1d(this->len, this->FID, this->SFID, FFTW_ESTIMATE);
	if(SNR<=0)
	{	for(int i=0; i<this->len; i++)
		{	
			this->FID[i]=(double)std::exp(-(float)i*this->dt/T2)*std::cos(2*pi*f*this->dt*(float)i);
			//this->FID.push_back(std::exp(-(float)i*dt/T2)*std::cos(2*pi*f*dt*(float)i));					
		}
		this->getSpec();
		return;	
	}	
	else
	{	//why can't I put this in the h file? why does it have to be written twice?
		//boost::mt19937 rng(82); 	
		//boost::normal_distribution<float> dist_noise(0.0f,1.0f);
		//boost::variate_generator<boost::mt19937, boost::normal_distribution<float> > noise(rng, dist_noise);
	
		for(int i = 0; i<this->len; i++)
		{	
			this->FID[i]=(double)std::exp(-(float)i*this->dt/T2)*std::cos(2*pi*f*this->dt*(float)i)+1.0/SNR*this->noise();
			//this->FID.push_back(std::exp(-(float)i*dt/T2)*std::cos(2*pi*f*dt*(float)i)+1/SNR*noise());		
			//this->FID->push_back(std::exp(-(float)i*dt/T2)*std::cos(2*pi*f*dt*(float)i)+1/SNR*noise());		
		}
		this->getSpec();
		return;
	}
}
void T2Decay::getSpec(){
	fftw_execute(this->p);
	return;
}
void T2Decay::update(float T2, float f, float SNR){
	if(SNR<=0)
	{	for(int i=0; i<this->len; i++)
		{	
			this->FID[i]=(double)std::exp(-(float)i*this->dt/T2)*std::cos(2*pi*f*this->dt*(float)i);
			//this->FID.push_back(std::exp(-(float)i*dt/T2)*std::cos(2*pi*f*dt*(float)i));					
		}
		this->getSpec();
		return;	
	}	
	else
	{
		for(int i = 0; i<this->len; i++)
		{	
			this->FID[i]=(double)std::exp(-(float)i*this->dt/T2)*std::cos(2*pi*f*this->dt*(float)i)+1.0/SNR*this->noise();
			//this->FID.push_back(std::exp(-(float)i*dt/T2)*std::cos(2*pi*f*dt*(float)i)+1/SNR*noise());		
			//this->FID->push_back(std::exp(-(float)i*dt/T2)*std::cos(2*pi*f*dt*(float)i)+1/SNR*noise());		
		}
		this->getSpec();
		return;
	}
}

