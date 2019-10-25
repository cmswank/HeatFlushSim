#include "cmsBField.h"

BField::BField(float B0, std::vector<float> params, std::vector<float> pos,int num){
	this->N=num;
	this->f=B0*20378.9;
	assert(params.size()==5);	
	assert(pos.size()==6);
	this->parameters=params;
	this->Bx=(double*) fftw_malloc(sizeof(double)*num*num*num);
	this->By=(double*) fftw_malloc(sizeof(double)*num*num*num);
	this->Bz=(double*) fftw_malloc(sizeof(double)*num*num*num);
	this->Btemp =(double*) fftw_malloc(sizeof(double)*num*num*num);
	this->Sx=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*num*num*num);
	this->Sz=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*num*num*num);
	this->Sy=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*num*num*num);
	this->out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*num*num*num);
	
	this->Barray(pos,this->N);
	
	this->Btemp=this->Bx;
	this->fieldp=fftw_plan_dft_r2c_3d(num, num, num, this->Btemp, out, FFTW_ESTIMATE);
	fftw_execute(fieldp);
	this->Sx=out;
	this->Btemp=this->By;
	fftw_execute(fieldp);
	this->Sy=out;
	this->Btemp=this->Bz;
	fftw_execute(fieldp);
	this->Sz=out;
	this->bxx();

}


void BField::Barray(std::vector<float> pos, int num){
	float xpos=0;
	float ypos=0;
	float zpos=0;
	int it=0;
	std::vector<float> Btemp;
	for(int ix=0; ix<num;ix++){
		xpos=pos[0]+ix*(pos[1]-pos[0])/(N-1);

		for(int iy=0; iy<num;iy++){
			ypos=pos[2]+iy*(pos[3]-pos[2])/(N-1);
			for(int iz=0; iz<num;iz++){
				zpos=pos[4]+iz*(pos[5]-pos[4])/(N-1);
				Btemp=this->Bpos(xpos,ypos,zpos);
				this->Bx[it]=(double)Btemp[0];
				this->By[it]=(double)Btemp[1];
				this->Bz[it]=(double)Btemp[2];
				it++;
			}

		}

	}	
}

std::vector<float> BField::Bpos(float x,float y,float z){
	std::vector<float> Btemp; 
	Btemp.push_back(this->parameters[0]*x+this->parameters[1]*y+this->parameters[2]*z);  //Bx
	Btemp.push_back(this->parameters[3]*y+this->parameters[1]*x+this->parameters[4]*z);  //By
	Btemp.push_back(-(this->parameters[0]+this->parameters[3])*z+this->parameters[2]*x+this->parameters[4]*y);  //Bz
	return Btemp;
}

void BField::bxx(){
	double tempT2=0.0;
	double tempdf=0.0;
	double numtemp=(double)this->N*this->N*this->N;
	//3D normalization
	double normtemp=8.0*8.0/numtemp/numtemp;
	for(int i=0; i<this->N*this->N*this->N;i++)
	{
		tempT2+=normtemp*this->Sx[i][0]*this->Sx[i][0];
		tempdf+=normtemp*(this->Sy[i][1]*this->Sy[i][1]+this->Sz[i][1]*this->Sz[i][1]);
	}
	
									
	this->T2=2.0/20378.9/20378.9/tempT2;
	this->df=20378.9*20378.9/2.0*5.3e2*tempdf;
	return;

}
