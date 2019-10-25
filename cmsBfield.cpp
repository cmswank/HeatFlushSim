#Include "cmsBField.h"

BField::BField(float B0, std::vector<float> params, std::vector<float> pos,int num){
	this->N=num;
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
}


void BField::Barray(std::vector<float> pos, int N){
	float xpos=0;
	float ypos=0;
	float zpos=0;
	int it=0;
	std::vector<float> Btemp;
	for(int ix=0; ix<N;ix++){
		xpos=pos[0]+ix*(pos[1]-pos[0])/(N-1);
		for(int iy=0; iy<N;iy++){
			ypos=pos[2]+iy*(pos[3]-pos[2])/(N-1);
			for(int iz=0; iz<N;iz++){
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
	Btemp.push_back(this->parameters[0]*x+this->parameters[1]*y+this->parameters[3]*z);  //Bx
	Btemp.push_back(this->parameters[4]*y+this->parameters[1]*x+this->parameters[5]*z);  //By
	Btemp.push_back(-(this->parameters[0]+this->parameters[4])*z+this->parameters[3]*x+this->parameters[5]*y);  //Bz
	return Btemp;
}
