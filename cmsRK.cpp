//g++ cmsRK.cpp cmsTransport.cpp cmsBoundary.cpp -I/data1/cmswank/BoostCodeSwank/boost_1_64_0 -o Trajtest -lboost_thread -lboost_system
//use above to compile and/or use to generate a makefile. 
#include "cmsBoundary.h"
#include "cmsTransport.h"
#include "cmsRK.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <boost/thread.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>

cmsRK::cmsRK(int N, double T, double dt, double error, double rkdt,  int startseed, int cores,std::string filename)
{

this->N=N;
this->T=T;
this->dt=dt;
this->error=error;
//this->maxangle=maxangle;
this->startseed=startseed;
this->cores=cores;
this->rkdt=rkdt;
if(this->rkdt>this->dt) this->rkdt=this->dt/5.0;
this->timesteps=(int)(T/dt)+1;
this->spinfile.open (filename.c_str(),std::ios::binary | std::ios::out);
	if(N%cores!=0)
	{
		this->Nactual=N-N%cores;
		std::cout<<"Warning number of spins has been changed to"<<Nactual<<std::endl;
		std::cout<<"If you want to have exactly N make mod(N,cores)=0"<<std::endl;

	}
	else Nactual=N;
//prepare file header. (3 doubles that are really ints but matlab and python can't tell the difference)
double ts=(double)this->timesteps;
double nvars=16.0;
double Nsave=(double)Nactual;
spinfile.write(reinterpret_cast <const char*> (&ts),sizeof(double));
spinfile.write(reinterpret_cast <const char*> (&Nsave),sizeof(double));
spinfile.write(reinterpret_cast <const char*> (&nvars),sizeof(double));

}

void cmsRK::conductorDH(int elem,std::vector<int> belement,\
							double* size,std::vector<double> hole,\
							void (*Temperature)(double*,double*),\
							void (*Convection)(double*,double*),\
							void (*Bfield)(double*,double*),\
							void (*Efield)(double*,double*,double*))
{
	//IO hungry method
	this->Spins=(double*)malloc(cores*timesteps*3*sizeof(double));		
	this->Position=(double*)malloc(cores*timesteps*3*sizeof(double));			
	this->Velocity=(double*)malloc(cores*timesteps*3*sizeof(double));					
	this->BField=(double*)malloc(cores*timesteps*3*sizeof(double));				
	this->EField=(double*)malloc(cores*timesteps*3*sizeof(double));

	for(int ii=0; ii<N-cores; ii+=cores)
	{

		for(int i = 0; i<cores; i++)
		{
			tr.push_back(new cmsTransport(startseed+i+cores*ii,0,belement,size,hole,Temperature,Convection,Bfield,Efield));
			tr[i]->rkdt=this->rkdt; //set rkdt to initial guess. 
			tr[i]->nnn=i; //file io index
			//std::cout<<"seed "<<startseed+i+cores*ii<<"\n";
		}	

		for(int i =0; i<cores; i++)
		{
			bthreads.push_back(new boost::thread(&cmsRK::run,this,tr[i]));
		}

		for(int i =0; i<cores; i++)
		{
			bthreads[i]->join();
			//std::cout<<"joined thread "<<i<<"\n";
		}

		this->saveSpinsDH();
		//std::cout<<"Saving to disk\n";

		tr.erase(tr.begin(), tr.end());
		bthreads.erase(bthreads.begin(), bthreads.end());
		//std::cout<<"running threads\n";		
		if(ii%cores==0) std::cout<<" :) "<<std::flush;
		else if(ii%828==0) std::cout<<" :( "<<std::flush;
		else if(ii%1313==0) std::cout<<" (_8(|) "<<std::flush; 
	
	}
	std::cout<<"\n"<<"Exiting\n";
	spinfile.close();
	return;
}

void cmsRK::conductorMH(bool spinsolve,int elem,std::vector<int> belement,\
							double* size,std::vector<double> hole,\
							void (*Temperature)(double*,double*),\
							void (*Convection)(double*,double*),\
							void (*Bfield)(double*,double*),\
							void (*Efield)(double*,double*,double*))
{
	//Memory Hog yummy pulled pork mmm
	this->Spins=(double*)malloc(Nactual*timesteps*3*sizeof(double));		
	this->Position=(double*)malloc(Nactual*timesteps*3*sizeof(double));			
	this->Velocity=(double*)malloc(Nactual*timesteps*3*sizeof(double));					
	this->BField=(double*)malloc(Nactual*timesteps*3*sizeof(double));				
	this->EField=(double*)malloc(Nactual*timesteps*3*sizeof(double));
	





	for(int i = 0; i<N; i++)
		{
			//initialization hmmm... maybe I did this wrong. but the constructor would be huge otherwise. 
			tr.push_back(new cmsTransport(startseed+i,0,belement,size,hole,Temperature,Convection,Bfield,Efield));
			tr[i]->rkdt=this->rkdt;
			tr[i]->nnn=i;
		}
		
	for(int ii=0; ii<Nactual; ii+=cores)
	{
		if (spinsolve)
		{
			for(int i =ii; i<ii+cores; i++)
			{	
				bthreads.push_back(new boost::thread(&cmsRK::run,this,tr[i]));
			}
		}
	    else
	    {
			for(int i =ii; i<ii+cores; i++)
			{	
				bthreads.push_back(new boost::thread(&cmsRK::runNS,this,tr[i]));
			}
		}
	
		
		for(int i = 0; i<cores; i++)
		{	
			
			bthreads[i]->join();
			//std::cout<<"joined core: "<<i<<"\n";
		}
		//std::cout<<"finished "<<ii<<std::endl;
		bthreads.erase(bthreads.begin(), bthreads.end());
		if(ii%cores==0) std::cout<<" :) "<<std::flush;
		else if(ii%82==0) std::cout<<" :( "<<std::flush;
		else if(ii%131==0) std::cout<<" (_8(|) "<<std::flush; 
	
	}
	std::cout<<"\n"<<"Saving\n";
	this->saveSpins();

	tr.erase(tr.begin(), tr.end()); 
	std::cout<<"exiting\n";
	return;
}


void cmsRK::saveSpins()
{
	double SSX,SSY,SSZ,XX,YY,ZZ,VVX,VVY,VVZ,BBX,BBY,BBZ,EEX,EEY,EEZ,TT;
	int numtemp;
	
	for(int i = 0; i<Nactual; i++)
	{
		for(int ii = 0; ii<timesteps; ii++)
		{
		//flattened array coordinate 
		numtemp=i*timesteps*3+3*ii;

		TT=dt*(double)(ii);
		SSX=this->Spins[numtemp];
		SSY=this->Spins[numtemp+1];
		SSZ=this->Spins[numtemp+2];
		XX=this->Position[numtemp];
		YY=this->Position[numtemp+1];
		ZZ=this->Position[numtemp+2];
		//if (ii == 0) std::cout<<this->Position[numtemp+2]<<"\n";
		VVX=this->Velocity[numtemp];
		VVY=this->Velocity[numtemp+1];
		VVZ=this->Velocity[numtemp+2];	
		BBX=this->BField[numtemp];
		BBY=this->BField[numtemp+1];
		BBZ=this->BField[numtemp+2];
		EEX=this->EField[numtemp];
		EEY=this->EField[numtemp+1];
		EEZ=this->EField[numtemp+2];

		spinfile.write(reinterpret_cast <const char*> (&SSX),sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&SSY),sizeof(double)); 		
		spinfile.write(reinterpret_cast <const char*> (&SSZ),sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&XX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&YY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&ZZ), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVZ), sizeof(double));
		spinfile.write(reinterpret_cast <const char*> (&BBX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&BBY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&BBZ), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&EEX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&EEY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&EEZ), sizeof(double));  
		spinfile.write(reinterpret_cast <const char*> (&TT), sizeof(double));

		}
	}
	spinfile.close();
   	std::cout<<"Saving probably Successful.\n";
	
	return;//format file
 	
}

void cmsRK::saveSpinsDH()
{
	double SSX,SSY,SSZ,XX,YY,ZZ,VVX,VVY,VVZ,BBX,BBY,BBZ,EEX,EEY,EEZ,TT;
	int numtemp;
	
	for(int i = 0; i<cores; i++)
	{
		for(int ii = 0; ii<timesteps; ii++)
		{
		//flattened array coordinate 
		numtemp=i*timesteps*3+3*ii;

		TT=dt*(double)ii;
		SSX=this->Spins[numtemp];
		SSY=this->Spins[numtemp+1];
		SSZ=this->Spins[numtemp+2];
		XX=this->Position[numtemp];
		YY=this->Position[numtemp+1];
		ZZ=this->Position[numtemp+2];
		VVX=this->Velocity[numtemp];
		VVY=this->Velocity[numtemp+1];
		VVZ=this->Velocity[numtemp]+2;	
		BBX=this->BField[numtemp];
		BBY=this->BField[numtemp+1];
		BBZ=this->BField[numtemp+2];
		EEX=this->EField[numtemp];
		EEY=this->EField[numtemp+1];
		EEZ=this->EField[numtemp+2];

		spinfile.write(reinterpret_cast <const char*> (&SSX),sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&SSY),sizeof(double)); 		
		spinfile.write(reinterpret_cast <const char*> (&SSZ),sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&XX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&YY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&ZZ), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVZ), sizeof(double));
		spinfile.write(reinterpret_cast <const char*> (&BBX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&BBY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&BBZ), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&EEX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&EEY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&EEZ), sizeof(double));  
		spinfile.write(reinterpret_cast <const char*> (&TT), sizeof(double));

		}
	}
	
   //	std::cout<<"Saving probably Successful.\n";
	
	return;//format file
 	
}





void cmsRK::run(cmsTransport* trp)
{
	//2D array starting place in flattened 3D array
	int n = trp->nnn*timesteps*3;
	double temp_t;
	int numtemp;


	    this->Spins[n]=trp->spin[0]; 
    	this->Spins[n+1]=trp->spin[1]; 
    	this->Spins[n+2]=trp->spin[2];
    	
    	this->Position[n]=trp->position[0];
    	this->Position[n+1]=trp->position[1];
    	this->Position[n+2]=trp->position[2];
    	//std::cout<<Position[n+2]<<"\n";
    	//std::cout<<"position "<<trp->position[0]<<" "<<trp->position[1]<<" "<<trp->position[2]<<"\n";
    	
    	this->Velocity[n]=trp->velocity[0];
    	this->Velocity[n+1]=trp->velocity[1];
    	this->Velocity[n+2]=trp->velocity[2];

    	this->BField[n]=trp->bfield[0];
    	this->BField[n+1]=trp->bfield[1];
    	this->BField[n+2]=trp->bfield[2];
    	
    	this->EField[n]=trp->efield[0];
    	this->EField[n+1]=trp->efield[1];
    	this->EField[n+2]=trp->efield[2];

	//std::cout<<"I Started "<<trp->nnn<<"\n";
	for(int i = 1 ; i < (int)(T/dt+1); i++)
    {
    	numtemp=n+3*i;
    	temp_t=dt*i;
    	//if(n==0)std::cout<<"time is pre "<<temp_t<<std::endl;
    	//if(n==0)std::cout<<"rkdt "<<trp->rkdt<<"\n";
    	this->integrate(trp,temp_t); //ugg.. this is going to be fun
    	//if(n==0)std::cout<<"time is post "<<temp_t<<std::endl;
		//Save point in time
    	this->Spins[numtemp]=trp->spin[0]; 
    	this->Spins[numtemp+1]=trp->spin[1]; 
    	this->Spins[numtemp+2]=trp->spin[2];
    	
    	this->Position[numtemp]=trp->position[0];
    	this->Position[numtemp+1]=trp->position[1];
    	this->Position[numtemp+2]=trp->position[2];
    	
    	this->Velocity[numtemp]=trp->velocity[0];
    	this->Velocity[numtemp+1]=trp->velocity[1];
    	this->Velocity[numtemp+2]=trp->velocity[2];
    	
    	this->BField[numtemp]=trp->bfield[0];
    	this->BField[numtemp+1]=trp->bfield[1];
    	this->BField[numtemp+2]=trp->bfield[2];
    	
    	this->EField[numtemp]=trp->efield[0];
    	this->EField[numtemp+1]=trp->efield[1];
    	this->EField[numtemp+2]=trp->efield[2];
    }

    //std::cout<<"I made it to the end "<<trp->nnn<<"\n";
    return;
}

void cmsRK::integrate(cmsTransport* tr,double temp_t)
{	//bad integrator testbed, almost as intensive as rk5 for speed tests.
	//will include RK 5th order cash karp in a few days. 
		//int n;

		double temp_t2 = tr->t;
		bool exit=false;
		double temperror;

		///double* tempang;

		
 	while(true) //I dislike while loops, but... it fits!
 	{
 		

 	

		if(temp_t2+tr->rkdt>=temp_t)
		{
			tr->rkdt=temp_t-temp_t2;
			exit=true;

		}
		
		//std::cout<<"time "<<temp_t2<<std::endl;



			 
		
		//for(int i=0; i<100; i++)
		//{	
		/*
		std::cout<<"PRE dt "<<tr->rkdt<<" sx "<<tr->spin[0]<<" sy "<<tr->spin[1]<<" sz "<<tr->spin[2]<<" z pos: "<<tr->position[2]<<"\n";
		std::cout<<"PRE time "<<tr->t<<" twall "<<tr->twall<<" t-twall "<<tr->t-tr->twall<<"\n";
		std::cout<<"PRE position "<<tr->position[0]<<" "<<tr->position[1]<<" "<<tr->position[2]<<"\n";	
		std::cout<<"PRE velocity "<<tr->velocity[0]<<" "<<tr->velocity[1]<<" "<<tr->velocity[2]<<"\n";
		/*	std::cout<<"PRE scatter prob "<<tr->pscatter<<" pcurrent "<<tr->pcurrent<<"\n";
		//std::cin>>n;
		//
 		if(std::isnan(tr->position[2]) || std::isnan(-tr->position[2]))
 		{	
 			std::cout<<"PRE we are outside of the box "<<tr->position[0]<<" "<<tr->position[1]<<" "<<tr->position[2]<<"\n";	
 			std::cout<<"PRE velocity "<<tr->velocity[0]<<" "<<tr->velocity[1]<<" "<<tr->velocity[2]<<"\n";
 			std::cin>>n;
 		}
 		*/
			//solve with the Cash-Karp RK method. 
			tr->solve();
		

		/*
		std::cout<<"POST dt "<<tr->rkdt<<" sx "<<tr->spin[0]<<" sy "<<tr->spin[1]<<" sz "<<tr->spin[2]<<" z pos: "<<tr->position[2]<<"\n";
		std::cout<<"POST time "<<tr->t<<" twall "<<tr->twall<<" t-twall "<<tr->t-tr->twall<<"\n";
		std::cout<<"POST position "<<tr->position[0]<<" "<<tr->position[1]<<" "<<tr->position[2]<<"\n";	
		std::cout<<"POST velocity "<<tr->velocity[0]<<" "<<tr->velocity[1]<<" "<<tr->velocity[2]<<"\n";
		/*std::cout<<"POST scatter prob "<<tr->pscatter<<" pcurrent "<<tr->pcurrent<<"\n";

	if(std::isnan(tr->position[2]) || std::isnan(-tr->position[2]))
 		{	
 			std::cout<<"POST we are outside of the box "<<tr->position[0]<<" "<<tr->position[1]<<" "<<tr->position[2]<<"\n";	
 			std::cout<<"POST velocity "<<tr->velocity[0]<<" "<<tr->velocity[1]<<" "<<tr->velocity[2]<<"\n";
 			std::cin>>n;
 		}

		*/
			//calculate error. 
			DiffMagV(&temperror,tr->cash,tr->karp);			
			//std::cout<<"cash "<<tr->cash[0]<<" "<<tr->cash[1]<<" "<<tr->cash[2]<<"\n";
			//std::cout<<"karp "<<tr->karp[0]<<" "<<tr->karp[1]<<" "<<tr->karp[2]<<"\n";
			//std::cout<<"error "<<temperror<<"\n";

			//
			double tempr= temperror/error; //temp error ratio
			if (tempr<=1.0) 
			{	
				if (tempr<0.82 && !exit) 
				{	
					//		
					temp_t2+=tr->rkdt;
					tr->propagate(tr->rkdt*(tr->a5-tr->a6));
					tr->rkdt*=std::pow(tempr,-.228);;
					Norm(tr->karp);
					tr->spin[0]=tr->karp[0];tr->spin[1]=tr->karp[1];tr->spin[2]=tr->karp[2];
					//break;
				}
				else
				{ 		
						
						temp_t2+=tr->rkdt;
						tr->propagate(tr->rkdt*(tr->a5-tr->a6));
						Norm(tr->karp); 
						tr->spin[0]=tr->karp[0];tr->spin[1]=tr->karp[1];tr->spin[2]=tr->karp[2];
						if(exit){
						
						 return;
						}
						
				}	
			}
			else{
			//error too big. 
			//reset state to beginning and reduce time step. 
			//ughh so many variables define a state? transport is probably too exact. 
			 tr->t=tr->temptime;
			 tr->position[0]=tr->temppos[0];tr->position[1]=tr->temppos[1];tr->position[2]=tr->temppos[2];
			 tr->velocity[0]=tr->tempvel[0];tr->velocity[1]=tr->tempvel[1];tr->velocity[2]=tr->tempvel[2];
			 tr->pscatpast=tr->tempppast;
			 tr->twall=tr->temptwall;
			 tr->pcurrent=tr->temppcurrent;
			 tr->bd->bounce_pos[0]=tr->tempbouncepos[0];tr->bd->bounce_pos[1]=tr->tempbouncepos[1];tr->bd->bounce_pos[2]=tr->tempbouncepos[2];
			 tr->bd->whatwall=tr->tempwhatwall;
			 tr->bd->current_element=tr->tempelement;
			 tr->rkdt*=std::pow(tempr,-0.282);
			 exit=false;
		
			}
		    
		//}
		//std::cout<<"warning: integrator could not converge, retrying\n";
	}
}


void cmsRK::runNS(cmsTransport* trp)
{
	//2D array starting place in flattened 3D array
	int n = trp->nnn*timesteps*3;
	double temp_t;
	int numtemp;


	    this->Spins[n]=0.0;
    	this->Spins[n+1]=0.0;
    	this->Spins[n+2]=0.0;
    	
    	this->Position[n]=trp->position[0];
    	this->Position[n+1]=trp->position[1];
    	this->Position[n+2]=trp->position[2];
    	//std::cout<<Position[n+2]<<"\n";
    	//std::cout<<"position "<<trp->position[0]<<" "<<trp->position[1]<<" "<<trp->position[2]<<"\n";
    	
    	this->Velocity[n]=trp->velocity[0];
    	this->Velocity[n+1]=trp->velocity[1];
    	this->Velocity[n+2]=trp->velocity[2];

    	this->BField[n]=0.0;
    	this->BField[n+1]=0.0;
    	this->BField[n+2]=0.0;
    	
    	this->EField[n]=0.0;
    	this->EField[n+1]=0.0;
    	this->EField[n+2]=0.0;

	//std::cout<<"I Started "<<trp->nnn<<"\n";
	for(int i = 1 ; i < (int)(T/dt+1); i++)
    {
    	numtemp=n+3*i;
    	temp_t=dt*i;
    	
    	//if(n==0)std::cout<<"time is pre "<<temp_t<<std::endl;
    	//if(n==0)std::cout<<"rkdt "<<trp->rkdt<<"\n";
    	
    	//at least ten scatters per mean free path. 
    	double dtp=1.0/trp->lam/10.0;

    	int stepsps=(int)(dt/dtp+1);
    	dtp=dt/((double)stepsps);
    	//std::cout<<"dtp "<<dtp<<"\n";
    	for(int ii=0; ii<stepsps; ii++)
    	{
    	trp->propagate(dtp);
    	}
    	//std::cout<<stepsps*dtp<<"\n";
    	
		//Save point in time
    	this->Spins[numtemp]=0.0; 
    	this->Spins[numtemp+1]=0.0;
    	this->Spins[numtemp+2]=0.0;
    	
    	this->Position[numtemp]=trp->position[0];
    	this->Position[numtemp+1]=trp->position[1];
    	this->Position[numtemp+2]=trp->position[2];
    	
    	this->Velocity[numtemp]=trp->velocity[0];
    	this->Velocity[numtemp+1]=trp->velocity[1];
    	this->Velocity[numtemp+2]=trp->velocity[2];
    	
    	this->BField[numtemp]=0.0;
    	this->BField[numtemp+1]=0.0;
    	this->BField[numtemp+2]=0.0;
    	
    	this->EField[numtemp]=0.0;
    	this->EField[numtemp+1]=0.0;
    	this->EField[numtemp+2]=0.0;
    }

    //std::cout<<"I made it to the end "<<trp->nnn<<"\n";
    return;
}














/* moved to Transport class. 
cmsTransport * cmsRK::getK1(cmsTransport* trp)
{
		trp->bmotional= 1/3E8*Cross(trp->velocity,trp->efield);
		trp->btotal= Add(trp->bfield,trp->bmotional);
		trp->k1 = Cross(trp->spin,trp->btotal)
		trp->k1 = MultiplyS(trp->k1,trp->gamma*trp->rdkt);
		
		return trp; 


}

cmsTransport * cmsRK::getK2(cmsTransport* trp)
{		
		trp->bmotional= 1/3E8*Cross(trp->velocity,trp->efield);
		trp->btotal= Add(trp->bfield,trp->bmotional);
		trp->k2 = Cross(Add(trp->spin,MultiplyS(trp->k1,trp->b21)),trp->btotal );
		trp->k2 = MultiplyS(trp->k2,trp->gamma*trp->rdkt*trp->a2);
		
		return trp; 


}

cmsTransport * cmsRK::getK3(cmsTransport* trp)
{		
		trp->bmotional= 1/3E8*Cross(trp->velocity,trp->efield);
		trp->btotal= Add(trp->bfield,trp->bmotional);
		trp->k3 = Cross(Add(Add(trp->spin,MultiplyS(trp->k1,trp->b31))\
								 ,MultiplyS(trp->k2,trp->b32))\
								 ,trp->btotal );
		
		trp->k3 = MultiplyS(trp->k3,trp->gamma*trp->rdkt*trp->a3);
		
		return trp; 


}
cmsTransport * cmsRK::getK4(cmsTransport* trp)
{		
		trp->bmotional= 1/3E8*Cross(trp->velocity,trp->efield);
		trp->btotal= Add(trp->bfield,trp->bmotional);
		trp->k4 = Cross(Add(Add(Add(trp->spin,MultiplyS(trp->k1,trp->b41))\
								 ,MultiplyS(trp->k2,trp->b42))\
								 ,MultiplyS(trp->k3,trp->b43))\
								 ,trp->btotal );
		
		trp->k4 = MultiplyS(trp->k4,trp->gamma*trp->rdkt*trp->a4);
		
		return trp; 
}
cmsTransport * cmsRK::getK5(cmsTransport* trp)
{		
		trp->bmotional= 1/3E8*Cross(trp->velocity,trp->efield);
		trp->btotal= Add(trp->bfield,trp->bmotional);
		trp->k5 = Cross(Add(Add(Add(Add(trp->spin,MultiplyS(trp->k1,trp->b51))\
								 ,MultiplyS(trp->k2,trp->b52))\
								 ,MultiplyS(trp->k3,trp->b53))\
								 ,MultiplyS(trp->k4,trp->b54))\
						,trp->btotal );
		
		trp->k5 = MultiplyS(trp->k5,trp->gamma*trp->rdkt*trp->a5);
		
		return trp; 
}


cmsTransport * cmsRK::getK6(cmsTransport* trp)
{		
		trp->bmotional= 1/3E8*Cross(trp->velocity,trp->efield);
		trp->btotal= Add(trp->bfield,trp->bmotional);
		trp->k6 = Cross(Add(Add(Add(Add(Add(trp->spin,MultiplyS(trp->k1,trp->b61))\
								 ,MultiplyS(trp->k2,trp->b62))\
								 ,MultiplyS(trp->k3,trp->b63))\
								 ,MultiplyS(trp->k4,trp->b64))\
								 ,MultiplyS(trp->k5,trp->b65))\
								 ,trp->btotal );
		
		trp->k6 = MultiplyS(trp->k6,trp->gamma*trp->rdkt*trp->a6);
		
		return trp; 
}*/