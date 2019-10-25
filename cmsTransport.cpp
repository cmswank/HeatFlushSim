#include "cmsBoundary.h"
#include "cmsTransport.h"
#include <assert.h>
#include <iostream>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>

cmsTransport::cmsTransport(int flat_seed,int elem,std::vector<int> belement,\
							double* size,std::vector<double> hole,\
							void (*Temperature)(double*,double*),\
							void (*Convection)(double*,double*),\
							void (*Bfield)(double*,double*),\
							void (*Efield)(double*,double*,double*))\
							:rng(flat_seed),dist_flat(0.0f,1.0f),rand_flat(rng,dist_flat),\
							dist_vel(0.0f,6.0f),rand_vel(rng,dist_vel)
{
	
	bd=new cmsBoundary(belement,size,hole);
	this->Temperature=Temperature;
	this->Convection=Convection;
	this->Bfield=Bfield;
	this->Efield=Efield;
	

	this->tempSpin= (double*)new double[3];
	this->velocity=  (double*)new double[3];
	this->spin=      (double*)new double[3];
	this->bfield=    (double*)new double[3];
	this->efield=    (double*)new double[3];
	this->bmotional= (double*)new double[3];
	this->btotal=    (double*)new double[3];
	this->cash=      (double*)new double[3];
	this->karp=      (double*)new double[3];
	
	this->temppos=  (double*) new double[3];
	this->tempvel=   (double*) new double[3];
	this->tempbouncepos=   (double*) new double[3];

	this->tempU=(double*)new double[3];
	this->tempT=(double*)new double[1];


	this->getInitial(elem);
	this->tempSpin=new double(*this->spin);
	this->tempPos=new double(*this->position);
	this->tempVel=new double(*this->velocity);
	this->tempWallt=this->twall;
	this->tempBouncePos=   (double*) new double[3];

	//RK constants. 
	this->k1 = (double*)new double[3];
	this->k2 = (double*)new double[3];
	this->k3 = (double*)new double[3];
	this->k4 = (double*)new double[3];
	this->k5 = (double*)new double[3];
	this->k6 = (double*)new double[3];


	this->pscatter=rand_flat();

}


void cmsTransport::getInitial(int elem){

	
	do
	{

	double randx=rand_flat();
	double randy=rand_flat();
	double randz=rand_flat();

	this->position=bd->getInitPos(elem,randx,randy,randz);
	
	
	}while(bd->element[elem]==0 && std::sqrt(position[0]*position[0]+position[1]*position[1])>bd->sizes[elem*3]);

	this->Temperature(tempT,position);
	lam=k/0.00016/m*std::pow(tempT[0],8);
	pscatnow=0.0;
	//std::cout<<"tau_c = "<<1.0/lam<<"\n";
	pscatpast=lam;
	this->t=0.0;

	double sigtemp=std::sqrt(k*tempT[0]/m);
	DIST_velocity dist_vel(0.0f,sigtemp);
	rand_vel.distribution()=dist_vel;
	
	this->Convection(tempU,position);
	this->velocity[0]=rand_vel()+tempU[0];
	this->velocity[1]=rand_vel()+tempU[1];
	this->velocity[2]=rand_vel()+tempU[2];
	this->twall=bd->nextBounce(this->position,this->velocity);
	this->Bfield(this->bfield,this->position);
	this->Efield(this->efield,this->position,this->velocity);
	this->spin[0]=0.0;
	this->spin[1]=0.0;
	this->spin[2]=1.0;
	//std::cout<<"Twall is = "<<twall<<"\n";
}

void cmsTransport::propagate(double dt){
	bool nohole=true;
	t+=dt;
	//std::cout<<"start time is "<<t<<" t wall is "<<twall<<"\n";
	if(t>=twall){

		if (bd->whatwall == 2 && std::sqrt(std::pow(bd->bounce_pos[0],2)+std::pow(bd->bounce_pos[1],2))<bd->holes[bd->current_element] && sgn(velocity[2])>0)
		{	//int ntest;
			bd->current_element++;
			//std::cout<<"we made it to the hole!\n";
			dt=t-twall;
			bd->benter=true;
			this->position[0]=bd->bounce_pos[0];this->position[1]=bd->bounce_pos[1];this->position[2]=bd->bounce_pos[2];
			twall+=bd->nextBounce(bd->bounce_pos,velocity);
			
			//std::cout<<"hole loop time is "<<t<<" dt is "<<dt<<" t wall is "<<twall<<"\n";
			nohole=false;
		}
		else if (bd->current_element>0)
		{
			if(bd->whatwall == 2 && std::sqrt(std::pow(bd->bounce_pos[0],2)+std::pow(bd->bounce_pos[1],2))<bd->holes[bd->current_element-1] && sgn(velocity[2])<0)
			{
				bd->current_element--;
				dt=t-twall;
				bd->bleave=true; //tell the boundary we are in a new boundary. 
				this->position[0]=bd->bounce_pos[0];this->position[1]=bd->bounce_pos[1];this->position[2]=bd->bounce_pos[2];
				twall+=bd->nextBounce(bd->bounce_pos,velocity);
				
				nohole=false;

			}
		}

		while (t>=twall){

	//if(this->nnn==0)std::cout<<"We made it to the wall! time= "<<t<<", vx = "<<velocity[0]<<"\n";
	
	this->Temperature(tempT,bd->bounce_pos);
	double sigtemp=std::sqrt(k*tempT[0]/m);
	DIST_velocity dist_vel(0.0,sigtemp);
	rand_vel.distribution()=dist_vel;
	double r1=rand_vel();
	double r2=rand_vel();
	double r3=rand_vel();
	double randv=std::sqrt(r1*r1+r2*r2+r3*r3);
	double rand1=rand_flat();
	double rand2=rand_flat();
	//std::cout<<"twall loop Before a bounce, t "<<t<<" twall "<<twall<<" dt "<< dt <<" z "<<"\n";
	dt=t-twall;
	
	twall+=bd->getBounce(rand1,rand2,randv,position,velocity);	
     //std::cout<<"twall loop, After a Bounce t "<<t<<" twall "<<twall<<" dt "<< dt <<" z "<<"\n";
		}
	//	std::cout<<"We made it PAST THE WALL!\n"; // posx = "<<position[0]<<", vx = "<<velocity[0]<<"\n";
	if (nohole){
		position[0]=position[0]+velocity[0]*dt;
		position[1]=position[1]+velocity[1]*dt;
		position[2]=position[2]+velocity[2]*dt;
		this->Temperature(tempT,position);
		lam=k/0.00016/m*std::pow(tempT[0],8.0);
		pscatnow=lam*std::exp(-dt*lam);
		pcurrent+=dt*(pscatnow);
		pscatpast=pscatnow;
	//twalltest=twall+bd->nextBounce(position,velocity);
		//double rmag=std::sqrt(std::pow(position[0],2.0)+std::pow(position[1],2.0));
//	std::cout<<"twall loop, r mag= "<<rmag<<" x "<<position[0]<<" y "<<position[1]<<" z "<<position[2]<<"\n";
	//std::cout<<"twall loop, t "<<t<<" twall "<<twall<<" dt "<< dt <<" z "<<"\n";

		return;
		}
	}
	

//std::cout<<" pre propagate loop time is "<<t<<" t wall is "<<twall<<"\n";
//std::cout<<"position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
//std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
//std::cout<<"element "<<bd->current_element<<"\n";
position[0]=position[0]+velocity[0]*dt;
position[1]=position[1]+velocity[1]*dt;
position[2]=position[2]+velocity[2]*dt;

this->Temperature(tempT,position);
lam=k/0.00016/m*std::pow(tempT[0],8.0);
pscatnow=lam*std::exp(-dt*lam);
pcurrent+=dt*(pscatpast+pscatnow);
pscatpast=pscatnow;
//std::cout<<pcurrent<<"\n";
//std::cout<<" post propagate loop time is "<<t<<" t wall is "<<twall<<"\n";
//std::cout<<"position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
//std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
if(pcurrent>pscatter  && nohole) {
//	std::cout<<"We scattered!\n";
	//double rmag=std::sqrt(std::pow(position[0],2.0)+std::pow(position[1],2.0));
	//std::cout<<"r mag= "<<rmag<<" x "<<position[0]<<" y "<<position[1]<<" z "<<position[2]<<"\n";
	pscatter=rand_flat();
	pcurrent=0.0;
	pscatpast=k/0.00016/m*std::pow(tempT[0],8.0);
	lam=k/0.00016/m*std::pow(tempT[0],8.0);
	//std::cout<<"lam = "<<lam<<"\n";
	//m/2/k/tempT
	double sigtemp=std::sqrt(k*tempT[0]/m);
	DIST_velocity dist_vel(0.0f,sigtemp);
	rand_vel.distribution()=dist_vel;
	this->Convection(tempU,position);
	velocity[0]=rand_vel()+tempU[0];
	velocity[1]=rand_vel()+tempU[1];
	velocity[2]=rand_vel()+tempU[2];
	//std::cout<<"scatter loop, before a Scatter, t "<<t<<" twall "<<twall<<" dt "<< dt <<" z "<<"\n";
	
	twall=t+bd->nextBounce(position,velocity);
	
//	std::cout<<"scatter loop, After a scatter, t "<<t<<" twall "<<twall<<" dt "<< dt <<" z "<<"\n";
	//std::cout<<"v_z = "<<velocity[2]<<"\n";
	//change velocity with convection current.
	}
	//std::cout<<" post scatter loop time is "<<t<<" t wall is "<<twall<<"\n";
  //	std::cout<<"position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
	//std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";

   return;
}

//I was planning on putting this in RK class, but it fits more naturally here with multi-threading in the end.
void cmsTransport::solve()
{	//int n;
		/*
		if(this->nnn==0 )
		{
			std::cout<<"Spin "<<this->spin[0]<<" "<<this->spin[1]<<" "<<this->spin[2]<<"\n";
			std::cout<<"RK pre position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
			std::cout<<"RK pre velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
		
		std::cout<<"RK pre dt "<<rkdt<<" sx "<<spin[0]<<" sy "<<spin[1]<<" sz "<<spin[2]<<" z pos: "<<position[2]<<"\n";
		std::cout<<"RK pre time "<<t<<" twall "<<twall<<"\n";
		if(std::isnan(rkdt))
			{
			
			std::cin>>n;
			}
		
		}*/
		
		/*
		std::cout<<"RK pre position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
		std::cout<<"RK pre scatter prob "<<pscatter<<" pcurrent "<<pcurrent<<"\n";

	if(std::isnan(position[2]) || std::isnan(-position[2]))
 		{	
 			std::cout<<"RK pre we are outside of the box "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
 			std::cout<<"RK pre velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
 			std::cin>>n;
 		}
		*/


		///this is the complete 'state' of the trajectory. 
		this->tempPos[0]=this->position[0];this->tempPos[1]=this->position[1];this->tempPos[2]=this->position[2];
		this->tempVel[0]= this->velocity[0];this->tempVel[1]= this->velocity[1];this->tempVel[2]= this->velocity[2];
		this->tempSpin[0]=this->spin[0];this->tempSpin[1]=this->spin[1];this->tempSpin[2]=this->spin[2];
		this->tempBouncePos[0]=this->bd->bounce_pos[0];this->tempBouncePos[1]=this->bd->bounce_pos[1];this->tempBouncePos[2]=this->bd->bounce_pos[2];
		this->tempWhatWall=this->bd->whatwall;
		tempWallt=this->twall;
		tempScatp=this->pcurrent;
		tempTime=this->t;
		this->tempppast=this->pscatpast;
		this->tempScatold=this->pscatpast;
		this->temptime= this->t;
		this->temppos[0]=this->position[0];this->temppos[1]=this->position[1];this->temppos[2]=this->position[2];
		this->tempvel[0]= this->velocity[0];this->tempvel[1]= this->velocity[1];this->tempvel[2]= this->velocity[2];
		this->tempbouncepos[0]=this->bd->bounce_pos[0];this->tempbouncepos[1]=this->bd->bounce_pos[1];this->tempbouncepos[2]=this->bd->bounce_pos[2];
		this->temptwall= this->twall;
		this->temppcurrent= this->pcurrent;
		this->tempwhatwall=this->bd->whatwall;
		this->tempelement=this->bd->current_element;
		this->tempElement=this->bd->current_element;
		//this gets the complete new 'state' of the trajectory. 
		
		//std::cout<<"pre k1 position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	

		this->propagate(this->rkdt);
		Bfield(this->bfield,this->position);
		Efield(this->efield,this->position,this->velocity);
		Cross(this->bmotional,this->velocity,this->efield);
		this->bmotional[0]*=1.0/9.0e16;
		this->bmotional[1]*=1.0/9.0e16;
		this->bmotional[2]*=1.0/9.0e16;
		//get k1

		this->btotal[0]=this->bmotional[0]+this->bfield[0];
		this->btotal[1]=this->bmotional[1]+this->bfield[1];
		this->btotal[2]=this->bmotional[2]+this->bfield[2];
		///write it all out, no functions except for cross please...
		Cross(this->k1,this->spin,this->btotal);
		this->k1[0]*=this->gamma*this->rkdt;
		this->k1[1]*=this->gamma*this->rkdt;
		this->k1[2]*=this->gamma*this->rkdt;
		//std::cout<<"Btotal "<<this->btotal[0]<<" "<<this->btotal[1]<<" "<<this->btotal[2]<<"\n";
		//std::cout<<"gamma "<<this->gamma<<std::endl;
		//std::cout<<"rkdt "<<this->rkdt<<std::endl;
		//std::cout<<"dt "<<this->rkdt<<"\n";
		//std::cout<<"k1: "<<" "<<k1[0]<<" "<<k1[1]<<" "<<k1[2]<<"\n";
		//std::cin>>n;
		//return to starting state to get k2 through k5. grab the values, not the pointer. 
		//this->spin[0]=tempSpin[0];this->spin[1]=tempSpin[1];this->spin[2]=tempSpin[2];
		this->position[0]=tempPos[0];this->position[1]=tempPos[1];this->position[2]=tempPos[2];
		this->velocity[0]=tempVel[0];this->velocity[1]=tempVel[1];this->velocity[2]=tempVel[2];
		this->twall=tempWallt;
		this->bd->bounce_pos[0]=this->tempBouncePos[0];this->bd->bounce_pos[1]=this->tempBouncePos[1];this->bd->bounce_pos[2]=this->tempBouncePos[2];
		this->bd->whatwall=this->tempWhatWall;
		this->bd->current_element=this->tempElement;
		//std::cout<<"twall "<<twall<<"\n";
		this->pcurrent=tempScatp;
		this->t=tempTime;
		this->pscatpast=this->tempScatold;
		//std::cout<<"pre k2 position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
		//propagate state to a2 and get fields. 
		this->propagate(this->rkdt*this->a2);
		Bfield(this->bfield,this->position);
		Efield(this->efield,this->position,this->velocity);
		Cross(this->bmotional,this->velocity,this->efield);
		this->bmotional[0]*=1.0/9.0e16;
		this->bmotional[1]*=1.0/9.0e16;
		this->bmotional[2]*=1.0/9.0e16;
		this->btotal[0]=this->bmotional[0]+this->bfield[0];
		this->btotal[1]=this->bmotional[1]+this->bfield[1];
		this->btotal[2]=this->bmotional[2]+this->bfield[2];

		//get k2
		
		tempSpin[0]=this->spin[0]+this->k1[0]*this->b21;
		tempSpin[1]=this->spin[1]+this->k1[1]*this->b21;
		tempSpin[2]=this->spin[2]+this->k1[2]*this->b21;

		//std::cout<<"Spin "<<this->spin[0]<<" "<<this->spin[1]<<" "<<this->spin[2]<<"\n";
		//std::cout<<"TempSpin "<<this->tempSpin[0]<<" "<<this->tempSpin[1]<<" "<<this->tempSpin[2]<<"\n";


		Cross(this->k2,this->tempSpin,this->btotal);
		
		this->k2[0]*=this->gamma*this->rkdt;
		this->k2[1]*=this->gamma*this->rkdt;
		this->k2[2]*=this->gamma*this->rkdt;
		/*
			std::cout<<"RK mid dt "<<rkdt<<" sx "<<spin[0]<<" sy "<<spin[1]<<" sz "<<spin[2]<<" z pos: "<<position[2]<<"\n";
		std::cout<<"RK mid time "<<t<<" twall "<<twall<<"\n";
		std::cout<<"RK mid position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
		std::cout<<"RK mid scatter prob "<<pscatter<<" pcurrent "<<pcurrent<<"\n";

	if(std::isnan(position[2]) || std::isnan(-position[2]))
 		{	
 			std::cout<<"RK mid we are outside of the box "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
 			std::cout<<"RK mid velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
 			std::cin>>n;
 		}
	
	*/
		//std::cout<<"k2: "<<" "<<k2[0]<<" "<<k2[1]<<" "<<k2[2]<<"\n";
		//std::cin>>n;
		//propagate state to a3 from a2;
		//propagate state to a3 and get fields. 
		//std::cout<<"pre k3 position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
		this->propagate(this->rkdt*(this->a3-this->a2));
		Bfield(this->bfield,this->position);
		Efield(this->efield,this->position,this->velocity);
		Cross(this->bmotional,this->velocity,this->efield);
		this->bmotional[0]*=1.0/9.0e16;
		this->bmotional[1]*=1.0/9.0e16;
		this->bmotional[2]*=1.0/9.0e16;
		this->btotal[0]=this->bmotional[0]+this->bfield[0];
		this->btotal[1]=this->bmotional[1]+this->bfield[1];
		this->btotal[2]=this->bmotional[2]+this->bfield[2];

		//get k3
		
		tempSpin[0]=this->spin[0]+this->k2[0]*this->b32+this->k1[0]*this->b31;
		tempSpin[1]=this->spin[1]+this->k2[1]*this->b32+this->k1[1]*this->b31;
		tempSpin[2]=this->spin[2]+this->k2[2]*this->b32+this->k1[2]*this->b31;

		Cross(this->k3,this->tempSpin,this->btotal);
		
		this->k3[0]*=this->gamma*this->rkdt;
		this->k3[1]*=this->gamma*this->rkdt;
		this->k3[2]*=this->gamma*this->rkdt;

		//std::cout<<"k3: "<<" "<<k3[0]<<" "<<k3[1]<<" "<<k3[2]<<"\n";
		//std::cin>>n;




		//propagate state to a4 from a3;
		//propagate state to a4 and get fields. 
		//std::cout<<"pre k4 position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
		this->propagate(this->rkdt*(this->a4-this->a3));
		Bfield(this->bfield,this->position);
		Efield(this->efield,this->position,this->velocity);
		Cross(this->bmotional,this->velocity,this->efield);
		this->bmotional[0]*=1.0/9.0e16;
		this->bmotional[1]*=1.0/9.0e16;
		this->bmotional[2]*=1.0/9.0e16;
		this->btotal[0]=this->bmotional[0]+this->bfield[0];
		this->btotal[1]=this->bmotional[1]+this->bfield[1];
		this->btotal[2]=this->bmotional[2]+this->bfield[2];

		//get k4
		
		tempSpin[0]=this->spin[0]+this->k3[0]*this->b43+this->k2[0]*this->b42+this->k1[0]*this->b41;
		tempSpin[1]=this->spin[1]+this->k3[1]*this->b43+this->k2[1]*this->b42+this->k1[1]*this->b41;
		tempSpin[2]=this->spin[2]+this->k3[2]*this->b43+this->k2[2]*this->b42+this->k1[2]*this->b41;

		Cross(this->k4,this->tempSpin,this->btotal);
		
		this->k4[0]*=this->gamma*this->rkdt;
		this->k4[1]*=this->gamma*this->rkdt;
		this->k4[2]*=this->gamma*this->rkdt;

		//std::cout<<"k4: "<<" "<<k4[0]<<" "<<k4[1]<<" "<<k4[2]<<"\n";
		//std::cin>>n;

		///save state for k6, (k6 is less than k5 for some crappy reason.)
		this->tempPos[0]=this->position[0];tempPos[1]=this->position[1];tempPos[2]=this->position[2];
		this->tempVel[0]=this->velocity[0];tempVel[1]=this->velocity[1];tempVel[2]=this->velocity[2];
		this->tempWhatWall=this->bd->whatwall;
		this->tempWallt=this->twall;
		this->tempScatp=this->pcurrent;
		this->tempScatold=this->pscatpast;
		this->tempTime=this->t;
		this->tempBouncePos[0]=this->bd->bounce_pos[0];this->tempBouncePos[1]=this->bd->bounce_pos[1];this->tempBouncePos[2]=this->bd->bounce_pos[2];
		this->tempElement=this->bd->current_element;



		//propagate state to a5 from a4;
		//propagate state to a5 and get fields. 
		this->propagate(this->rkdt*(this->a5-this->a4));
		Bfield(this->bfield,this->position);
		Efield(this->efield,this->position,this->velocity);
		Cross(this->bmotional,this->velocity,this->efield);
		this->bmotional[0]*=1.0/9.0e16;
		this->bmotional[1]*=1.0/9.0e16;
		this->bmotional[2]*=1.0/9.0e16;
		this->btotal[0]=this->bmotional[0]+this->bfield[0];
		this->btotal[1]=this->bmotional[1]+this->bfield[1];
		this->btotal[2]=this->bmotional[2]+this->bfield[2];

		/*
		if(this->nnn==0){
		std::cout<<"RK post dt "<<rkdt<<" sx "<<spin[0]<<" sy "<<spin[1]<<" sz "<<spin[2]<<" z pos: "<<position[2]<<"\n";
		std::cout<<"RK post time "<<t<<" twall "<<twall<<"\n";
		std::cout<<"RK post position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
		std::cout<<"RK post scatter prob "<<pscatter<<" pcurrent "<<pcurrent<<"\n";

	if(std::isnan(position[2]) || std::isnan(-position[2]))
 		{	
 			std::cout<<"RK post we are outside of the box "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";	
 			std::cout<<"RK post velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
 			//	std::cin>>n;
 		}
		}*/



		//get k5
		
		tempSpin[0]=this->spin[0]+this->k4[0]*this->b54+this->k3[0]*this->b53+this->k2[0]*this->b52+this->k1[0]*this->b51;
		tempSpin[1]=this->spin[1]+this->k4[1]*this->b54+this->k3[1]*this->b53+this->k2[1]*this->b52+this->k1[1]*this->b51;
		tempSpin[2]=this->spin[2]+this->k4[2]*this->b54+this->k3[2]*this->b53+this->k2[2]*this->b52+this->k1[2]*this->b51;

		Cross(this->k5,this->tempSpin,this->btotal);
		
		this->k5[0]*=this->gamma*this->rkdt;
		this->k5[1]*=this->gamma*this->rkdt;
		this->k5[2]*=this->gamma*this->rkdt;

		//std::cout<<"k5: "<<" "<<k5[0]<<" "<<k5[1]<<" "<<k5[2]<<"\n";
		//std::cin>>n;



		//return state to a4
		this->position[0]=tempPos[0];this->position[1]=tempPos[1];this->position[2]=tempPos[2];
		this->velocity[0]=tempVel[0];this->velocity[1]=tempVel[1];this->velocity[2]=tempVel[2];
		this->bd->bounce_pos[0]=this->tempBouncePos[0];this->bd->bounce_pos[1]=this->tempBouncePos[1];this->bd->bounce_pos[2]=this->tempBouncePos[2];
		this->bd->whatwall=this->tempWhatWall;
		this->twall=this->tempWallt;
		this->pcurrent=this->tempScatp;
		this->pscatpast=this->tempScatold;
		this->t=tempTime;
		this->bd->current_element=this->tempElement;
		//propagate to a6 from a4 and get fields
		this->propagate(this->rkdt*(this->a6-this->a4));
		Bfield(this->bfield,this->position);
		Efield(this->efield,this->position,this->velocity);
		Cross(this->bmotional,this->velocity,this->efield);
		this->bmotional[0]*=1.0/9.0e16;
		this->bmotional[1]*=1.0/9.0e16;
		this->bmotional[2]*=1.0/9.0e16;
		this->btotal[0]=this->bmotional[0]+this->bfield[0];
		this->btotal[1]=this->bmotional[1]+this->bfield[1];
		this->btotal[2]=this->bmotional[2]+this->bfield[2];

		//get k6
		
		tempSpin[0]=this->spin[0]+this->k5[0]*this->b65+this->k4[0]*this->b64+this->k3[0]*this->b63+this->k2[0]*this->b62+this->k1[0]*this->b61;
		tempSpin[1]=this->spin[1]+this->k5[1]*this->b65+this->k4[1]*this->b64+this->k3[1]*this->b63+this->k2[1]*this->b62+this->k1[1]*this->b61;
		tempSpin[2]=this->spin[2]+this->k5[2]*this->b65+this->k4[2]*this->b64+this->k3[2]*this->b63+this->k2[2]*this->b62+this->k1[2]*this->b61;

		Cross(this->k6,this->tempSpin,this->btotal);
		
		this->k6[0]*=this->gamma*this->rkdt;
		this->k6[1]*=this->gamma*this->rkdt;
		this->k6[2]*=this->gamma*this->rkdt;

		//std::cout<<"k6: "<<" "<<k6[0]<<" "<<k6[1]<<" "<<k6[2]<<"\n";
		//std::cin>>n;

		this->cash[0]=this->spin[0]+this->k1[0]*this->c1+this->k3[0]*this->c3+this->k4[0]*this->c4+this->k6[0]*this->c6;
		this->cash[1]=this->spin[1]+this->k1[1]*this->c1+this->k3[1]*this->c3+this->k4[1]*this->c4+this->k6[1]*this->c6;
		this->cash[2]=this->spin[2]+this->k1[2]*this->c1+this->k3[2]*this->c3+this->k4[2]*this->c4+this->k6[2]*this->c6;

		this->karp[0]=this->spin[0]+this->k1[0]*this->d1+this->k3[0]*this->d3+this->k4[0]*this->d4+this->k5[0]*this->d5+this->k6[0]*this->d6;
		this->karp[1]=this->spin[1]+this->k1[1]*this->d1+this->k3[1]*this->d3+this->k4[1]*this->d4+this->k5[1]*this->d5+this->k6[1]*this->d6;
		this->karp[2]=this->spin[2]+this->k1[2]*this->d1+this->k3[2]*this->d3+this->k4[2]*this->d4+this->k5[2]*this->d5+this->k6[2]*this->d6;
		
		return; 


}









/*

//temperature function!
double Temperature(std::vector<double> pos){

		return 0.45;
}

std::vector<double> Bfield(std::vector<double> pos){
	std::vector<double> temp;
	temp.push_back(0);
	temp.push_back(0);
	temp.push_back(3e-6);

	return temp;  //30mG in Tesla
}

std::vector<double> Efield(std::vector<double> pos, std::vector<double> vel){
	std::vector<double> temp;
	temp.push_back(0);
	temp.push_back(0);
	temp.push_back(75E7);

	return temp; //75 kV/cm in kV/m
}

std::vector<double> Convection(std::vector<double> pos){
	std::vector<double> temp;
	temp.push_back(0.0f);
	temp.push_back(0.0f);
	temp.push_back(-0.5f);

	return temp; //75 kV/cm 
}

//double sgn function. 
double sgn(double x) {
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}


/*
int main(){

	std::vector<int> belem;
	std::vector<double*> bsize;
	std::vector<double> hole;
	belem.push_back(0);
	double tempsize[1][3]= {1.0,1.0,1.0};
	//bsize.push_back((double*)tempsize);
	hole.push_back(-1.0);
	//double temp23=(double)bsize[0];	
	double* temp23=(double*)tempsize;
	double temp24=temp23[1];

	cmsTransport * tr= new cmsTransport(13,0,belem,(double*)tempsize,hole,Temperature,Convection,Bfield,Efield);
	
	
	
	for(int i = 0 ; i < 1E7; ++i)
    {	

    	tr->propagate(2e-6);
 		if (i%(int)1E6 ==0) std::cout<<"hello we made it to "<<i/((int)1E6)<<" million time steps\n";

    }
    return 0;

}
*/