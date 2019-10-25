#include "cmsBoundary.h"
#include <assert.h>
#include <iostream>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>





cmsBoundary::cmsBoundary(std::vector<int> belement,double* bsize,std::vector<double> hole){
	this->element=belement;
	this->holes=hole;
	this->sizes=bsize;
	this->tempv=(double*)new double[3];
	this->bounce_pos=(double*)new double[3];
	this->benter=false;
	this->bleave=false;
	//TODO: add bou
	//c++ is essentially broken...  :.(    its pretty fast though.
	int i=0;
	double tempsize=0.0;
	for(std::vector<int>::iterator it = belement.begin(); it != belement.end(); ++it) {
    
    bstart.push_back(tempsize);

    tempsize+=bsize[3*i+2];
    i++;
	
	}

	bstart.push_back(0.0);
	//bounce_pos[0]push_back(0.0);
	//bounce_pos.push_back(0.0);
	//bounce_pos.push_back(0.0);

}

double cmsBoundary::getBounce(double rand0, double rand1, double speed,double *position, double *velocity){
	//there is no slip if it actually hits the wall. doesn't make sense. 
		//int n; //debugging stuff...
	if(whatwall==0){
		//x wall. 
		double thetap=std::asin(std::sqrt(rand0));
		velocity[0]=-sgn(velocity[0])*speed*std::cos(thetap); //make negative from reflection.
		velocity[1]=speed*std::sin(thetap)*std::sin(2*pi*rand1);
		velocity[2]=speed*std::sin(thetap)*std::cos(2*pi*rand1);

	}
	else if(whatwall==1){
		//ywall
		
		double thetap=std::asin(std::sqrt(rand0));
		velocity[1]=-sgn(velocity[1])*speed*std::cos(thetap); //make negative from reflection.
		velocity[2]=speed*std::sin(thetap)*std::sin(2*pi*rand1);
		velocity[0]=speed*std::sin(thetap)*std::cos(2*pi*rand1);
	}
	else if(whatwall==2){
		//double theta=std::asin(std::sqrt(rand0));
		//double phi=2*pi*rand1;
		//cos(theta) dist, (with sin(theta) area element required in 3D). 
		double thetap=std::asin(std::sqrt(rand0));
		velocity[2]=-sgn(velocity[2])*speed*std::cos(thetap); //make negative from reflection.
		velocity[0]=speed*std::sin(thetap)*std::sin(2*pi*rand1);
		velocity[1]=speed*std::sin(thetap)*std::cos(2*pi*rand1);


		//debugging....
		//std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
		//std::cout<<"Theta = "<<theta<<"\n";
		//std::cout<<"position x "<<bounce_pos[0]<<" y "<<bounce_pos[1]<<" z "<<bounce_pos[2]<<"\n";
		
		//velocity[2]=-sgn(velocity[2])*speed*std::cos(theta);
		//velocity[0]=speed*std::sin(theta)*std::cos(phi);
		//velocity[1]=speed*std::sin(theta)*std::sin(phi);
		
		//zwall
		//std::cout<<"we hit the z wall !!!!\n";
		//std::cout<<"velocity vx "<<velocity[0]<<" vy "<<velocity[1]<<" vz "<<velocity[2]<<"\n";
		//std::cin>> n;
	}
	else{
		//whatwall ==4 cylinder radius bounce. 

		
		
		//	std::cout<<"we hit the radius with velocity!!!!\n";
		//std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
		double radius = sizes[3*current_element];
		double xbounce = bounce_pos[0];
		double ybounce = bounce_pos[1];
		
		//cos(theta)^2*sin(theta) distribution (cos^2)
		//double thetap=std::acos(std::pow(rand0,0.333333333333333));

		//cos(theta) dist
		double thetap=std::asin(std::sqrt(rand0));
		
		//sqrt(cos(theta)) dist
		//double thetap=std::acos(std::pow(rand0,0.66666666666666667));
		double vr=-speed*std::cos(thetap); //make negative from reflection.
		double vp=speed*std::sin(thetap)*std::sin(2*pi*rand1);
		velocity[2]=speed*std::sin(thetap)*std::cos(2*pi*rand1);

		//double speed2d=std::sqrt(std::pow(speed,2)-std::pow(speed1d,2));
		//double theta0 = std::atan(bounce_pos[1]/bounce_pos[0]);
		//double vp=speed2d*(2*rand0-1);  //parallel speed
		//double vr=-std::sqrt(std::pow(speed2d,2)-std::pow(vp,2)); //radial speed (inwards now - (minus sign)). 
		
        velocity[0]=(-vp*ybounce/radius+vr*xbounce/radius);
		velocity[1]=vr*ybounce/radius-vp*xbounce/radius; 
		
		//double phi=2*pi*rand1;
	
		//double theta=thetap+theta0;
		
		//std::cout<<"theta= "<<theta<<" theta0 = "<<theta0<<" thetap = " <<thetap <<"\n";
		//std::cout<<"position "<<bounce_pos[0]<<" "<<bounce_pos[1]<<" "<<bounce_pos[2]<<"\n";
		//velocity[0]=-speed2d*std::cos(theta);
		//velocity[1]=-speed2d*std::sin(theta);
		//velocity[2]=speed1d;
		//assert(sgn(velocity[0])!=sgn(bounce_pos[0]));
	//	std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
		//std::cin>> n;
	}

	
	

	position[0]=bounce_pos[0];
	position[1]=bounce_pos[1];
	position[2]=bounce_pos[2];
	dtwall=nextBounce(bounce_pos,velocity);

	
	return dtwall;

}
double cmsBoundary::nextBounce(double* position,double* velocity){
	if(element[current_element]==0){
			//int n;
		/*std::cout<<"boundary checking for cylinder! dt wall  "<<dtwall<<"\n";
		std::cout<<"position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
		std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
		std::cout<<"bstart "<<bstart[current_element]<<"\n";
		std::cout<<"benter "<<benter<<"\n";
		std::cout<<"bleave "<<bleave<<"\n";
		std::cout<<"current element "<<current_element<<"\n";
		*/



		//then do cylinder	
		//the following line is probably the biggest slow down in the code.
		//it is worth considering a less exact but faster solution. 
		double twalltemp1=(velocity[0]==0 && velocity[1]==0)  ? 1.0E32 :(-velocity[0]*position[0]-velocity[1]*position[1]\
				+std::sqrt(std::pow(sizes[3*current_element],2)*((std::pow(velocity[0],2)\
				+std::pow(velocity[1],2)))-std::pow(position[1]*velocity[0]-position[0]*velocity[1],2)))\
		    	/(std::pow(velocity[0],2)+std::pow(velocity[1],2));  //that took a while to get right...
   		
   		double tempz=position[2]-bstart[current_element]-sizes[3*current_element+2]/2.0;
  
   		double twalltemp2= velocity[2]==0 ? 1.0E32 :(sizes[3*current_element+2]/2.0-sgn(velocity[2])*tempz)/(sgn(velocity[2])*velocity[2]);
   		



   		//std::cout<<"circle wall time "<<twalltemp1<<" end wall time "<<twalltemp2<<"\n";
   		dtwall=twalltemp2<twalltemp1 ? twalltemp2:twalltemp1; 
   		whatwall=twalltemp2<twalltemp1 ? 2:4;
   			bounce_pos[0]=position[0]+velocity[0]*dtwall;
   			bounce_pos[1]=position[1]+velocity[1]*dtwall;
   			bounce_pos[2]=position[2]+velocity[2]*dtwall;
   			//normalization is necessary for this way of doing the trajectories... maybe slow... maybe do a margin like R.Schmid. 
   			if(whatwall==4){
   					bounce_pos[0]=bounce_pos[0]/std::sqrt(bounce_pos[0]*bounce_pos[0]+bounce_pos[1]*bounce_pos[1])*sizes[3*current_element];
   					bounce_pos[1]=bounce_pos[1]/std::sqrt(bounce_pos[0]*bounce_pos[0]+bounce_pos[1]*bounce_pos[1])*sizes[3*current_element];
   			} else{ 
   				bounce_pos[2]= (benter)?  bstart[current_element]: ( (bleave)? bstart[current_element+1] : (sgn(velocity[2])<0 ? bstart[current_element] : bstart[current_element]+sizes[3*current_element+2]));
   			}

			//		else std::cout<<"twall >0 whoa!!! = "<<dtwall<<" twall cylinder whatwall? "<<whatwall<<", z position "<<position[2]<<", z vel "<<velocity[2]<<"\n";
		
		benter=false;
		bleave=false;

 			//std::cout<<"dtwall cylinder "<<dtwall<<"\n";
		//if(dtwall<-1.0)std::cin>>n;
		return dtwall;
	} 
	else { 	//do rectangle
		int n;
		//This whole boudary code must drastically change if we want an offset axis for any elements. 
		//std::cout<<"boundary checking for rectangle! dt wall  "<<dtwall<<"\n";
		//std::cout<<"in position "<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
		//std::cout<<"velocity "<<velocity[0]<<" "<<velocity[1]<<" "<<velocity[2]<<"\n";
		//std::cout<<"bstart "<<bstart[current_element]<<"\n";
		//std::cout<<"benter "<<benter<<"\n";
		//std::cout<<"bleave "<<bleave<<"\n";
		//std::cout<<"current element "<<current_element<<"\n";
		
		double tempz=position[2]-bstart[current_element]-sizes[3*current_element+2]/2.0;
		//std::cout<<"tempz "<<tempz<<"\n";
		double twalltemp1=velocity[0]==0 ? 1.0E32 : (sizes[3*current_element+0]/2.0-sgn(velocity[0])*position[0])/(sgn(velocity[0])*velocity[0]);
		double twalltemp2=velocity[1]==0 ? 1.0E32 : (sizes[3*current_element+1]/2.0-sgn(velocity[1])*position[1])/(sgn(velocity[1])*velocity[1]);
   		double twalltemp3=twalltemp3=velocity[2]==0 ? 1.0E32 :(sizes[3*current_element+2]/2.0-sgn(velocity[2])*tempz)/(sgn(velocity[2])*velocity[2]);
   		
   		dtwall=(twalltemp1<twalltemp2) ?  ((twalltemp1<twalltemp3) ? twalltemp1:twalltemp3) : ((twalltemp2<twalltemp3) ? twalltemp2:twalltemp3); 
   		whatwall=(twalltemp1<twalltemp2) ?  ((twalltemp1<twalltemp3) ? 0:2) : ((twalltemp2<twalltemp3) ? 1:2); 
   		bounce_pos[0]=position[0]+velocity[0]*dtwall;
   		bounce_pos[1]=position[1]+velocity[1]*dtwall;
   		bounce_pos[2]=position[2]+velocity[2]*dtwall;
   		//std::cout<<"dtwall "<<dtwall<<" whatwall "<<whatwall<<"\n";
   		//normalization of bounce position is necessary
   		if (whatwall==0) 	  bounce_pos[0]= sgn(velocity[0])*sizes[3*current_element]/2.0;
   		else if (whatwall==1) bounce_pos[1]= sgn(velocity[1])*sizes[3*current_element+1]/2.0;
   		else   bounce_pos[2]= (benter)?  bstart[current_element]: ( (bleave)? bstart[current_element+1] : (velocity[2]<0 ? bstart[current_element] : bstart[current_element]+sizes[3*current_element+2]));
   		
   		benter=false;
   		bleave=false;
 		//std::cout<<"dtwall rect "<<dtwall<<"\n";
 		
		//std::cout<<"bounce pos "<<bounce_pos[0]<<" "<<bounce_pos[1]<<" "<<bounce_pos[2]<<"\n";
		//std::cout<<"dtwall "<<dtwall<<" whatwall "<<whatwall<<"\n";
		//if(std::isnan(position[0]) || std::isinf(position[0]))std::cin>>n;
   		return dtwall;
	}
	
	


}

		//randomly bunched up at the low z position 1 mm of the specified volume.
double* cmsBoundary::getInitPos(int elem, double randx,double randy, double randz){
	double* temp_pos= (double*)new double[3];//this is where the pointer for position is created.
	current_element=elem;
	if (element[elem]==0)
	{
		temp_pos[0]=(2*randx-1.0)*sizes[elem*3+0];
		temp_pos[1]=(2*randy-1.0)*sizes[elem*3+1];
		temp_pos[2]=randz*0.001+bstart[elem];//(double)sizes[elem*3+2]
	}
	else
	{
		temp_pos[0]=(randx-0.5)*sizes[elem*3+0];
		temp_pos[1]=(randy-0.5)*sizes[elem*3+1];
		temp_pos[2]=randz*0.001+bstart[elem];//(double)sizes[elem*3+2]	
	}
	return temp_pos;
}
