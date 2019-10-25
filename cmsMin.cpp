//g++ -I/data1/cmswank/BoostCodeSwank/boost_1_64_0 cmsMin.cpp -o mintest
#include "cmsMin.h"

std::vector<double> test(std::vector <double> x){
	
	for(std::vector<double>::iterator it = x.begin(); it != x.end(); ++it) {
    *it=82.2;
	}
	std::vector<double> answer=x;
	return answer;
}

std::vector<double> cmMin::df(std::vector<double> in){
	std::vector<double> inup;
	std::vector<double> indn;
	std::vector<double> temp1;
	std::vector<double> temp2;
	std::vector<double> temp;
	int ii = 0;
	for(std::vector<double>::iterator it = in.begin(); it != in.end(); ++it) {
    inup.push_back(*it+gstep.at(ii)*0.5);
    indn.push_back(*it-gstep.at(ii)*0.5);
	ii++;
	}
	temp1=functor(inup);
	temp2=functor(indn);
	ii = 0;
	for(std::vector<double>::iterator it = temp1.begin(); it != temp1.end(); ++it){
		temp.push_back(*it-temp2.at(ii)); 
		ii++;
	}
	return temp;
}

std::vector<double> cmMin::d2f(std::vector<double> in){
	std::vector<double> inup;
	std::vector<double> indn;
	std::vector<double> temp1;
	std::vector<double> temp2;
	std::vector<double> temp;
	int ii = 0;
	for(std::vector<double>::iterator it = in.begin(); it != in.end(); ++it) {
    inup.push_back(*it+gstep.at(ii)*0.5);
    indn.push_back(*it-gstep.at(ii)*0.5);
	ii++;
	}
	temp1=functor(inup);
	temp2=functor(indn);
	ii = 0;
	for(std::vector<double>::iterator it = temp1.begin(); it != temp1.end(); ++it){
		temp.push_back(*it-temp2.at(ii)); 
		ii++;
	}	


	return in;

}

int main(){
	double x = 2.0;
	double y = 2.0;

	//this is how you define a function address type. 
	//std::vector<double> (*pfun)(std::vector<double>);
	//pfun=test;
	//std::cout<<&p_fun;
	
	std::vector<double> tester;
	std::vector<double> stp;
	std::vector<double> gtp;

	for(int i = 0; i<10; i++)
	{
		tester.push_back(82.0);
		gtp.push_back(0.1);	
		stp.push_back(1.0);
	}
	cmMin *cmin;
	cmin=new cmMin(test,stp,gtp);
	std::vector<double> farts;
	std::cout<<cmin->functor<<"\n";
	farts=cmin->df(tester);
	double temp;
	temp=farts[1];
	std::cout<<temp<<"\n";
	return 0;

}
