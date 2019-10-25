//g++ -o extractrootbinary TimeSxSySzBinBinaryArgumentsROOT534.cpp `root-config --cflags --glibs`
#include <fstream>
#include <iostream>
#include <sstream>
//#include <stdlib.h>
#include <stdio.h>
//#include <unistd.h> 
//#include <stdlib.h>
//#include <math.h>
//#include <cmath> 
#include "TFile.h"
#include "TTree.h"
//#include "TFormula.h"
//#define mem 1000000
using namespace std;

//#define Bin 1001 //this should add one to the par file. 
//#define Event 1000 

 	
int main(int argc, char *argv[] )
{	
	
	
	if(argc!=2) cout<<"Only input argument is run #";

	string run=argv[1];	
	
	//usefull code for conversions	
	/*	
	string s2=argv[2];
    
    stringstream convert(s2);//object from the class stringstream.
    convert>>Bin; //the object has the value 12345 and stream it to the integer x.
    //now the variable x holds the value 12345.

	Bin=Bin+1;
	*/

	
	
    const double pi = 3.1415926535897;  //pi
   

	//get filename from run number
	stringstream ss1;	
	ss1<<"/home/ezra/spin-sim-SD/ExtractData/run"<<run<<".root";
	

	string filename1=ss1.str();  
	  

	//open root file
	TFile *myFile = new TFile(filename1.c_str(), "UPDATE");	    
	TTree *dat; 
	myFile->GetObject("dat", dat);;

	//get needed run parameters     
	Int_t linum=dat->GetEntries();//The number of events?
	
	int Bin;	
	dat->SetBranchAddress("n",&Bin);
	dat->GetEntry(0);	 	
	
	//init tree variables
    Double_t Sx[Bin];
    Double_t Sy[Bin];     
    Double_t Sz[Bin];
    Double_t t[Bin];
	Double_t X[Bin];
    Double_t Y[Bin];     
    Double_t Z[Bin];
    Double_t VX[Bin];
    Double_t VY[Bin];
    Double_t VZ[Bin];


   
   	//init export variables
	double SSX; 
    double SSY;
	double SSZ;
	double XX;
	double YY;
	double ZZ;
	double VVX;
	double VVY;
	double VVZ;
	double TT;
	
	//create file with a name
	stringstream ss2;
	ss2<<"/home/ezra/spin-sim-SD/ExtractData/OutputHeatflush"<<run<<".dat";
	string filename2=ss2.str();
	ofstream  myfile;
   
	//format file
 	myfile.open (filename2.c_str(),ios::binary | ios::out);
	
	//set branches	
    dat->SetBranchAddress("t",&t);
    dat->SetBranchAddress("sx",&Sx);
    dat->SetBranchAddress("sy",&Sy);
	dat->SetBranchAddress("sz",&Sz);
	dat->SetBranchAddress("x",&X);
    dat->SetBranchAddress("y",&Y);
	dat->SetBranchAddress("z",&Z);
	dat->SetBranchAddress("vx",&VX);
	dat->SetBranchAddress("vy",&VY);
	dat->SetBranchAddress("vz",&VZ);    


	for(int ii=0;ii<linum;ii++)
    {
      
		
      
    
		//fill double array with tree branches
		dat->GetEntry(ii);
 		
		for(int j=0;j<Bin;j++){

		//double from doulbe array entry
		TT=t[j];
		SSX=Sx[j];
		SSY=Sy[j];
		SSZ=Sz[j];
		XX=X[j];
		YY=Y[j];
		ZZ=Z[j];
		VVX=VX[j];
		VVY=VY[j];
		VVZ=VZ[j];

		//Export double
		myfile.write(reinterpret_cast <const char*> (&SSX),sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&SSY),sizeof(double)); 		
		myfile.write(reinterpret_cast <const char*> (&SSZ),sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&XX), sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&YY), sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&ZZ), sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&VVX), sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&VVY), sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&VVZ), sizeof(double)); 
		myfile.write(reinterpret_cast <const char*> (&TT), sizeof(double)); 	


     }
    } 
  
	myfile.close();    
	cout<<"Bin #: "<<Bin<<"   Event #: "<<linum<<"\n";
	return 0;
}
