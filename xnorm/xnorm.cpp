//*****************************************************************************
//                                                         
//		Xnorm							 
//														 
//		v. 4.3 - 20151021
//
//		2007-2015 - Nicola Ferralis
//                                                        
//		Renormalize X axis from Cartesian data in ASCII format		
//
//		This program (source code and binaries) is free software; 
//		you can redistribute it and/or modify it under the terms of the
//		GNU General Public License as published by the Free Software 
//		Foundation, in version 3 of the License.
//
//		This program is distributed in the hope that it will be useful,
//		but WITHOUT ANY WARRANTY; without even the implied warranty of
//		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//		GNU General Public License for more details.
//
//		You can find a complete copy of the GNU General Public License at:
//		http://www.gnu.org/licenses/gpl.txt
//												             
//**********************************************************************************


#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <iomanip>
using namespace std;


#define OPENMP

int operate(char *namein);
int ReadKey();
float ReadKeyF();
int decimal(double a, int maxDec);

int maxDec = 6;
int pX = 1; //decimal precision X-axis
int pY = 5; //decimal precision Y-axis

char version[]="4.3 - 20151021";
char extension[]="xnorm.";

int main(int argc, char *argv[])
{
    
    if(argc<2)
	{cout<< "\n xnorm v."<<version<< "\n\n";
        cout<< " Usage:\n xnorm <filename> \n\n";
        }
    
    if(argc>=2)
	{  operate(argv[1]);
        
        return 0;
    }
}


//OPERATE 
int operate(char *namein)

{	double Egrid;
	unsigned long i,j,npoints, npoints2;

	
   	ifstream infile(namein);
	
	if(!infile)
        {
        cout<<"\n file '"<< namein<<"' not found\n";
        return 0;
        }
	string tmp;
    unsigned long ind = 0;
    while (!infile.eof()) {
        
        getline(infile, tmp);
        ind++;
    }
    infile.close();
    
////
    
    
    double * energy1 = new double[ind];
    double * intensity1 = new double[ind];
    
    double enIn, enEnd = 0.0;
    
    j=0;
	
	string a, b;

   	infile.open(namein);
	while(!infile.eof()| !infile.eof())
		{infile>>a>>b;
		
		replace(a.begin(),a.end(), ',','.');
		replace(b.begin(),b.end(), ',','.');
			
		energy1[j] = atof(a.c_str());
		intensity1[j] = atof(b.c_str());
		j++;}
    
    infile.close();

    cout<<"\nPress 1) fixed x step; 2) for fixed number of data points:   ";
	int type2;
	type2=ReadKey();
    
    
	npoints = ind;
	npoints2 = 0;
	Egrid = 1;
    
	enIn=fabs(energy1[0]);
	enEnd = fabs(energy1[npoints-1]);
    
	if (type2==2)
		{
		cout<<"\n number of points: ";
        npoints2=ReadKey();
        Egrid=(enEnd-enIn)/npoints2;
            cout<<Egrid<<"\n";
            //if(decimal(Egrid, 6)>=6)
              //  {cout<< "\n Precision too high. Reduce number of points and try again. \n\n";
              //  return 0;}
        }

	if (type2==1)
        {
        cout<<"\n x step:  ";
        Egrid= ReadKeyF();
            if(decimal(Egrid, 6)>=6)
                {cout<< "\n Precision too high. Increase x step and try again. \n\n";
                return 0;}
            
        npoints2=(long) ((enEnd-enIn)/Egrid)+1;}
    
    double * energy = new double[npoints2];
    double * intensity = new double[npoints2];
    double * intensity2 = new double[npoints2];
    
    
#ifdef OPENMP
#pragma omp parallel for shared(Egrid)
#endif

	for(unsigned long f=0;f<npoints2; f++)
		{energy[f]=enIn+f*Egrid;
		intensity[f]=0;
		}
    
#ifdef OPENMP
#pragma omp parallel for
#endif

    
	for (i=0;i<npoints2;i++)
		{j=0;
		while (true)
            {j++;
                if ((energy1[j]>= energy[i])|(j>=npoints))
                    {break;}
        }
            
	j=j-1;
            
    double I0, I1, x0, x1,x;
	I0=intensity1[j];
	I1=intensity1[j+1];
	x0=energy1[j];
	x1=energy1[j+1];
	x=energy[i];
	intensity[i]=(I0+((I1-I0)/(x1-x0))*(x-x0));
	}
    
    
    double Idiff = 100.0;
    double min = 0.0;
    
    if(intensity[0]<min)
        {min=intensity[0];}
    
#ifdef OPENMP
#pragma omp parallel for shared(Idiff)
#endif
    for (i=1;i<npoints2;i++)
        {if ((intensity[i]-intensity[i-1]) < Idiff)
            {Idiff=intensity[i]-intensity[i-1];}
            
        if(intensity[i]<min)
            {min=intensity[i];}
        }
    
    
//******************************************************
//Y value renormalization (for negative values)
    
    if(min>=0)
        {cout<<"\n No negative values detected (min Y = "<<min<<")\n";}
    
    else
        {int type3;
        cout<<"\n Min Y = "<<min<<" - Renormalize negative values? (1:NO - 2:YES):  ";
        type3 = ReadKey();
            
            if (type3==2) {
                
            cout<<"\n Correcting negative values\n";
            #ifdef OPENMP
            #pragma omp parallel for shared(min)
            #endif
                for (unsigned long h=0; h<npoints2; h++)
                {intensity[h] += fabs(min);}
                }
            
        }
    

    delete[] energy1;
    delete[] intensity1;
    
    char* outname = new char[strlen(namein)+strlen(extension)+1];
    sprintf(outname, "%.*s%s", (int)  strlen(namein)+5, extension, namein);
	
    ofstream outfile(outname);


    double precGrid = decimal(Egrid, maxDec);
    double precInt = decimal(Idiff,maxDec);

	for(i=0; i<npoints2; i++)
		{outfile<<fixed<<setprecision(precGrid+pX)<<energy[i]<<"\t";
        outfile<<fixed<<setprecision(precInt+pY)<<intensity[i]<<"\n";
        }

	delete[] energy;
    delete[] intensity;
    delete[] intensity2;
    
	outfile.close();
	cout<<"\n Saved in: "<<outname<<"\n";
    delete[] outname;
	return 0;}

int ReadKey()
{	string tkc;
	int tk;
	cin>>tkc;
		
    tk=(int) atof(tkc.c_str());
	if(tk<0)
		{return 10;}
	else
		{}

	return tk;
}

float ReadKeyF()
{	string tkc;
	float tk=0.0;
	cin>>tkc;
	tk=atof(tkc.c_str());
    
	return tk;
}

int decimal(double a, int maxDec) {
    a -= round(fabs(a));
    int t=0;
    double b=0;
    
    while(a!= 0 || a>1e-10)    {
        if(fabs(a/b)<0.1)
            {//cout<<"break: ";
               // cout<<t<<"\t"<<fixed<<setprecision(20)<<a<<"\t"<<b<<"\t"<<a/b<<endl;
                break;}
        else
            {//cout<<t<<"\t"<<fixed<<setprecision(20)<<a<<"\t"<<b<<"\t"<<a/b<<endl;
            b=a;
            a*=10;
            a = a - round(a);
            t++;}
    }
    
    if (t>maxDec)
        {t=maxDec;}
    return t;
}
