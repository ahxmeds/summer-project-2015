#include<iostream>
#include<stdio.h>
#include<time.h>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<string>
#include<sys/stat.h>


using namespace std;

//structure for Lattice type
struct LATTICE
	{
	  int x;
	  int y;
	};
	

//declaring and initializing global variables 
const int L = 60; //length of the lattice
const int N = L*L; //number of spins
const int TempCounts = 51; //number of discrete temperatures used

double Lattice[L][L]; //2D array to store the lattice
float V[TempCounts][16]; //2D array to store the number of congfigurations at each temperature
float C[TempCounts][2]; //2D array to store specific heat for the purpose of calculating entropy
float Tri[TempCounts][2];//2D array to store the number of triangles with 3 unhappy and 1 unhappy bonds
float T;

int monte_carlo_steps = 10000; //number of monte carlo steps
int transient_steps = 50000; //number of steps to achieve equilibrium
float T_begin = 5.0; //value of temperature at which we begin simulation
float T_steps = 0.1; //value by which we decrement the temperature at each iteration
long int seed = 2147483647; //seed value for random number generator

double J1 = 0.89; //coupling constant for diagonal bonds
double J2 = 1.0; //coupling constant for nearest neighbour interactions, always set to 1
double p = floor(float(J1/J2)*100 + 0.5)/100; //value of J1/J2 to two decimal places


//function to print the lattice, if need be
void print_lattice()
	{ 
	  
	  for(int i=0; i<L; i++)
	   {
	     for(int j=0; j<L; j++)
	     	{
	     	   if(Lattice[i][j]==1)
	     	   cout<<"+"<<" ";
	     	   else
	     	   cout<<"-"<<" ";
	     	}
	     	cout<<"\n";
	   }
	 cout<<"\n";  
	}


//function to initialize the lattice randomly using an almost uniform random number generator
void initialize_lattice(double Lattice[L][L])
	{  
	  int i,j;
	  double RN;
	  for(i=0; i<L; i++)
	  	{
	  	  for(j=0; j<L; j++)
	  	    {
	  	      RN = (double)rand()/(double)RAND_MAX;
	  	      if(RN<0.5)
	  	      	Lattice[i][j] = 1;
	  	      if(RN>=0.5)
	  	      	Lattice[i][j] = -1;
	  	    }
	  	}
  }


//modulus function to return positive values always
int mod(int a, int b)
	{
	   int r;
	   r = a%b;
	   if(r<0)
	      r = r+b;
	   return r;   
	}
	
	
//function to choose a point on the lattice randomly
void choose_random_position(LATTICE &pos)
	{  
	   pos.x= mod(rand(),L);
	   pos.y= mod(rand(),L);
	   
	}


//function to return the 2xenergy at point (x,y) 
double Energy_at_xy(LATTICE &position)
	{  
	   double Energy;
	   Energy = Lattice[position.x][position.y]*(J2*(Lattice[position.x][mod(position.y+1,L)]+Lattice[position.x][mod(position.y-1, L)] + Lattice[mod(position.x-1,L)][position.y] + Lattice[mod(position.x+1,L)][position.y]) + J1*(Lattice[mod(position.x-1,L)][mod(position.y-1,L)] + Lattice[mod(position.x+1, L)][mod(position.y+1, L)]));
	 
	  return Energy;
	}	
	
	
//function to check whether flip is valid or not	
bool test_flip(LATTICE position, double &Del_E)
	{
	  Del_E = -2*Energy_at_xy(position);
	  
	  if(Del_E<0)
	  	return true;
	  
	  else if ((float)rand()/(float)((unsigned)RAND_MAX) < exp(-(float)Del_E/T))
	    return true;
	   
	  else 
	    return false;
	}	
	
	
//function to flip the spin at site (x,y)
void flip_at_xy(LATTICE position)
	{
	  Lattice[position.x][position.y] = -Lattice[position.x][position.y];
	}	
	
	
//function to attain equilibrium at each temperature for 'transient steps' number of iterations
void transient(void)
	{
	  LATTICE position;
	  double de = 0, t,i;
	  
	  for(t=1; t<=transient_steps; t++)
	  	{
	  	  for (i=1; i<=N; i++)
	  	  {
	  	   choose_random_position(position);
	  	   
	  	   if(test_flip(position, de))
	  	   {
	  	  	flip_at_xy(position);
	  	 	 }
	  	 	}
	  	}
	}
	

//function to initialize the array V
void initialize_V_to_zero(float V[TempCounts][16])
	{
	  int i,j;
	  
	  for(i=0; i<TempCounts; i++)
	  	{
	  	  for(j=0; j<16; j++)
	  	   V[i][j] = 0;
	  	}
	}
	
	
//function to initialize the array Tri
void initialize_Tri_to_zero(float Tri[TempCounts][2])
	{
	  int i,j;
	  
	  for(i=1; i<TempCounts; i++)
	  	{
	  	  Tri[i][0] = 0;
	  	  Tri[i][1] = 0;
	    }
	}	


//function to return the total magnetisation of the system
double Sum_Magnetisation()
	{
	  int i,j;
	  double M=0.0;
	  
	  for(i=0; i<L; i++)
	   {
	    for(j=0; j<L; j++)
	     {
	      M = M + pow(1, abs(i-j))*Lattice[i][j];
	     }
	   }
	  
	  return M; 
	}
	

//function to return 2x(total energy) of the system
double Sum_Energy()
	{
	  LATTICE position;
	  int x,y;
	  double E=0.0; 	  
	  for(x=0; x<L; x++)
	  	{  
	  	  position.x = x;
	  	  for(y=0; y<L; y++)
	  	   { 
	  	     position.y = y;
	  	     E = E + Energy_at_xy(position);
	  	   }
	  	}
	  	
	  return E;
	}	
	

//function to calculate and store the number of squares with various configuration at different temperatures
void HistogramVector()
	{ 
	  int i,j;
	  int t = round(T*int(1/T_steps));
	  
	  double Lat[L][L];
	  
	  //copying lattice to lat
	  for(i=0; i<L; i++)
	  {
	   for(j=0; j<L; j++)
	   	{
	   	  Lat[i][j] = Lattice[i][j];
	   	}
	  }
	  
	  
	  for(i=0; i<L; i++)
	    {
	     for(j=0; j<L; j++)
	      { 
	        if((i+j)%2==0) //if squares are unshaded
	        {
	         if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][0]++;
	         else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][1]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][2]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][3]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][4]++;
	         else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][5]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][6]++;
	       	 else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][7]++;
	       	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][8]++;
	      	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][9]++;
	      	 else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][10]++;
	       	 else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][11]++;
	       	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][12]++;
        	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][13]++;
		     	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][14]++;
	         else 
	        	V[t][15]++;	
	        }		 
	        
	        else //if squares are shaded
	        {
	          Lat[i][j] = -1*Lat[i][j];
	          Lat[i][mod(j+1,L)] = -1*Lat[i][mod(j+1,L)];
	          Lat[mod(i+1,L)][j] = -1*Lat[mod(i+1,L)][j];
	          Lat[mod(i+1,L)][mod(j+1,L)] = -1*Lat[mod(i+1,L)][mod(j+1,L)];
 	         
	          if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][0]++;
	         else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][1]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][2]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][3]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][4]++;
	         else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][5]++;
	         else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][6]++;
	       	 else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][7]++;
	       	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][8]++;
	      	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][9]++;
	      	 else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][10]++;
	       	 else if(Lat[i][j]==-1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][11]++;
	       	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==-1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][12]++;
        	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==1 && Lat[mod(i+1,L)][mod(j+1,L)]==-1)
	        	V[t][13]++;
		     	 else if(Lat[i][j]==1 && Lat[i][mod(j+1,L)]==1 && Lat[mod(i+1,L)][j]==-1 && Lat[mod(i+1,L)][mod(j+1,L)]==1)
	        	V[t][14]++;
	         else 
	        	V[t][15]++;
	        	
	          Lat[i][j] = -1*Lat[i][j];
	          Lat[i][mod(j+1,L)] = -1*Lat[i][mod(j+1,L)];
	          Lat[mod(i+1,L)][j] = -1*Lat[mod(i+1,L)][j];
	          Lat[mod(i+1,L)][mod(j+1,L)] = -1*Lat[mod(i+1,L)][mod(j+1,L)]; 
 	         }       	
	      }
	   } 

	}


//function to calculate and store the value of number of triangles with 3 unhappy/1 unhappy bond(s) at different temperatures
void HistoTri()
	{ 
	   int i,j;
	   
	   int t = round(T*int(1/T_steps));
	   
	   for(i=0; i<L; i++)
	    {
	      for(j=0; j<L; j++)
	      {
	        if(Lattice[i][j] == Lattice[i][mod(j+1,L)] && Lattice[i][mod(j+1,L)] == Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][0]++;
	        if(Lattice[i][j] == Lattice[i][mod(j+1,L)] && Lattice[i][mod(j+1,L)] == -1*Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][1]++;
	        if(Lattice[i][j] == -1*Lattice[i][mod(j+1,L)] && Lattice[i][mod(j+1,L)] == Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][1]++;
	        if(Lattice[i][j] == -1*Lattice[i][mod(j+1,L)] && Lattice[i][mod(j+1,L)] == -1*Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][1]++;
	        
	        if(Lattice[i][j] == Lattice[mod(i+1,L)][j] && Lattice[mod(i+1,L)][j] == Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][0]++;
	        if(Lattice[i][j] == Lattice[mod(i+1,L)][j] && Lattice[mod(i+1,L)][j] == -1*Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][1]++;
	        if(Lattice[i][j] == -1*Lattice[mod(i+1,L)][j] && Lattice[mod(i+1,L)][j] == Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][1]++;
	        if(Lattice[i][j] == -1*Lattice[mod(i+1,L)][j] && Lattice[mod(i+1,L)][j] == -1*Lattice[mod(i+1,L)][mod(j+1,L)])
	         Tri[t][1]++;
	     }
	   }
	}


//driver function
int main()

{  
   //seeding the srand function
   srand(seed);
   
   //declaring the file pointers and variables to store file names 
   ofstream file1, file2, file3, file4, file5, file6, file7, file8, file9, file10;
   string fileName1,fileName2,fileName3,fileName4,fileName5,fileName6,fileName7,fileName8,fileName9, fileName10;
   
   //intializing the file names 
    fileName1 = "E_L="+ to_string(L) +"_p="+ to_string(p)+".csv";  
    fileName2 = "E2_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName3 = "M_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName4 = "Mabs_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName5 = "Cv_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName6 = "Chi_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName7 = "BiCu_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName8 = "Hist_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName9 = "S_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    fileName10 = "Tri_L="+ to_string(L) +"_p="+ to_string(p)+".csv";
    
    //opening files for storing various data
    file1.open(fileName1.c_str(), ios::out);
    file2.open(fileName2.c_str(), ios::out);
    file3.open(fileName3.c_str(), ios::out); 
    file4.open(fileName4.c_str(), ios::out);
    file5.open(fileName5.c_str(), ios::out);
    file6.open(fileName6.c_str(), ios::out);
    file7.open(fileName7.c_str(), ios::out);
    file8.open(fileName8.c_str(), ios::out);
    file9.open(fileName9.c_str(), ios::out);
    file10.open(fileName10.c_str(), ios::out);
   
   //declaring and initializing various variables for storing the value of observables
   double Energy=0, Energy_avg=0, Energy_sq_avg=0, Energy_sum=0, Energy_sq_sum=0;
   double Mag=0, Mag_avg=0, Mag_sq_avg=0, Mag_sum=0, Mag_sq_sum=0;
   double Mag_abs=0, Mag_abs_avg=0, Mag_q_avg=0, Mag_abs_sum=0, Mag_q_sum=0;
   double Specific_Heat, Mag_Susceptibility, Binder_Cumulant, BiCuE,Correl;
   double Entropy;
   
   int t;
   double de=0;
   int i,j;
   LATTICE position;
   
   initialize_Tri_to_zero(Tri);
   initialize_V_to_zero(V);
   initialize_lattice(Lattice);
   
   //this will write the first line to the file10, which stores the number of triangles with 3 unhappy/1 unhappy bonds
   file10<<"T"<<"\t"<<"No. of tri with 3 UH bonds"<<"\t"<<"No. of tri with 1 UH bond"<<"\n";
   
   //Temperature loop
   for(T=T_begin; T>0.001; T = (round(int(1/T_steps)*(T - T_steps))/int(1/T_steps)))
   	{ 
   	  cout<<"Beginning Temp="<<T<<"\n"; 
   	  
   	  transient();
   	  
   	  Energy = Sum_Energy();
   	  Mag = Sum_Magnetisation();
   	  Mag_abs = fabs(Mag);
   	  
   	  //initializing the summation variables of each observable to zero
   	  Energy_sum=0;
   	  Energy_sq_sum=0;
   	  Mag_sum=0; 
   	  Mag_sq_sum=0;
   	  Mag_abs_sum=0;
   	  Mag_q_sum=0;
   	  
   	  //start of MonteCarlo loop
   	  for(i=1; i<=monte_carlo_steps; i++)
   	  	{
   	  	   //start of Metropolis loop
   	  	   for(j=1; j<=N; j++)
   	  	    {
   	  	      choose_random_position(position);
   	  	     
   	  	      if(test_flip(position,de))
   	  	      	{
   	  	      	  flip_at_xy(position);
   	  	      	  
   	  	      	  //updating the value of observables
   	  	      	  Energy = Energy + de;
   	  	      	  Mag = Mag + 2*pow(1, abs(position.x-position.y))*Lattice[position.x][position.y];
   	  	      	  Mag_abs = Mag_abs + fabs(Lattice[position.x][position.y]);
   	  	      	}
   	  	    }
   	        
   	        //adding the observables to the summation variables
   	        Energy_sum += Energy/2.0;
   	        Energy_sq_sum += (Energy/2.0)*(Energy/2.0);
   	        Mag_sum += Mag; 
      	        Mag_sq_sum += Mag*Mag;
   	        Mag_q_sum += Mag*Mag*Mag*Mag;  	        
   	        Mag_abs_sum += sqrt(Mag*Mag);

		
		if(i%100==0)
		HistogramVector();
   	       	        
   	    }

     //averaging the observables per unit spin
     Energy_avg = Energy_sum/((float)monte_carlo_steps*N);
     Energy_sq_avg = Energy_sq_sum/((float)monte_carlo_steps*N);
     Mag_avg = Mag_sum/((float)monte_carlo_steps*N);
     Mag_sq_avg = Mag_sq_sum/((float)monte_carlo_steps*N);
     Mag_abs_avg = Mag_abs_sum/((float)monte_carlo_steps*N);
     Mag_q_avg = Mag_q_sum/((float)monte_carlo_steps*N);
     
     //calling functions that are used for storing the values to plot histograms
    // HistogramVector();
     HistoTri();
     
     //writing to file8, which stores the number of squares with various configurations
     file8<<"Config"<<"\t"<<"T = "<<T<<"\n";
     
     t = round(T*int(1/T_steps));
     
     //writing the values to file8 again, the number of squares in the lattice pertaining to particular configurations, at a particular temp
     for(int k=0; k<=15; k++)
     	{
     	 file8<<k<<"\t"<<V[t][k]<<"\n";
     	}  
     file8<<"\n";
     
     //writing the values to file10, whose first column stores the value of temperature, T, second column stores number of triangles with 
     //3 unhappy bonds and third column stores number of triangles with  1 unhappy bond
     file10<<T<<"\t"<<Tri[t][0]<<"\t"<<Tri[t][1]<<"\n";
     
     //calculating other observables based of the values obtained so far
     Specific_Heat = (Energy_sq_avg - Energy_avg*Energy_avg*N)/(T*T);
     Mag_Susceptibility = (Mag_sq_avg - Mag_avg*Mag_avg*N)/(T);
     Binder_Cumulant = 1.0 - (Mag_q_avg)/(3.0*Mag_sq_avg*Mag_sq_avg*N);
     
     //putting values into the C array, which will be used latter to calculate entropy
     C[t][0] = T;
     C[t][1] = Specific_Heat;
     
     //writing data to various files
     file1<<T<<"\t"<<Energy_avg<<"\n";
	 file2<<T<<"\t"<<Energy_sq_avg<<"\n";
     file3<<T<<"\t"<<Mag_avg<<"\n"; 
     file4<<T<<"\t"<<Mag_abs_avg<<"\n";
	 file5<<T<<"\t"<<Specific_Heat<<"\n";
	 file6<<T<<"\t"<<Mag_Susceptibility<<"\n";
	 file7<<T<<"\t"<<Binder_Cumulant<<"\n";
	  	
     //printing a flag on the screen after the run for specific T is over
     cout<<"L="<<L<<"\t"<<"p="<<p<<"\t"<<"Temp = "<<T<<" is done\n";
     
     //function call to print the lattice obtained at the end of each iteration. This also becomes the input lattice for the next interation.
     //comment out this function call if it takes too much time while running the code
     print_lattice();
  
  }
   
   //calculating and writing data to file9, which stores the entropy at different temperatures.
   //the array C has been used for the purpose of calculation
   for(i=1; i<=TempCounts-1; i++)
   	{ 
   	  Entropy = 0.0;
   	  
   	  file9<<C[i][0]<<"\t";
   	  
   	  for(j=i; j<=TempCounts-1; j++)
   	  {
   	    Entropy += T_steps*(double((C[j][1])/(C[j][0]))); 
   	  }
   	  
   	  file9<<log(2.0)-Entropy<<"\n";
   	}
    
   //closing the files
   file1.close();
   file2.close();
   file3.close();
   file4.close();
   file5.close();
   file6.close();   
   file7.close();
   file8.close();
   file9.close();
   file10.close();
   
 return 0;
   
}


	
	
		
		
