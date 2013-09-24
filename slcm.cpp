#include<iostream>
#include<cmath>
#include<stdio.h>
#include<cassert>
#include<fstream>
using namespace std;

struct loop_output {
  double* c;
  double* D;
  double s;
  int* mass_distribution;
} output, Simultaneous_Loop, out;

double randomdouble (double min, double max);
int  mod (int a, int L);
int*  add (int *u, int *v, int L);
int* Initial (double rho, int L);
int* Sim_Dyn (int* u2, int* vinn, double omega, int L, double dt);
int* Initial_Run (int* vin, int initial_max, double omega, int L, double dt);
int* MinusMass (int* vector, int sub, int L);
double* cx (int* vector, int L, double* clocal);
double sx (int* vector, int L);
double* Dx (int* vector, int L, double* DLocal);
loop_output Simulation_Loop (double rho, int* vin, int kmax, double omega, int L, double dt, double rhostep, loop_output output);
double** cDs (double** vector, int L, double rhostep, double rhoint, int imax);
double*** data (int initial_max, int imax, int kmax, double rho, int L, double omega, double dt, double rhostep, loop_output output, int trial_count);
double*  dadd (double *u, double* v, int L);


int main (int argc, char * argv[])
{
    srand ((double) time(NULL));
    int L;
    double omega;
    const double dt = 0.1;
    // here the starting value of rho is defined
    double rho = 4.3;
    int kmax = 100000;
    //rhostep gives the change in rho and imax will let you state how many values of rho are used in a single simulation.
    double rhostep = 1.0;
    int initial_max = 200000;
    int imax = 4;
    int trial_count = 20;
    char buffer [50];
    int n;
       //  Here we define an array fin which data will be written to.  sprintf is then used so that when fin is written to a saved document, the name will depend on our input variables L and omega.  Below you want to define the range of values of L and omega to be used in the simulation.  Once these are chosen, compile through terminal and run (no inputs are needed).
    double*** fin = new double** [imax];
    for (int i = 0; i<imax; i++)
      {
	fin[i] = new double* [4];
      }
    for (int i = 0; i < imax; i++)
      {
	for (int j = 0; j < 4; j++)
	  {
	    int length = L + 1;
	    fin[i][j] = new double [length];
	  }
      }
    
     for (omega = 2; omega < 2.5; omega += 0.5)
       { 
	for (L = 100; L < 102; L += 50)
	{ 
	  fin = data (initial_max, imax, kmax, rho, L, omega, dt, rhostep, output, trial_count);
	   n = sprintf (buffer,"8_8_Data_L_%d_omega_%4.2f.tsv", L, omega);
	  std::ofstream write_output (buffer);
	  assert (write_output.is_open ());
	  for (int j = 0; j<imax; j++)
	    {
	      for (int k = 0; k < L+1; k++)
		{
		  write_output << fin[j][0][k] << "\t" << fin[j][1][k]  << "\t" << fin[j][2][k] << "\t" << fin[j][3][k] << "\n";
		}
	    }
	  write_output.close();
	  cout << "{omega, L} =" << "  {" << omega << ", " << L << "} done" << "\n";
      	}
       }
      for (int i = 0; i < imax; i++)
       {
	 for (int k = 0; k < 4; k++)
	   {
	     delete[] fin[i][k];
	   }
       }
    for (int i = 0; i < imax; i++)
      {
	delete[] fin[i];
	}
    delete[] fin;
    return 0;
}

/*
data runs the simulation over many trials and averages s, C, and D over the trials and generates the standard deviation for each value.  The final data is stored in a triple pointer. This can then be written to a text file to be imported into mathematica.  The way the data is stored is a little funky, so a little bit of manipulation is required to get all of the data separated again.
 */

double*** data (int initial_max, int imax, int kmax, double rho, int L, double omega, double dt, double rhostep, loop_output output, int trial_count)
{
    double* sumrho = new double [imax];
    double** sumc = new double* [imax];
    double** sumD = new double* [imax];
    double* sums = new double [imax];
    for (int i = 0; i < imax; i++)
      {
	sumc[i] = new double [L];
      }
    
    for (int i = 0; i < imax; i++)
      {
	sumD[i] = new double [L];
      }


    for (int i = 0; i<imax;i++)
      {
      for (int j = 0; j < L; j++)
	{
	  sumc[i][j]=0;
	}
      }
    
    for (int i = 0; i<imax;i++) 
      {
	for (int j = 0; j < L; j++)
	  {
	    sumD[i][j]=0;
	  }
      }
    
      for (int i = 0; i<imax;i++) {
        sums[i]=0;
      }

    for (int i = 0; i<imax;i++) {
        sumrho[i]=0;
    }



    double** stand_dev_squarec = new double* [imax];
    // double** stand_devc = new double* [imax];
    double** stand_dev_squareD = new double* [imax];
    // double** stand_devD = new double* [imax];
    double* stand_dev_squares = new double [imax];
    // double* stand_devs = new double [imax];
    for (int i = 0; i < imax; i++)
      {
	stand_dev_squarec[i] = new double [L];
      }
   
     for (int i = 0; i < imax; i++)
      {
	stand_dev_squareD[i] = new double [L];
      }

    for (int i = 0; i<imax;i++) 
      {
	for (int j = 0; j < L; j++)
	  {
	    stand_dev_squarec[i][j]=0;
	  }
      }
    
    
    for (int i = 0; i<imax;i++)
      {
	for (int j = 0; j < L; j++)
	  {
	    stand_dev_squareD[i][j]=0;
	  }
      }
    

    for (int i = 0; i<imax;i++) {
      stand_dev_squares[i]=0;
    }
    
    int* vin = new int [L];
    double*** ave_cD = new double** [imax];
    for (int i = 0; i < imax; i++)
      {
	ave_cD[i] = new double* [2];
      }
    for (int i = 0; i < imax; i++)
      {
	for (int j = 0; j < 2; j++)
	  {
	    ave_cD[i][j] = new double [L];
	  }
      }
    double** ave_rs = new double* [imax];
    for (int i = 0; i < imax; i++)
      {
	ave_rs[i] = new double [2];
      }
    double rhoint = rho;
    for (int i = 0; i < trial_count; i++)
      {
	vin = Initial_Run (Initial (rho, L), initial_max, omega, L, dt);
	for ( int j = 0; j< imax; j++)
	  {
	    cout << j << "\n";
	    out = Simulation_Loop (rho, vin, kmax, omega, L, dt, rhostep, output);
	    ave_cD[j][0] = out.c;
	    ave_cD[j][1] = out.D;
	    ave_rs[j][0] = rhoint - rhostep *(j+1);
	    ave_rs[j][1] = out.s;
	    vin = out.mass_distribution;
	    sumrho[j] += ave_rs[j][0]/double(trial_count);
	    sums[j] += ave_rs[j][1]/double(trial_count);
	    for (int k = 0; k < L; k++)
	      {
		sumc[j][k] += ave_cD[j][0][k]/double(trial_count);
		sumD[j][k] += ave_cD[j][1][k]/double(trial_count);
		stand_dev_squarec[j][k] += (1/double(trial_count))*(ave_cD[j][0][k] - sumc[j][k])*(ave_cD[j][0][k] - sumc[j][k]);
		stand_dev_squareD[j][k] += (1/double(trial_count))*(ave_cD[j][1][k] - sumD[j][k])*(ave_cD[j][1][k] - sumD[j][k]);
	      }
	    stand_dev_squares[j] += (1/double(trial_count))*(ave_rs[j][1] - sums[j])*(ave_rs[j][1] - sums[j]);
	  }
      }
    
    
    for (int i = 0; i<imax; i++)
      {
	for (int k = 0; k < L; k++)
	  {
	    stand_dev_squarec[i][k] = sqrt (stand_dev_squarec[i][k]);
	    stand_dev_squareD[i][k] = sqrt (stand_dev_squareD[i][k]);
	  }
	stand_dev_squares[i] = sqrt (stand_dev_squares[i]);
      }
    
     double*** local_data = new double** [imax];
     for (int i = 0; i < imax; i++) {
       local_data [i] = new double* [4];
     }
     for (int i = 0; i < imax; i++)
       {
	 for (int j = 0; j < 4; j++)
	   {
	     int length = L + 1;
	     local_data [i][j] = new double [length];
	   }
       }
    for (int i = 0; i<imax;i++)
    {
      for (int k = 0; k < L; k++)
	{
        local_data[i][0][L] = sumrho[i];
        local_data[i][0][k] = sumc[i][k];
	local_data[i][1][k] = sumD[i][k];
	local_data[i][1][L] = sums[i];
        local_data[i][2][k] = stand_dev_squarec[i][k];
	local_data[i][2][L] = 0.0;
	local_data[i][3][k] = stand_dev_squareD[i][k];
	local_data[i][3][L] = stand_dev_squares[i];
	}
    }


     for (int i = 0; i < imax; i++)
      {
	for (int j = 0; j < 2; j++)
	  {
	    delete[] ave_cD[i][j]; 
	  }
      }
    for (int i = 0; i < imax; i++)
      {
	delete[] ave_cD[i];
      }
    delete[] ave_cD;
    for (int i = 0; i < imax; i++)
      {
	delete[] ave_rs[i];
      }
      delete[] ave_rs;
    
    delete[] vin;
    delete[] sumrho;
    for (int i = 0; i < imax; i++)
      {
	delete[] sumc[i];
      }
    delete[] sumc;
    for (int i = 0; i < imax; i++)
      {
	delete[] stand_dev_squarec[i];
      }
    delete[] stand_dev_squarec;
    for (int i = 0; i < imax; i++)
      {
	delete[] sumD[i];
      }
    delete[] sumD;
    for (int i = 0; i < imax; i++)
      {
	delete[] stand_dev_squareD[i];
      }
    delete[] stand_dev_squareD;
    
    delete[] sums;
    delete[] stand_dev_squares;
    
    return local_data;
}



/*
Simulation_Loop will take in a mass distribution and first subtract rhostep * L mass.  It will then define the vecor u2 and run Sim_Dyn for a set number of iterations kmax.  Starting at half way through the loop, the loop will begin to sum up c, D, and s after each iteration and averaging it to give the average values after the full loop.  The output will then be an object consisting of the mass distribution after the final iteration (to be used as an input for the next value of rho) and c, D, and s. I have placed a cout before and after MinusMass is called upon in order to find where the loop gets caught.
 */


loop_output Simulation_Loop (double rho, int* vin, int kmax, double omega, int L, double dt, double rhostep, loop_output output)
{
  int* locall = new int [L];
  int rs = rhostep;
  int blah = rs*L;
  locall =  MinusMass (vin, blah, L);
  double* local_c = new double [L];
  double* local_D = new double [L];
  double local_s;
  double s = 0.0;
  double* clocal = new double [L];
  double* Dlocal = new double [L];
  double* D = new double [L];
  for (int i = 0; i < L; i++)
    {
      D[i] = 0.0;
    }
  double* c = new double [L];
  for (int i = 0; i < L; i++)
    {
      c[i] = 0.0;
    }
  int* u2 = new int [L];
  for (int i = 0; i < kmax; i++) {
    locall = Sim_Dyn (u2,locall, omega, L, dt);
    if (i > kmax / 2)
      {
	for (int j = 0; j < L; j++)
	  {
	    local_c[j] = 2.0* ( cx (locall, L, clocal)[j] )/(double (kmax));
	    local_D[j] = 2.0* ( Dx (locall, L, Dlocal)[j] )/(double (kmax));
	    c = dadd(c, local_c, L);
	    D = dadd(D, local_D, L);
	  }
	local_s = 2.0* ( sx (locall, L) )/(double (kmax));
	s += local_s;
      }
  }
  output.c = c;
  output.s = s;
  output.D = D;
  output.mass_distribution = locall;
  delete[] local_c;
  delete[] local_D;
  delete[] Dlocal;
  delete[] clocal;
  delete[] u2;
  return output;
}


/*
cx takes in a mass distribution, and finds c(x) between the origin and every other point and adds it to a running sum.  It is then averaged over the size of the lattice.
 */

double* cx (int* vector, int L, double* clocal)
{
    for (int i = 0; i < L; i++)
      {
	clocal[i] = 0.0;
      }
    for (int i = 0; i<L; i++)
    {
      for (int j = 0; j<L; j++)
	{	 
	   int t2 = i - j;
	  t2 = abs(t2);
	  
	  int tt2 = j - L -i;
	  tt2 = abs(tt2);
	  
	  if (tt2 < t2)
	    {
	      t2 = tt2;
	      }
	  int ttt2 = j + L -i;
	  ttt2 = abs(ttt2);
	  if (ttt2 < t2)
	    {
	      t2 = ttt2;
	    }
	  

	  clocal[t2] += vector[j] * vector[i];
	}
      
    }
    for (int i = 0; i < L; i++)
      {
	clocal[i] = 0.5 * clocal[i]/double(L);
      }	
    return clocal;
}

/*
Dx finds the expectation value D(x) for each point on the lattice and adds it to a running sum.  It is then averaged over the size of the lattice.
*/

double* Dx (int* vector, int L, double* Dlocal)
{
   for (int i = 0; i < L; i++)
     {
       Dlocal[i] = 0.0;
     }
    for (int i = 0; i<L; i++)
      {
	for (int j = 0; j<L; j++)
	  {
	    if (vector[j] == 0)
	      {

		int t = i - j;
		t = abs(t);

		int tt = j - L - i;
		tt = abs(tt);
		if (tt < t)
		  {
		    t = tt;
		  }
		int ttt = j + L -i;
		ttt = abs(ttt);
		if (ttt < t)
		  {
		    t = ttt;
		  }
	  


		Dlocal[t] += vector[i];
	      }
	  }
      }
    for (int i = 0; i < L; i++)
      {
      	Dlocal[i] = 0.5 * Dlocal[i]/double(L);
      }
    return Dlocal;
}

/*
sx finds the probability of a site having mass s.  This is either a one or a zero depending on whether or not the origin has mass.  Averaging this over all of the sites and time gives the probability.
 */

double sx (int* vector, int L)
{
  double slocal;
  for (int j = 0; j < L; j++)
    {
      if (vector[j] == 0)
	{
	  slocal += 0.0;
	}
      else
	{
	  slocal += 1.0;
	}
    }
  slocal = slocal / double(L);
  return slocal;
}
      
      


/*
The main simulation works by taking a steady-state distribution, removing a set amount of mass and using it as the input vector for a new value of rho so that less iterations are required.  MinusMass is a function that will randomly remove a set amount of mass off of the lattice structure.  As inputs we have a mass distribution, the lattice length, and the amount of mass we want removed.  As long as this amount of mass is greater than zero, a random integer will be generated to find a random lattice site, and as long as the site has non-zero mass then a single mass unit will be removed.  After which, one will be subtracted from the original determined amount of mass to be subtracted.  This reoccurs until the amount hits zero.
 */


int* MinusMass (int* vector, int sub, int L)
{
  //cout << "mtest" << "\n";
  double templ = double(L);
  //  cout << "mtest2" << "\n";
  while (sub > 0) {
    int b = floor(randomdouble (0.0, templ));
    // cout << "mtest" << "\t" << sub << "\t" << b << "\n";
    if (vector[b] != 0) {
      vector[b] -= 1;
      sub -= 1;
    }
  }
  return vector;
}

/*
Initial_Run will generate an steady-state distribution of mass.  This is generally run for around 1 to 1.5 million iterations.  It simply runs Sim_Dyn through a for loop and keeps the distribution of the final iteration.
 */


int* Initial_Run (int* vin, int initial_max, double omega, int L, double dt)
{
  int* u2 = new int [L];
  for (int k=1; k<initial_max; k++) {
    vin = Sim_Dyn (u2, vin, omega, L, dt);
  }
  delete[] u2;
  return vin;
}


/*
This is the main set of dynamics.  The given inputs are the mass distribution, an empty vector used to store all of the mass movement which will then be combined with the input distribution, the frequency of chipping to hopping, the lenght, and the time step.  We first define a few different constants based off of omega and dt used as the probabilities.  Then the empty vector u2 is defined to have all zero elements.  The main for loop is then stated.  On the condition that a lattice site has mass, a random number will be generated between zero and one.  A couple of modulus variables are then defined so that the periodic boundary conditions can be used.  For a simple lattice, there is chipping and hopping onto either neighboring sites, so we define variables that will identify a lattice site one site away in either direction.  The random number will then be run through a series of if statements, each considering a certain range of the random number, which acts as a probability.  The first two statments consider the range of probabilities corresponding to hopping in either direction.  In this case, the site's mass will be transported either one spot left or right.  The second two if statements relate to the chipping of a single unit of mass.  Once every site has been acted on, the vector u2 will be added onto the original vector vin and returned.
 */


int* Sim_Dyn (int* u2, int* vinn, double omega, int L, double dt)
{
    double dthalf = 0.5 * dt;
    double dtwhalf = dt * (1.0 + omega*0.5);
    double dtw = dt * (1.0 + omega);
    for (int k = 0; k<L; k++)
      {
	u2[k] = 0;
      }
    for (int j = 0; j < L; j++)
    {
      if (vinn[j] > 0)
        {
	  double  x = randomdouble (0.0, 1.0);
	  int y = mod (j + 1, L);
	  int z = mod (j - 1, L);
	  if ((0.0 < x) && (x < dthalf)) {
	    u2[y] += vinn[j];
	    u2[j] -=  vinn[j];
	  };
	  if ((dthalf < x) && (x < dt)) {
	    u2[z] += vinn[j];
	    u2[j] -= vinn[j];
	  };
	  if ((dt < x) && (x < dtwhalf)) {
	    u2[y] += 1.0;
	    u2[j] -= 1.0;
	  };
	  if ((dtwhalf < x) && (x < dtw)) {
	    u2[z] += 1.0;
	    u2[j] -= 1.0;
	  };
        }
    }
    vinn =  add (vinn, u2, L);
    return vinn;
}


/*
This generates the intial vector mass distribution used in this simulation.  It takes all of the mass and places it on the middle site.  It can be easily modified to randomly distribute mass throughout the lattice structure.
 */

int* Initial (double rho, int L)
{
    int* out_local = new int [L];
    for (int i=0; i< L; i++)
    {
 // Below the middle element has all of the mass placed on it, the rest of the sites are defined to be zero.
        if (i== floor (L/2)) {
            out_local[i] = rho * L;
        }
        else {
            out_local[i] = 0;
        }
    }
    return out_local;
}


/*
This is just a function to add two vectors, one a for integer vectors and one for double vectors.  Given an input of the two vectors and the vector length, each element of the second vector is added onto the relevent element of the first vector.  After the first vector is replaced, it is given as the output.
 */


double*  dadd (double *u, double* v, int L)
{
    for (int i=0; i<L; i++)
    {
       u[i] += v[i];
    }
    return u;
}

int*  add (int *u, int* v, int L)
{
    for (int i=0; i<L; i++)
    {
       u[i] += v[i];
    }
    return u;
}

/*
This is a modulus that will enforce the periodic boundary conditions of the lattice structure. While considering all of the elements of the vector representing the lattice, if there is an attempt to find the element right before the first or right after the last element, the element will be reassigned to, respectively, the last or first element.
 */


int  mod (int a, int L)
{
    if (a==-1) {
        a = L-1;
    }
    else if (a == L) {
        a = 0;
    }
    return a;
}

/* 
randomdouble will generate a random number in between a designated minimum and maximum value.
This is achieved using the standard library's rand() and RAND_MAX, where a random number between zero and RAND_MAX will be generated.  This random number, which is a double, is then scaled by multiplying it by the difference of the input max and min values.  By adding the minimum to this scaled random number, we now have a random number in between the minimum and maximum values.
 */


double randomdouble (double min, double max)
{
    double random = ((double) rand()) / (double) RAND_MAX;
    // returns a random double number in [0 , 1]
    double diff = max - min;
    double r = random * diff;
    //expands [0,1] to [0,max-min]
    return min+r;
    // adjusts to [min, max] and returns
}
