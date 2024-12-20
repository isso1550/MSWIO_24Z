// SA algorithm implementation (c) Dariusz Baczy�ski - Warsaw University Of Technology
// 
// 2020-11-27


#include <cstring>
#include <math.h>
#include <ctime>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "sa.H"
#include "problem.h"



extern int    NTRANSF;
extern int    NVARS;

extern SA_data dSA;

extern problem_limits z_z [1000];

extern double Best_solution[1000]; /* Data of best solution */









// MENU AND CONTROL COMMAND
int SA_optimization(void)
{
FILE *auto_SA, *log;
char buforek[2000];
int i;
    
data_read();
 
printf("\n\n\n Data loading \n");    
    
    dSA.pop_size=1;     // simulated annealing works using one solution of the problem 
    dSA.num_iter=10000;   
    
    dSA.num_dim=NVARS;   
    dSA.p_rand=0;       
    
    dSA.best_objective=0;
    dSA.sim_num=0;


   if((auto_SA=fopen("auto_SA.bxt","rt"))!=NULL)
    {
     if((log=fopen("report_SA_auto.txt","wt"))==NULL)
       {
        printf("\n Can't open output file report_SA_auto.txt");
        exit(1);
       }


        fgets(buforek,1998,auto_SA); // omit 2 starting lines
        fgets(buforek,1998,auto_SA);

        while(fgets(buforek,1998,auto_SA))
         {
           if(strlen(buforek)<5) continue;
           
/*  
    float start_temp;  // temperature at the start of annealing process
    float gamma;       // gamma coeficient
    float boltz;       // Boltzman's coeficient
*/
           sscanf(buforek,"%f|%f|%f|%li|%i|", &dSA.start_temp, &dSA.gamma, &dSA.boltz, &dSA.num_iter, &dSA.p_rand);
           
           dSA.sim_num++;
            
           SA_simulation();

           fprintf(log,"%f|%f|", (((double)dSA.t_best)/CLOCKS_PER_SEC), (((double)dSA.t_all)/CLOCKS_PER_SEC) );
           fprintf(log,"%li|%li|", dSA.fe_best, dSA.fe_curr );

           fprintf(log,"%f|%f|%f|%li|%i|", dSA.start_temp, dSA.gamma, dSA.boltz, dSA.num_iter, dSA.p_rand);

           fprintf(log, "%e|" , dSA.best_objective);

           for(i=0;i<NVARS; i++)
            {
             fprintf(log, "%e|" , Best_solution[i]);
            }

           fprintf(log, "\n" );
           fflush(log);
         }

     
     fclose(log);
     fclose(auto_SA);
    }
   else
    {
    SA_simulation();

    }


    return 0;
}



 

  void SA::set_problem_limits(void)
   {
   	int i;
   	
     for(i=0; i<num_dim; i++) 
	   { 
	      pro_lim[i].tn=z_z[i].tn; 
	      pro_lim[i].min=z_z[i].min;
	      pro_lim[i].max=z_z[i].max;
       }   	
   }


/*
void individual::calc_fitness(void) //Paraboloid function
 {
  double sum1=0.0;
  int i;
  
   
   for(i=0;i<dim;i++) // paraboloid Y=-X^2 (n - dimensional) - in this aplication one dimensional
    {
    	sum1 = sum1 - x[i]*x[i];
    }
  
  
  objective=sum1;
  fitness=sum1;
  dSA.fe_curr++;
 }
*/



/*

void individual::calc_fitness(void) //Schwefel's function
 {
  double sum1=0.0, sum2=0.0;
  double sum=0.0;
  int i;
  const double PI = 3.141592653589793;
   
   
   for(i=0;i<dim;i++) //Schwefel's
    {
    	sum1 = sum1 + x[i]*sin(sqrt(fabs(x[i])));
    }
  
  sum=-sum1/dim + 500.0;
  
  objective=sum;
  fitness=1.0/(sum+1.0);
  dSA.fe_curr++;
 }
*/
/*
void individual::calc_fitness(void) //Rastrigin's function
 {
  double sum1=0.0, sum2=0.0;
  double sum=0.0;
  int i;
  const double PI = 3.141592653589793;
   
   
   for(i=0;i<dim;i++) //Rastrigin's
    {
    	sum1 = sum1 + pow(x[i],2.0) - 10.0*cos(2.0*PI*x[i]);
    }
  
  sum=10.0*dim+sum1;
  
  objective=sum;
  fitness=1.0/(sum+1.0);
  dSA.fe_curr++;
 }
*/

/*
void individual::calc_fitness(void) //Griewank's function
 {
  double sum1=0.0, sum2=0.0;
  double sum=0.0;
  int i;
  const double PI = 3.141592653589793;
   
   sum2=1.0;
   
   for(i=0;i<dim;i++) //Griewank's
    {
    	sum1 = sum1 + pow(x[i],2.0);
    	sum2 = sum2*cos((x[i])/sqrt((i+1)));	
    }
  
  sum1 = sum1/4000.0;
  sum=sum1-sum2+1.0;
  
  objective=sum;
  fitness=1.0/(sum+1.0);
  dSA.fe_curr++;
 }
*/

/*
void individual::calc_fitness(void) //De Jong's function
 {
  double sum1=0.0, sum2=0.0;
  double sum=0.0;
  int i;
  const double PI = 3.141592653589793;
   
   for(i=0;i<dim;i++) //De Jong's
    {
    	sum = sum + pow(x[i],2.0);
    }
  
  objective=sum;
  fitness=1.0/(sum+1.0);
  dSA.fe_curr++;
 }
*/

/*
void individual::calc_fitness(void) //Rosenbrock's function
 {
  double sum1=0.0, sum2=0.0;
  double sum=0.0;
  int i;
  const double PI = 3.141592653589793;
   
   for(i=0;i<(dim-1);i++) // Rosenbrock's
    {
    	sum = sum + pow((x[i]-1.0),2.0) + 100.0*pow((x[i+1]-x[i]*x[i]),2.0);
    }
  
  objective=sum;
  fitness=1.0/(sum+1.0);
  dSA.fe_curr++;
 }
*/



 void individual::calc_fitness(void) //Ackley�s function
 {
  double sum1=0.0, sum2=0.0;
  double sum;
  int i;
  const double PI = 3.141592653589793;
   
   for(i=0;i<dim;i++)
    {
    	sum1=sum1+x[i]*x[i];
    	sum2=sum2+cos(x[i]*2*PI);
    }

   sum=-20.0*exp(-1.0/5.0*sqrt(sum1/((double)dim))) - exp(sum2/((double)dim)) + 20.0 + exp(1.0);                       // Standard Ackley's function
//   sum=-20.0*exp(-1.0/5.0*sqrt(sum1/((double)dim))) - exp(sum2/((double)dim)) + 20.0 + exp(1.0) + fabs(tan(sum1));     // Ackley's function with discontinuities
  
  objective=sum;
  fitness=1.0/(sum+1.0);
  dSA.fe_curr++;
 }
 
 
 
 
 
