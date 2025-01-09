// HS algorithm implementation (c) Dariusz Baczy�ski - Warsaw University Of Technology
// 
// 2020-12-18

#include <cstring>
#include <math.h>
#include <ctime>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "HS.H"
#include "problem.h"




extern int    NVARS;

extern HS_data dHS;

extern problem_limits z_z [1000];
extern problem_data p_d[200]; // not used

extern double Best_solution[1000]; /* Data of best solution */









// MENU AND CONTROL COMMAND
int HS_optimization(void)
{
FILE *auto_HS, *log;
char buforek[2000];
int i;
    
data_read();
 
printf("\n\n\n Data loading \n");    
    
  
    dHS.acc_rate=0.95; // accept rate
    dHS.pa_rate=0.7; // pitch (note) adjustement rate 
    dHS.pab=200; // pitch adjustment coeficient - the value to divide the range of n-th dimension and have the pitch adjustment range

    dHS.pop_size=50; // number of harmonies at each iteration of algorithm
    dHS.num_iter=100000; // number of iterations
    dHS.num_dim=NVARS;

    dHS.p_rand=0;
    
    dHS.best_objective=0;
    dHS.sim_num=0;


   if((auto_HS=fopen("auto_HS.bxt","rt"))!=NULL)
    {
     if((log=fopen("report_HS_auto.txt","wt"))==NULL)
       {
        printf("\n Can't open output file report_HS_auto.txt");
        exit(1);
       }

        fgets(buforek,1998,auto_HS); // omit 2 starting lines
        fgets(buforek,1998,auto_HS);

        while(fgets(buforek,1998,auto_HS))
         {
           if(strlen(buforek)<5) continue;
           sscanf(buforek,"%f|%f|%f|%li|%li|%li|", &dHS.acc_rate, &dHS.pa_rate, &dHS.pab , &dHS.pop_size, &dHS.num_iter, &dHS.p_rand);

           dHS.sim_num++;
            
           HS_simulation();
           fprintf(log,"%f|%f|", (((double)dHS.t_best)/CLOCKS_PER_SEC), (((double)dHS.t_all)/CLOCKS_PER_SEC) );
           fprintf(log,"%li|%li|", dHS.fe_best, dHS.fe_curr );           

           fprintf(log,"%f|%f|%f|%li|%li|%li|", dHS.acc_rate, dHS.pa_rate, dHS.pab , dHS.pop_size, dHS.num_iter, dHS.p_rand);

           fprintf(log, "%e|" , dHS.best_objective);

           for(i=0;i<NVARS; i++)
            {
             fprintf(log, "%e|" , Best_solution[i]);
            }

           fprintf(log, "\n" );
           fflush(log);

         }

     
     fclose(log);
     fclose(auto_HS);
    }
   else
    {
    HS_simulation();
    }


    return 0;
}



 

  void HS::set_problem_limits(void)
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
void harmony::calc_fitness(void) //Paraboloid function
 {
  double sum1=0.0;
  int i;
  
   
   for(i=0;i<dim;i++) // paraboloid Y=-X^2 (n - dimensional) - in this aplication one dimensional
    {
    	sum1 = sum1 - x[i]*x[i];
    }
  
  
  objective=sum1;
  fitness=sum1;
  dHS.fe_curr++;
 }
 */
 
/*
  void harmony::calc_fitness(void) //Ackley�s function
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
  dHS.fe_curr++;
 }
*/

/*
 void harmony::calc_fitness(void) //Rastrigin's function
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
  dHS.fe_curr++;
 }
 */

  void harmony::calc_fitness(void) //De Jong's function
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
  dHS.fe_curr++;
 }
 
