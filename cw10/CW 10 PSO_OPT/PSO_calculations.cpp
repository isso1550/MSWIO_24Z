// PSO algorithm implementation (c) Dariusz Baczy�ski - Warsaw University Of Technology
// 
// 2020-12-12

#include <cstring>
#include <math.h>
#include <ctime>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "PSO.H"
#include "problem.h"




extern int    NVARS;

extern PSO_data dPSO;

extern problem_limits z_z [1000];
extern problem_data p_d[200]; // not used

extern double Best_solution[1000]; /* Data of best solution */









// MENU AND CONTROL COMMAND
int PSO_optimization(void)
{
FILE *auto_PSO, *log;
char buforek[2000];
int i;
    
data_read();
 
printf("\n\n\n Data loading \n");    
    
  
    dPSO.c0=1.7; // coef for current speed 
    dPSO.c1=0.1; // coef for best position of this particle (found so far)
    dPSO.c2=0.1; // coef for best position among all particles (found so far)
    dPSO.c3=0.1; // coef for best position among neighbor particles (found so far)
    dPSO.iter_danger=160; // number of iterations between dangerous situation for particles(animals) - scattering/dispersion of swarm
    dPSO.neighborhood=3; // neighborhood width

    dPSO.pop_size=50;
    dPSO.num_iter=100000;
    dPSO.num_dim=NVARS;

    dPSO.p_rand=0;
    
    dPSO.best_objective=0;
    dPSO.sim_num=0;


   if((auto_PSO=fopen("auto_PSO.bxt","rt"))!=NULL)
    {
     if((log=fopen("report_PSO_auto.txt","wt"))==NULL)
       {
        printf("\n Can't open output file report_PSO_auto.txt");
        exit(1);
       }

        fgets(buforek,1998,auto_PSO); // omit 2 starting lines
        fgets(buforek,1998,auto_PSO);

        while(fgets(buforek,1998,auto_PSO))
         {
           if(strlen(buforek)<5) continue;
           sscanf(buforek,"%f|%f|%f|%f|%li|%i|%li|%li|%li|", &dPSO.c0, &dPSO.c1, &dPSO.c2, &dPSO.c3, &dPSO.iter_danger, &dPSO.neighborhood, &dPSO.pop_size, &dPSO.num_iter, &dPSO.p_rand);

           dPSO.sim_num++;
            
           PSO_simulation();
           fprintf(log,"%f|%f|", (((double)dPSO.t_best)/CLOCKS_PER_SEC), (((double)dPSO.t_all)/CLOCKS_PER_SEC) );
           fprintf(log,"%li|%li|", dPSO.fe_best, dPSO.fe_curr );           

           fprintf(log,"%f|%f|%f|%f|%i|%i|%i|%i|%i|", dPSO.c0, dPSO.c1, dPSO.c2, dPSO.c3, dPSO.iter_danger, dPSO.neighborhood, dPSO.pop_size, dPSO.num_iter, dPSO.p_rand);

           fprintf(log, "%e|" , dPSO.best_objective);

           for(i=0;i<NVARS; i++)
            {
             fprintf(log, "%e|" , Best_solution[i]);
            }

           fprintf(log, "\n" );
           fflush(log);

         }

     
     fclose(log);
     fclose(auto_PSO);
    }
   else
    {
    PSO_simulation();

    }


    return 0;
}



 

  void PSO::set_problem_limits(void)
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
void particle::calc_fitness(void) //Schwefel's function
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
  dPSO.fe_curr++;
 }
*/


 void particle::calc_fitness(void) //Rastrigin's function
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
  dPSO.fe_curr++;
 }



/*
 void particle::calc_fitness(void) //Griewank's function
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
  dPSO.fe_curr++;
 }
*/

/*
 void particle::calc_fitness(void) //De Jong's function
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
  dPSO.fe_curr++;
 }
*/
 
/*
 void particle::calc_fitness(void) //Ackley�s function
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
  dPSO.fe_curr++;
 }
 */
 
 
