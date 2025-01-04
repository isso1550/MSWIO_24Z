// EA algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-04


#include <cstring>
#include <math.h>
#include <ctime>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "ea.H"
#include "problem.h"



extern int    NTRANSF;
extern int    NVARS;

extern EA_data dEA;

extern problem_limits z_z [1000];
extern struct dane_o_sieci dane[200];

extern double Best_solution[1000]; /* Data of best solution */









// MENU AND CONTROL COMMAND
int EA_optimization(void)
{
FILE *auto_EA, *log;
char buforek[2000];
int i;
    
data_read();
 
printf("\n\n\n Data loading \n");    
    
    dEA.pop_size=50;    
    dEA.num_iter=10000;   
    dEA.px=0.7;          
    dEA.pm=0.1;         

    dEA.scal=0;           
    dEA.elit=0;           
    
    dEA.aryt_xo=0;        
    dEA.num_dim=NVARS;   
    dEA.p_rand=0;       
    
    dEA.best_objective=0;
    dEA.sim_num=0;


   if((auto_EA=fopen("auto_EA.bxt","rt"))!=NULL)
    {
     if((log=fopen("report_EA_auto.txt","wt"))==NULL)
       {
        printf("\n Can't open output file report_EA_auto.txt");
        exit(1);
       }


        fgets(buforek,1998,auto_EA); // omit 2 starting lines
        fgets(buforek,1998,auto_EA);

        while(fgets(buforek,1998,auto_EA))
         {
           if(strlen(buforek)<5) continue;
           sscanf(buforek,"%i|%li|%f|%f|%i|%i|%i|%i|", &dEA.pop_size, &dEA.num_iter, &dEA.px, &dEA.pm, &dEA.aryt_xo, &dEA.scal, &dEA.elit, &dEA.p_rand);
           
           dEA.sim_num++;
            
           EA_simulation();

           fprintf(log,"%f|%f|", (((double)dEA.t_best)/CLOCKS_PER_SEC), (((double)dEA.t_all)/CLOCKS_PER_SEC) );
           fprintf(log,"%li|%li|", dEA.fe_best, dEA.fe_curr );

           fprintf(log,"%i|%li|%f|%f|%i|%i|%i|%i|", dEA.pop_size, dEA.num_iter, dEA.px, dEA.pm, dEA.aryt_xo, dEA.scal, dEA.elit, dEA.p_rand);

           fprintf(log, "%e|" , dEA.best_objective);

           for(i=0;i<NVARS; i++)
            {
             fprintf(log, "%e|" , Best_solution[i]);
            }

           fprintf(log, "\n" );
           fflush(log);
         }

     
     fclose(log);
     fclose(auto_EA);
    }
   else
    {
    EA_simulation();

    }


    return 0;
}



 

  void EA::set_problem_limits(void)
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
  dEA.fe_curr++;
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
  dEA.fe_curr++;
 }
/*

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
  dEA.fe_curr++;
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
  dEA.fe_curr++;
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
  dEA.fe_curr++;
 }
*/



 void individual::calc_fitness(void) //Ackley’s function
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
  dEA.fe_curr++;
 }
 
 
 
 
 
