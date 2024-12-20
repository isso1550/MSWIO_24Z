// SA algorithm implementation (c) Dariusz Baczy�ski - Warsaw University Of Technology
// 
// 2020-11-27

#include <math.h>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "sa.H"
#include "problem.h"

extern int    NTRANSF;

extern double Best_solution[1000]; /* Data of best solution */

 SA_data dSA;

 individual::individual(void)
 {
     
 }

 individual::~individual(void)
 {
   delete [] x;
 }


void individual::init(int dim_of_problem, problem_limits *v_l)
 {

   dim=dim_of_problem;
   x=new double[dim];
   
   if(x==NULL)  { printf("Memory allocation error 001"); exit(0);}

   int i;
    for(i=0;i<dim;i++)
     {
      x[i]=0;
     }
     
    fitness=-1e300;
    objective=-1e300;

   // genome initialization   
   for(i=0;i<dim;i++) 
    {
     if(v_l[i].tn==1) 
     x[i]=RandVal(v_l[i].min,v_l[i].max);
     else
     x[i]=0;
    }

 }


     

 SA::~SA()
 {
     delete [] pop;
     delete [] pop2;
 }


 SA::SA(void)
 {
     int i;
     
    num_iter=dSA.num_iter;     // number of iterations MAXGENS

	start_temp=dSA.start_temp;  // temperature at the start of annealing process
	gamma=dSA.gamma;       // gamma coeficient
	boltz=dSA.boltz;       // Boltzman's coeficient
    pop_size=dSA.pop_size;

    num_dim=dSA.num_dim;      // number of problem dimensions
     
     pop=new individual[pop_size];           if(pop==NULL)  { printf("Memory allocation error 004"); exit(0);}
     pop2=new individual[pop_size];          if(pop2==NULL)  { printf("Memory allocation error 005"); exit(0);}
     pro_lim=new problem_limits[num_dim];    if(pro_lim==NULL)  { printf("Memory allocation error 006"); exit(0);}
     set_problem_limits();
 }     


void SA::init_population(void)
 {
   int i;
   
     for(i=0;i<pop_size;i++)
      {
          pop[i].init(num_dim, pro_lim);
          pop2[i].init(num_dim, pro_lim);
      }
     p_best.init(num_dim, pro_lim);

  p_best.fitness=-1e300; 
  p_best.objective=-1e300; 
 }



 void SA::evaluate_population(void)
 {
  int i; 
     
  for(i=0;i<pop_size;i++)
   {
    pop[i].calc_fitness();
   // if(pop[i].fitness<0.0) { printf("\n ****** Error - fitness function lower than 0.0");    }
   }
 }



 void SA::remember_best_of_population(void)
 {
  int i,j;  


   if(pop[0].fitness>p_best.fitness)
    {
       p_best.fitness=  pop[0].fitness;
       p_best.objective=pop[0].objective;
       dSA.t_best = clock() - dSA.ts;
       dSA.fe_best = dSA.fe_curr;
	      
       for(j=0;j<num_dim;j++)
        {
          p_best.x[j]=pop[0].x[j];
        }
    }
 }





/*********************** create new solution ********************************/
void SA::create_new_solution(void)
{
	
  int i,j;

     
  for(i=0;i<pop_size;i++)
   {
    for(j=0;j<num_dim;j++) 
     {
     	//pop2[i].x[j]= pop[i].x[j] + RandVal (-1.0, 1.0);
    	pop2[i].x[j]= pop[i].x[j] + RandVal (-50.0, 50.0);   // create a new solution just by adding some small number to current position
    	// but you can develope a another way to create new solutions RandVal(pro_lim[var].min,pro_lim[var].max)
    	//this way does not take into account the boundaries of the problem so we can go out of the box
     }
   	
    pop2[i].calc_fitness(); // evaluate new solution
   }


/*
 if(pro_lim[var].tn==1) 
     pop[member].x[var]=RandVal(pro_lim[var].min,pro_lim[var].max);
   else
     pop[member].x[var]=0;
*/

}



/*********************** decide if new solution should be accepted as the current solution ********************************/
int SA::decide (void)
{
	
  int i,j;
  double p, del_E, r;

	if(pop2[0].fitness>pop[0].fitness)
	{
     for(j=0;j<num_dim;j++) 
    	{
    		pop[0].x[j] = pop2[0].x[j];   // new solution is better than the old one so we set it as current
     	}
	 return(1);
	}
	else
	{
	 del_E=gamma*(pop2[0].fitness-pop[0].fitness);	// eg. -7 and -5 = -2   or    10 and 100 = -90
	 p=exp(del_E/(boltz*Temp));
	 r=RandVal(0,1);
	 
	 if(p>r)
	 	{
	 		return (0); // do nothing 
		}
		else
		{
		    for(j=0;j<num_dim;j++) 
    		{
    		  pop[0].x[j] = pop2[0].x[j];   // new solution is not better but it is set as current
     		}	
     		return(1);
		}
	}


}













// Main function of annealing simulation

int SA_simulation(void)
{
   
int i,j,k,z;
FILE *galog;

char buforek[100]; 
  double avg=0;
  double stddev=0;
  double sum_square=0;
  double square_sum=0;
  double sum=0;

  double T_step;
  
 

  if((galog=fopen("log.txt","wt"))==NULL)
    {
     printf("\n Can't open output file log.txt");
     exit(1);
    }


  srand(dSA.p_rand);

  dSA.ts = clock();
  dSA.fe_curr=0;


  SA sa;


  fprintf(galog,"%i\n",sa.pop_size);
  fprintf(galog,"%i\n",sa.num_iter);
  fprintf(galog,"%f|%f|%f|%li|\n", sa.start_temp, sa.gamma, sa.boltz, sa.num_iter);
  
  fprintf(galog,"%i|\n", sa.num_dim ); 

  z=0;

  sa.Temp=sa.start_temp;
  T_step=sa.start_temp/sa.num_iter; // we assume that the system should cool down to 0 degrees during num_iter


  sa.init_population();
  sa.evaluate_population();  // evaluate population of solution - in case of SA - it is just one solution
  		sa.remember_best_of_population(); // remember best solution for report file




  printf("\n Starting simulation \n\n");


  
 for(k=0;k<sa.num_iter;k++)
  {
   
   
    sa.create_new_solution();
	sa.decide();
    
	sa.Temp=sa.Temp-T_step;
    
    sa.evaluate_population();
    		sa.remember_best_of_population();  // remember best solution for report file
 


/*
// not useful when one solution is processed at each iteration
    sum=0.0;
    sum_square=0.0;
    avg=0.0; 

    for(i=0;i<sa.pop_size;i++)
     {
      sum=sum+sa.pop[i].objective;
      sum_square+=sa.pop[i].objective*sa.pop[i].objective;
     }

    avg=sum/(double)sa.pop_size;
    square_sum=sum*sum/(double)sa.pop_size;
    stddev=sqrt((1.0/(double)(sa.pop_size-1))*(sqrt((sum_square-square_sum)*(sum_square-square_sum))));
*/

    dSA.t_all = clock() - dSA.ts;


    fprintf(galog,"%d,%9e,%9e,%9e,%9e,%f,%li\n",k,sa.p_best.objective,sa.pop[0].objective,sum,stddev,(((double)dSA.t_all)/CLOCKS_PER_SEC),dSA.fe_curr);

    
    printf("\r %i z %i (%i)",k+1,sa.num_iter, dSA.sim_num);

	if (sa.Temp<0) break;


 } //end of for(k=0;k<sa.num_iter;k++)






   fprintf(galog,"\n\n SA simulation finished\n");
  
   fprintf(galog,"\nBest Fitness=%e,",sa.p_best.objective);
   for(j=0;j<sa.num_dim;j++)
    {
      fprintf(galog,"x[%i]=%e|,",j,sa.p_best.x[j]);
    }

   dSA.best_objective=sa.p_best.objective;

   for(i=0;i<sa.num_dim; i++)
   {
    Best_solution[i]=sa.p_best.x[i];
   }

  fclose(galog);

  return 0;
}



// Additional functions

/************************ RandVal *****************************************/
double RandVal(double low,double high)
{
  return(((double) (rand() % RAND_MAX)/((double) RAND_MAX))*(high-low)+low);
}


/************************ iRandVal *****************************************/
long iRandVal(long low,long high) 
{
  high++;
  return (long) (((double) (rand() % RAND_MAX)/((double) RAND_MAX))*(high-low)+low);
}



