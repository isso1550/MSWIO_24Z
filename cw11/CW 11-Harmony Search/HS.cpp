// HS algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-18


#include <math.h>

#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "HS.H"
#include "problem.h"



extern double Best_solution[1000]; /* Data of best solution */

 HS_data dHS;

 harmony::harmony(void)
 {
     
 }

 harmony::~harmony(void)
 {
   delete [] x;
 }


void harmony::init(int dim_of_problem, problem_limits *v_l)
 {

   dim=dim_of_problem;
   x=new double[dim];
   
   if(x==NULL)  { printf("Memory allocation error 001"); exit(0);}

   int i;
    for(i=0;i<dim;i++)
     {
      x[i]=0;
     }
    fitness=0;
    objective=0;

   // notes initialization   
   for(i=0;i<dim;i++) 
    {
     if(v_l[i].tn==1) 
     x[i]=RandVal(v_l[i].min,v_l[i].max);
     else
     x[i]=0;
    }

 }

    

 HS::~HS()
 {
     delete [] p;
     delete [] p_copy;
     delete [] pro_lim;
 }


 HS::HS(void)
 {
     int i;
     num_dim=dHS.num_dim;
     pop_size=dHS.pop_size;
     num_iter=dHS.num_iter;

	 acc_rate=dHS.acc_rate; // accept rate
	 pa_rate=dHS.pa_rate; // pitch (note) adjustement rate 
     pab=dHS.pab; // pitch adjustment coeficient - the value to divide the range of n-th dimension and have the pitch adjustment range

     p=new harmony[pop_size];                if(p==NULL)  { printf("Memory allocation error 004"); exit(0);}
     p_copy=new harmony[pop_size];           if(p_copy==NULL)  { printf("Memory allocation error 005"); exit(0);}
     pro_lim=new problem_limits[num_dim];    if(pro_lim==NULL)  { printf("Memory allocation error 006"); exit(0);}
     set_problem_limits();
 }     


void HS::init_harmonies(void)
 {
   int i;
   
     for(i=0;i<pop_size;i++)
      {
          p[i].init(num_dim, pro_lim);
          p_copy[i].init(num_dim, pro_lim);
      }
     p_best.init(num_dim, pro_lim);

  p_best.fitness=-1e300; 
 }

 void HS::evaluate_harmonies(void)
 {
  int i;   
  for(i=0;i<pop_size;i++)
   {
    p[i].calc_fitness();
   }
 }

 void HS::remember_best_of_harmonies(void)
 {
  int i,j;   
  for(i=0;i<pop_size;i++)
   {
    if(p[i].fitness>p_best.fitness) 
     {
       p_best.fitness=p[i].fitness;
       p_best.objective=p[i].objective;
       dHS.t_best = clock() - dHS.ts;
       dHS.fe_best = dHS.fe_curr;
	   	      
       for(j=0;j<num_dim;j++)
        {
          p_best.x[j]=p[i].x[j];
        }
      }
   }
 }

 void HS::generate_new_harmonies(void)
 {
  int i,j, k;
 
  double r;

 
  for(i=0;i<pop_size;i++)
   {
       
   for(j=0;j<num_dim;j++)
    {
     r=RandVal(0.0,1.0); 
     
     if(r<acc_rate)  // accept current note in the harmony from among rest of harmonies
     {
     	 k=iRandVal(0, pop_size-1);
     	 p_copy[i].x[j]=p[k].x[j];
	 }
	 else
	 {
	    r=RandVal(0.0,1.0); 
	 	if(r<pa_rate) // pitch adjust - a small move aroound current position
	 	{
	 		p_copy[i].x[j]=p[i].x[j]  +  ((pro_lim[j].max-pro_lim[j].min)/pab) * (RandVal(-0.5, 0.5));
		}
		else
		{
	 		p_copy[i].x[j]=	RandVal(pro_lim[j].min,pro_lim[j].max);
		}
	 }

     if(p_copy[i].x[j]<pro_lim[j].min) p_copy[i].x[j]=pro_lim[j].min;//keeping the notes in the limits of the problem being optimized
     if(p_copy[i].x[j]>pro_lim[j].max) p_copy[i].x[j]=pro_lim[j].max;
     if(pro_lim[j].tn==0) p_copy[i].x[j]=0;
    }

	p_copy[i].calc_fitness();
   }


  for(i=0;i<pop_size;i++)
   {

    if(p_copy[i].fitness > p[i].fitness)
     {
       p[i].fitness  =p_copy[i].fitness;   
       p[i].objective=p_copy[i].objective;   
       
       for(j=0;j<num_dim;j++)
        {
          p[i].x[j]=p_copy[i].x[j];
        }     	
	 }

   }
 }


 void HS::make_harmony_copy(void)
 { 
  int i,j;   
  for(i=0;i<pop_size;i++)
   {
       p_copy[i].fitness=p[i].fitness;   
       p_copy[i].objective=p[i].objective;   
       
       for(j=0;j<num_dim;j++)
        {
          p_copy[i].x[j]=p[i].x[j];
        }
   }
 }






// Main function of swarm simulation

int HS_simulation(void)
{
   
int i,j,k,z;
FILE *galog;

long status;
char buforek[100]; 
  double avg;
  double stddev;
  double sum_square;
  double square_sum;
  double sum;


  if((galog=fopen("log.txt","wt"))==NULL)
    {
     printf("\n Can't open output file log.txt");
     exit(1);
    }


  srand(dHS.p_rand);
  
  dHS.ts = clock();
  dHS.fe_curr=0;

  HS hs;


  fprintf(galog,"%i\n",hs.pop_size);
  fprintf(galog,"%i\n",hs.num_iter);
  fprintf(galog,"%f|%f|%f\n",hs.acc_rate,hs.pa_rate,hs.pab);
  fprintf(galog,"%i|%i| %i\n",hs.iter_danger,hs.neighborhood, hs.num_dim ); 

  z=0;


  hs.init_harmonies();
 
  hs.evaluate_harmonies();
  hs.remember_best_of_harmonies();
 
  
  printf("\n Starting simulation \n\n");


  
 for(k=0;k<hs.num_iter;k++)
  {
    hs.generate_new_harmonies();     
    hs.remember_best_of_harmonies();


    sum=0.0;
    sum_square=0.0;

    for(i=0;i<hs.pop_size;i++)
     {
      sum=sum+hs.p[i].objective;
      sum_square+=hs.p[i].objective*hs.p[i].objective;
     }

    avg=sum/(double)hs.pop_size;
    square_sum=sum*sum/(double)hs.pop_size;
    stddev=sqrt((1.0/(double)(hs.pop_size-1))*(sqrt((sum_square-square_sum)*(sum_square-square_sum))));

    dHS.t_all = clock() - dHS.ts;

    fprintf(galog,"%d,%9e,%9e,%9e,%f,%li\n",k,hs.p_best.objective,sum,stddev,(((double)dHS.t_all)/CLOCKS_PER_SEC),dHS.fe_curr);

    status++;
   
    if(status==100)
     {
          status=0;
          // update status on dialog                  
     }
     
     printf("\r %i z %i (%i)",k+1,hs.num_iter, dHS.sim_num);

 } //end of for(k=0;k<hs.num_iter;k++)

   fprintf(galog,"\n\n HS simulation finished\n");
  
   fprintf(galog,"\nBest Fitness=%e,",hs.p_best.objective);
   for(j=0;j<hs.num_dim;j++)
    {
      fprintf(galog,"x[%i]=%f|,",j,hs.p_best.x[j]);
    }

   dHS.best_objective=hs.p_best.objective;

   for(i=0;i<hs.num_dim; i++)
   {
    Best_solution[i]=hs.p_best.x[i];
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



