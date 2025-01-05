// PSO algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-12

#include <math.h>

#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "PSO.H"
#include "problem.h"



extern double Best_solution[1000]; /* Data of best solution */

 PSO_data dPSO;

 particle::particle(void)
 {
     
 }

 particle::~particle(void)
 {
   delete [] x;
   delete [] v;
   delete [] x_best;
 }


void particle::init(int dim_of_problem, problem_limits *v_l)
 {

   dim=dim_of_problem;
   x=new double[dim];
   v=new double[dim];
   x_best=new double[dim];
   
   if(x==NULL)  { printf("Memory allocation error 001"); exit(0);}
   if(v==NULL)  { printf("Memory allocation error 002"); exit(0);}
   if(x_best==NULL)  { printf("Memory allocation error 003"); exit(0);}

   int i;
    for(i=0;i<dim;i++)
     {
      x[i]=0;
      v[i]=0;
      x_best[i]=0;
     }
    fitness=0;
    objective=0;
    fitness_best=-1e300;  

   // position initialization   
   for(i=0;i<dim;i++) 
    {
     if(v_l[i].tn==1) 
     x[i]=RandVal(v_l[i].min,v_l[i].max);
     else
     x[i]=0;
    }

   // speed initialization 
   for(i=0;i<dim;i++) // dim = NVARS
    {
     if(v_l[i].tn==1) 
     v[i]=x[i]-RandVal(v_l[i].min,v_l[i].max);
     else
     v[i]=0;
    }
 }

 
 void particle::remember_best_p(void) 
 {
     int i;
     if(fitness>fitness_best)
      {
       fitness_best=fitness;   
       for(i=0;i<dim;i++)
        {
          x_best[i]=x[i];
        }
      }
 }

    

 PSO::~PSO()
 {
     delete [] p;
     delete [] p_copy;
     delete [] pro_lim;
 }


 PSO::PSO(void)
 {
     int i;
     num_dim=dPSO.num_dim;
     pop_size=dPSO.pop_size;
     num_iter=dPSO.num_iter;

     c0=dPSO.c0; 
     c1=dPSO.c1; 
     c2=dPSO.c2; 
     c3=dPSO.c3; 
     iter_danger=dPSO.iter_danger;
     neighborhood=dPSO.neighborhood;

     p=new particle[pop_size];               if(p==NULL)  { printf("Memory allocation error 004"); exit(0);}
     p_copy=new particle[pop_size];          if(p_copy==NULL)  { printf("Memory allocation error 005"); exit(0);}
     pro_lim=new problem_limits[num_dim];    if(pro_lim==NULL)  { printf("Memory allocation error 006"); exit(0);}
     set_problem_limits();
 }     


void PSO::init_swarm(void)
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

 void PSO::evaluate_swarm(void)
 {
  int i;   
  for(i=0;i<pop_size;i++)
   {
    p[i].calc_fitness();
    p[i].remember_best_p();
   }
 }

 void PSO::remember_best_of_swarm(void)
 {
  int i,j;   
  for(i=0;i<pop_size;i++)
   {
    if(p[i].fitness>p_best.fitness) 
     {
       p_best.fitness=p[i].fitness;
       p_best.objective=p[i].objective;
       dPSO.t_best = clock() - dPSO.ts;
       dPSO.fe_best = dPSO.fe_curr;
	   	      
       for(j=0;j<num_dim;j++)
        {
          p_best.x[j]=p[i].x[j];
        }
      }
   }
 }

 void PSO::calc_speed_place(void)
 {
  int i,j,k,neig_best;
  int n1,n2; // neighborhood bounds
  double r0,r1,r2,r3;
  double neig_fitness=-1e300; 
 
  for(i=0;i<pop_size;i++)
   {
    n1=i;
    n2=i;
    neig_best=i;

    for(k=0;k<neighborhood;k++)//setting neighborhood bounds
     {
       if(n1-1>=0) n1--; else {if(n2+1<pop_size) n2++;}
       if(n2+1<pop_size) n2++; else {if(n1-1>=0) n1--;}
     }

    neig_fitness=-1e300; 

    if(neighborhood==0)
     {
       neig_fitness=p[i].fitness;
       neig_best=i;
     }
     else
     {
      for(k=n1;k<=n2;k++) // finding best particle in neighborhood
       {
         if(p[k].fitness>neig_fitness && k!=i) 
          {
            neig_fitness=p[k].fitness;
            neig_best=k;
          }
       }
     }
       
   for(j=0;j<num_dim;j++)
    {
     r0=RandVal(0.0,1.0); 
     r1=RandVal(0.0,1.0);
     r2=RandVal(0.0,1.0);
     r3=RandVal(0.0,1.0);        
 
     p[i].v[j] = c0*r0*p[i].v[j] + c2*r1*(p_best.x[j]-p[i].x[j]) + c1*r2*(p[i].x_best[j]-p[i].x[j]) + c3*r3*(p[neig_best].x[j]-p[i].x[j]);
     
     p[i].x[j]=p[i].x[j]+p[i].v[j];

     if(p[i].x[j]<pro_lim[j].min) p[i].x[j]=p[i].x[j]-p[i].v[j]; //keeping the particle in the limits of the problem being optimized
     if(p[i].x[j]>pro_lim[j].max) p[i].x[j]=p[i].x[j]-p[i].v[j];

     if(p[i].x[j]<pro_lim[j].min) p[i].x[j]=pro_lim[j].min;
     if(p[i].x[j]>pro_lim[j].max) p[i].x[j]=pro_lim[j].max;

     if(pro_lim[j].tn==0) p[i].x[j]=0;
    }

   }

 }


 void PSO::make_particle_copy(void)
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

 void PSO::danger(void)
 {
  //danger simulation - in case of danger particles escape/scatter in all directions (random directions)
  int i,j;   
  for(i=0;i<pop_size;i++)
   {
   for(j=0;j<num_dim;j++)    
    {
     if(pro_lim[j].tn==1) 
     p[i].v[j]=p[i].x[j]-RandVal(pro_lim[j].min,pro_lim[j].max);
     else
     p[i].v[j]=0;
    }
   }
 }







// Main function of swarm simulation

int PSO_simulation(void)
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


  srand(dPSO.p_rand);


  
  dPSO.ts = clock();
  dPSO.fe_curr=0;

  PSO pso;


  fprintf(galog,"%i\n",pso.pop_size);
  fprintf(galog,"%i\n",pso.num_iter);
  fprintf(galog,"%f|%f|%f|%f\n",pso.c0,pso.c1,pso.c2,pso.c3);
  fprintf(galog,"%i|%i| %i\n",pso.iter_danger,pso.neighborhood, pso.num_dim ); 

  z=0;


  pso.init_swarm();
 
  pso.evaluate_swarm();
  pso.remember_best_of_swarm();
  // pso.make_particle_copy(); not used
  
  
  printf("\n Starting simulation \n\n");


  
 for(k=0;k<pso.num_iter;k++)
  {

    pso.calc_speed_place();     
    pso.evaluate_swarm();
    pso.remember_best_of_swarm();
    // pso.make_particle_copy(); // not used
    if(z==pso.iter_danger) { pso.danger(); z=0;} else z++;

    sum=0.0;
    sum_square=0.0;

    for(i=0;i<pso.pop_size;i++)
     {
      sum=sum+pso.p[i].objective;
      sum_square+=pso.p[i].objective*pso.p[i].objective;
     }

    avg=sum/(double)pso.pop_size;
    square_sum=sum*sum/(double)pso.pop_size;
    stddev=sqrt((1.0/(double)(pso.pop_size-1))*(sqrt((sum_square-square_sum)*(sum_square-square_sum))));

    dPSO.t_all = clock() - dPSO.ts;


    fprintf(galog,"%d,%9e,%9e,%9e,%f,%li\n",k,pso.p_best.objective,sum,stddev,(((double)dPSO.t_all)/CLOCKS_PER_SEC),dPSO.fe_curr);

    status++;
   
    if(status==100)
     {
          status=0;
          // update status on dialog                  
     }
     
     printf("\r %i z %i (%i)",k+1,pso.num_iter, dPSO.sim_num);

 } //end of for(k=0;k<pso.num_iter;k++)

   fprintf(galog,"\n\n PSO simulation finished\n");
  
   fprintf(galog,"\nBest Fitness=%e,",pso.p_best.objective);
   for(j=0;j<pso.num_dim;j++)
    {
      fprintf(galog,"x[%i]=%f|,",j,pso.p_best.x[j]);
    }

   dPSO.best_objective=pso.p_best.objective;

   for(i=0;i<pso.num_dim; i++)
   {
    Best_solution[i]=pso.p_best.x[i];
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



