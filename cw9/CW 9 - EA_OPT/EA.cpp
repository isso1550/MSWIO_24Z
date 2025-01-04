// EA algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-04

#include <math.h>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "ea.H"
#include "problem.h"

extern int    NTRANSF;

extern double Best_solution[30]; /* Data of best solution */

 EA_data dEA;

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
     
    fitness=0.0;
    objective=0.0;

   // genome initialization   
   for(i=0;i<dim;i++) 
    {
     if(v_l[i].tn==1) 
     x[i]=RandVal(v_l[i].min,v_l[i].max);
     else
     x[i]=0;
    }

 }


     

 EA::~EA()
 {
     delete [] pop;
     delete [] pop2;
 }


 EA::EA(void)
 {
     int i;
     
    pop_size=dEA.pop_size;     // number of individuals
    num_iter=dEA.num_iter;     // number of iterations MAXGENS
    px=dEA.px;          // probability of crossover
    pm=dEA.pm;          // probability of mutation

    scal=dEA.scal;           // is fitness function scaling on
    elit=dEA.elit;           // is elitist strategy on
    
    aryt_xo=dEA.aryt_xo;        // is arythmetic crossover on

    num_dim=dEA.num_dim;      // number of problem dimensions
    fmultiple=2.0;
     
     pop=new individual[pop_size];           if(pop==NULL)  { printf("Memory allocation error 004"); exit(0);}
     pop2=new individual[pop_size];          if(pop2==NULL)  { printf("Memory allocation error 005"); exit(0);}
     pro_lim=new problem_limits[num_dim];    if(pro_lim==NULL)  { printf("Memory allocation error 006"); exit(0);}
     set_problem_limits();
 }     


void EA::init_population(void)
 {
   int i;
   
     for(i=0;i<pop_size;i++)
      {
          pop[i].init(num_dim, pro_lim);
          pop2[i].init(num_dim, pro_lim);
      }
     p_best.init(num_dim, pro_lim);

  p_best.fitness=-1e300; 
  p_best.objective=-1234; 
 }



 void EA::evaluate_population(void)
 {
  int i; 
     
  for(i=0;i<pop_size;i++)
   {
    pop[i].calc_fitness();
    if(pop[i].fitness<0.0) { printf("\n ****** Error - fitness function lower than 0.0");    }
   }
 }



 void EA::remember_best_of_population(void)
 {
  int i,j;  
  double previous_best;


   u_avg=0;
   u_min=1e300;
   u_max=-1e300;

  for(i=0;i<pop_size;i++)
   {

	u_avg=u_avg+pop[i].fitness;
    
    if(pop[i].fitness<u_min)
     {
     	u_min=pop[i].fitness;
     	Worst=i;
     }
    
    if(pop[i].fitness>u_max)
     {
     	u_max=pop[i].fitness;
     	Best=i;
     }
   }

   u_avg=u_avg/pop_size;     

   if(pop[Best].fitness>p_best.fitness)   best_survived=1;  // new best is better than old one, so we can assume that the best have survived ;-)
                                                            // check how it works - it is not the must 
         
   // in case of elit strategy make sure that the best individual survived 
   if(elit==1 && best_survived==0) 
    {
       pop[Worst].fitness  =p_best.fitness;
       pop[Worst].objective=p_best.objective;

       for(j=0;j<num_dim;j++)
        {
          pop[Worst].x[j]=p_best.x[j];
        }
    }
   
   if(pop[Best].fitness>p_best.fitness)
    {
       p_best.fitness=  pop[Best].fitness;
       p_best.objective=pop[Best].objective;
       dEA.t_best = clock() - dEA.ts;
       dEA.fe_best = dEA.fe_curr;
	      
       for(j=0;j<num_dim;j++)
        {
          p_best.x[j]=pop[Best].x[j];
        }
    }
 }


/****************************** scale ****************************************/
void EA::scale_fitness (void)
{
 double delta,a,b;
 int mem;

 if(scal==1)
  {

    if(u_min>((fmultiple*u_avg-u_max)/(fmultiple-1)))
     {
      delta=u_max-u_avg;
      a=(fmultiple-1)*u_avg/delta;
      b=u_avg*(u_max-fmultiple*u_avg)/delta;
     }
    else
     {
      delta=u_avg-u_min;
      a=u_avg/delta;
      b=-u_min*u_avg/delta;
     }

   for(mem=0;mem<pop_size;mem++)
    {
     pop[mem].fitness=a*pop[mem].fitness+b;
    }
    
  }

}


/****************************** select ***************************************/
void EA::select(void)
{
  int mem,i,j,bm;
  double sum,r;

  best_survived=0;

  sum=0.0;

  for(mem=0;mem<pop_size;mem++)
  {
    sum+=pop[mem].fitness;
  }

  for(mem=0;mem<pop_size;mem++)
  {
   pop[mem].rfitness=pop[mem].fitness/sum;
  }

  pop[0].cfitness=pop[0].rfitness;

  for(mem=1;mem<pop_size;mem++)
  {
   pop[mem].cfitness=pop[mem-1].cfitness + pop[mem].rfitness;
  }

  for(mem=0;mem<pop_size;mem++)
  {
   r=RandVal(0.0,1.0);

   for(i=0;i<pop_size;i++)
    {
     if( (i>0 && pop[i].cfitness>=r && pop[i-1].cfitness<r) ||
             (i==0 && r<=pop[i].cfitness) )
      {
       if(i==Best)
        {
         bm=mem;
         best_survived=1;
        }
       for(j=0;j<num_dim;j++) pop2[mem].x[j]=pop[i].x[j];
      }
    }
  }

Best=bm;

   for(i=0;i<pop_size;i++)
    {
       for(j=0;j<num_dim;j++) pop[i].x[j]=pop2[i].x[j];
    }
}


/************************* crossover *****************************************/

void EA::crossover(void)
{
  int i,count;
  int mem;
  int a;
  int b;
  int point;
  double xnum;
  double fa=0.5f;
  

  xnum=pop_size*px;

  for(mem=0;mem<xnum;mem++)
  {

    a=(rand() % pop_size);
    b=(rand() % pop_size);

    if(a==Best || b==Best) best_survived=0;

    point = (rand() %num_dim)+1;

    if(aryt_xo==0)
      {
        for(i=0;i<num_dim;i++)
        {
         pop2[0].x[i]=pop[a].x[i];
        }
    
        for(i=0;i<num_dim;i++)
        {
         pop2[1].x[i]=pop[b].x[i];
        }
    
        for(i=point;i<num_dim;i++)
        {
         pop[a].x[i]=pop2[1].x[i];
        }
    
        for(i=point;i<num_dim;i++)
        {
         pop[b].x[i]=pop2[0].x[i];
        }
      }
      else
      {
        for(i=0;i<num_dim;i++)
        {
         pop2[0].x[i]=pop[a].x[i];
        }
    
        for(i=0;i<num_dim;i++)
        {
         pop2[1].x[i]=pop[b].x[i];
        }
    
        for(i=0;i<num_dim;i++)
        {
         pop[a].x[i]= fa*pop2[0].x[i] + (1-fa)*pop2[1].x[i];
        }
    
        for(i=0;i<num_dim;i++)
        {
         pop[b].x[i]= fa*pop2[1].x[i] + (1-fa)*pop2[0].x[i];
        }
      }    
  }
}



/*********************** mutate ********************************/
void EA::mutate (void)
{
  int i;
  int nmutations;
  int member;
  int var;

  nmutations=(int)((double)(pop_size*num_dim)*pm);

  for(i=0;i<nmutations;i++)
  {
   member=rand()%pop_size;
   var=rand()%num_dim;
   if(pro_lim[var].tn==1) 
     pop[member].x[var]=RandVal(pro_lim[var].min,pro_lim[var].max);
   else
     pop[member].x[var]=0;
   if(member==Best) best_survived=0;
  }
}















// Main function of evolution simulation

int EA_simulation(void)
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


  srand(dEA.p_rand);

  dEA.ts = clock();
  dEA.fe_curr=0;


  EA ea;


  fprintf(galog,"%i\n",ea.pop_size);
  fprintf(galog,"%i\n",ea.num_iter);
  fprintf(galog,"%f|%f|%i|%i\n",ea.pm,ea.px,ea.scal,ea.elit);
  fprintf(galog,"%i| %i|\n",ea.aryt_xo, ea.num_dim ); 

  z=0;


  ea.init_population();
  ea.evaluate_population();
  ea.remember_best_of_population();
  ea.scale_fitness();
  
  printf("\n Starting simulation \n\n");


  
 for(k=0;k<ea.num_iter;k++)
  {
   
    ea.select();
    ea.crossover();
    ea.mutate();
    ea.evaluate_population();
    ea.remember_best_of_population();
    ea.scale_fitness();
    
    sum=0.0;
    sum_square=0.0;

    for(i=0;i<ea.pop_size;i++)
     {
      sum=sum+ea.pop[i].objective;
      sum_square+=ea.pop[i].objective*ea.pop[i].objective;
     }

    avg=sum/(double)ea.pop_size;
    square_sum=sum*sum/(double)ea.pop_size;
    stddev=sqrt((1.0/(double)(ea.pop_size-1))*(sqrt((sum_square-square_sum)*(sum_square-square_sum))));

    dEA.t_all = clock() - dEA.ts;


    fprintf(galog,"%d,%9e,%9e,%9e,%f,%li\n",k,ea.p_best.objective,sum,stddev,(((double)dEA.t_all)/CLOCKS_PER_SEC),dEA.fe_curr);

    status++;
   
    if(status==100)
     {
          status=0;
          // update status on dialog                  
     }
     
     printf("\r %i z %i (%i)",k+1,ea.num_iter, dEA.sim_num);

 } //end of for(k=0;k<ea.num_iter;k++)

   fprintf(galog,"\n\n EA simulation finished\n");
  
   fprintf(galog,"\nBest Fitness=%e,",ea.p_best.objective);
   for(j=0;j<ea.num_dim;j++)
    {
      fprintf(galog,"x[%i]=%e|,",j,ea.p_best.x[j]);
    }

   dEA.best_objective=ea.p_best.objective;

   for(i=0;i<ea.num_dim; i++)
   {
    Best_solution[i]=ea.p_best.x[i];
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



