// EA algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-04

#ifndef _EA_H_
#define _EA_H_



#include <stdio.h>      /* printf, fopen */
#include <ctime>

typedef struct 
        {
            int tn;
            float min;
            float max;
        }problem_limits ;
          
          

struct EA_data
{

    long pop_size;     // number of individuals
    long num_iter;     // number of iterations MAXGENS
    float px;          // probability of crossover
    float pm;          // probability of mutation

    int scal;           // is fitness function scaling on
    int elit;           // is elitist strategy on
    
    int aryt_xo;        // is arythmetic crossover on

    long num_dim;      // number of problem dimensions

    long p_rand;       // starting point of pseudo random generator
    double best_objective; // best objective found
    int sim_num; //number of simulation calculated currently (from auto file)
    clock_t t_best,t_all; //time of geting the baest solution and all iterations
    clock_t ts; //time of starting 
    long fe_best, fe_curr; //number of function evaluations for the best and all current 
};
          

/***************************  individual ********************************************************/          
    
class individual
 {
  public:
  int dim;                // number of problem dimensions

  double *x;              // individual genome
 
  double fitness;         // individual fitness - the function witch will be optimised (max)
  double objective;       // individual objective - aim function, problem criterium  
  
  double rfitness;        // relative fitness
  double cfitness;        // cumulative rfitness for select function

  individual(void);         // object constructor
  ~individual();            // object destructor

  void init(int dim, problem_limits *);     // function for individual data initialisation 
  void calc_fitness(void); // function for individual fitness calculation 
};


/************************************************ EA - Evolutionary Algorithm  ********************************/

class EA
{
public:    
individual *pop;         // population table
individual *pop2;        // population 2 table - copy
individual p_best;     // best individual from population
problem_limits *pro_lim; // problem limits - limits of variables in each dimension

long pop_size;     // number of individuals
long num_iter;     // number of iterations MAXGENS
float px;          // probability of crossover
float pm;          // probability of mutation

int scal;           // is fitness function scaling on
int elit;           // is elitist strategy on
    
int aryt_xo;        // is arythmetic crossover on
long num_dim;      // number of problem dimensions

int Best, Worst;  // best and worst member 
int best_survived; // flag if best survived during recombinantion nad selection

double u_avg,u_min,u_max;
double fmultiple;

EA(void);
~EA();
void init_population(void);   
void evaluate_population(void);  
void remember_best_of_population(void);  
void scale_fitness(void);
void select(void);
void crossover(void);
void mutate (void);

virtual void set_problem_limits(void);
};



int EA_simulation(void);          
          
double RandVal(double low,double high);
long iRandVal(long low,long high);
          
          
          
          
#endif          
