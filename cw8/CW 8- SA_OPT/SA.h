// SA algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-11-27

#ifndef _SA_H_
#define _SA_H_



#include <stdio.h>      /* printf, fopen */
#include <ctime>



typedef struct  //structure to make a table with problem limits for each problem variable
        {
            int tn; // indicator if variable is used
            float min;
            float max;
        }problem_limits ;
          
          

struct SA_data
{

    long num_iter;     // number of iterations - in EA - MAXGENS
    long pop_size;
    
    float start_temp;  // temperature at the start of annealing process
    float gamma;       // gamma coeficient
    float boltz;       // Boltzman's coeficient

    long num_dim;      // number of problem dimensions

    long p_rand;       // starting point of pseudo random generator
    double best_objective; // best objective found
    int sim_num; //number of simulation calculated currently (from auto file)
    clock_t t_best,t_all; //time of geting the baest solution and all iterations
    clock_t ts; //time of starting 
    long fe_best, fe_curr; //number of function evaluations for the best and all current 
};
          

/***************************  individual - representing solution to the problem solved ********************************************************/          
    
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


/************************************************ SA - Simulated Annealing  ********************************/

class SA
{
public:    
individual *pop;         // population table - in case of SA - number of individuals is equal to 1
individual *pop2;        // population 2 table - copy
individual p_best;     // best individual from population - to remember the best solution
problem_limits *pro_lim; // problem limits - limits of variables in each dimension

long pop_size;     // number of individuals
long num_iter;     // number of iterations MAXGENS

float start_temp;  // temperature at the start of annealing process
float gamma;       // gamma coeficient
float boltz;       // Boltzman's coeficient
float Temp;		   // current temperature

    
long num_dim;      // number of problem dimensions

SA(void);
~SA();
void init_population(void);   
void evaluate_population(void);  
void remember_best_of_population(void);  
void create_new_solution(void);
int  decide(void);

virtual void set_problem_limits(void);
};



int SA_simulation(void);          
          
double RandVal(double low,double high); // simple functions to obtain random numbers
long iRandVal(long low,long high);
          
          
          
          
#endif          
