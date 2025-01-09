// HS algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-18

#ifndef _HS_H_
#define _HS_H_



#include <stdio.h>      /* printf, fopen */
#include <ctime>

typedef struct 
        {
            int tn;
            float min;
            float max;
        }problem_limits ;
          
          

struct HS_data
{

    float acc_rate; // accept rate
	float pa_rate; // pitch (note) adjustement rate 
    float pab; // pitch adjustment coeficient - the value to divide the range of n-th dimension and have the pitch adjustment range

    long pop_size;     // number of harmonies
    long num_iter;     // number of iterations

    long num_dim;      // number of problem dimensions

    long p_rand;       // starting point of pseudo random generator
    double best_objective; // best objective found
    int sim_num; //number of simulation calculated currently (from auto file)
    clock_t t_best,t_all; //time of geting the baest solution and all iterations
    clock_t ts; //time of starting 
    long fe_best, fe_curr; //number of function evaluations for the best and all current 
};
          

/***************************  harmony ********************************************************/          
    
class harmony
 {
  public:
  int dim;                // number of problem dimensions

  double *x;              // notes - location of solution in n-dimensional space  - number of harmony notes

  double fitness;         // harmony fitness - the function which will be optimised (max)
  double objective;       // harmony objective - aim function, problem criterium  

  harmony(void);         // object constructor
  ~harmony();            // object destructor

  void init(int dim, problem_limits *);     // function for harmony data initialisation 
  void calc_fitness(void); // function for harmony fitness calculation 
};


/************************************************ HS - Harmony Search  ********************************/

class HS
{
public:    
harmony *p;         // harmonies table
harmony *p_copy;    // harmonies table copy
harmony p_best;     // best harmony
problem_limits *pro_lim; // problem limits - limits of variables in each dimension

double acc_rate; // accept rate
double pa_rate; // pitch (note) adjustement rate 
double pab; // pitch adjustment coeficient - the value to divide the range of n-th dimension and have the pitch adjustment range

int num_dim;         
int pop_size;        // number of harmonies at each iteration of algorithm
int num_iter;        // number of iterations
int neighborhood;      
int iter_danger;     

HS(void);
~HS();


void init_harmonies(void);   
void evaluate_harmonies(void);  
void remember_best_of_harmonies(void);  
void generate_new_harmonies(void);

void make_harmony_copy(void); 
virtual void set_problem_limits(void);
};



int HS_simulation(void);          
          
double RandVal(double low,double high);
long iRandVal(long low,long high);
          
          
          
          
#endif          
