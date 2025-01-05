// PSO algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-12-12

#ifndef _PSO_H_
#define _PSO_H_



#include <stdio.h>      /* printf, fopen */
#include <ctime>

typedef struct 
        {
            int tn;
            float min;
            float max;
        }problem_limits ;
          
          

struct PSO_data
{

    float c0; // coef for current speed 
    float c1; // coef for best position of this particle (found so far)
    float c2; // coef for best position among all particles (found so far)
    float c3; // coef for best position among neighbor particles (found so far)
    long iter_danger;  // number of iterations between dangerous situation for particles(animals) - scattering/dispersion of swarm
    int  neighborhood; // neighborhood width

    long pop_size;     // number of particles
    long num_iter;     // number of iterations

    long num_dim;      // number of problem dimensions

    long p_rand;       // starting point of pseudo random generator
    double best_objective; // best objective found
    int sim_num; //number of simulation calculated currently (from auto file)
    clock_t t_best,t_all; //time of geting the baest solution and all iterations
    clock_t ts; //time of starting 
    long fe_best, fe_curr; //number of function evaluations for the best and all current 
};
          

/***************************  particle ********************************************************/          
    
class particle
 {
  public:
  int dim;                // number of problem dimensions

  double *x;              // location of particle in n-dimensional space  
  double *v;              // particle speed in n-dimensional space
  double *x_best;         // best position of organism in n-dimansional space

  double fitness;         // particle fitness - the function which will be optimised (max)
  double fitness_best;    // best fitness of particle found so far
  double objective;       // particle objective - aim function, problem criterium  

  particle(void);         // object constructor
  ~particle();            // object destructor

  void init(int dim, problem_limits *);     // function for particle data initialisation 
  void calc_fitness(void); // function for particle fitness calculation 
  void remember_best_p(void); // remembering of particle best position 
};


/************************************************ PSO - Particle Swarm Optimization  ********************************/

class PSO
{
public:    
particle *p;         // particles table
particle *p_copy;    // particles table copy
particle p_best;     // best particle
problem_limits *pro_lim; // problem limits - limits of variables in each dimension

double c0,c1,c2,c3;  
int num_dim;         
int pop_size;        
int num_iter;        
int neighborhood;      
int iter_danger;     

PSO(void);
~PSO();
void init_swarm(void);   
void evaluate_swarm(void);  
void remember_best_of_swarm(void);  
void calc_speed_place(void);
void make_particle_copy(void); 
void danger(void); // simulation of danger appearing in the swarm - swarm scattering
virtual void set_problem_limits(void);
};



int PSO_simulation(void);          
          
double RandVal(double low,double high);
long iRandVal(long low,long high);
          
          
          
          
#endif          
