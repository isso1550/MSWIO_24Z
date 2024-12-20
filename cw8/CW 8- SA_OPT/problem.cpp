// SA algorithm implementation (c) Dariusz Baczyñski - Warsaw University Of Technology
// 
// 2020-11-27


#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include "problem.h"
#include "sa.h"




int    NVARS;



// Functions distinctive for problem solved

problem_data p_d[200]; // not used
problem_limits z_z [1000];
double Best_solution[1000]; /* Data of best solution */


// Reading data for problem solved
void data_read(void)
{
FILE *ga;
int i;



// Reading additional data concerning the problem 



// Bounds reading 
  if((ga=fopen("PROTOTYP.GOI","r"))==NULL)
    {
     printf("\n Can't open input file PROTOTYP.GOI for reading");
     exit(1);
    }

  fscanf(ga,"%i\n",&NVARS);

  for(i=0;i<NVARS;i++)
   {
    fscanf(ga,"%i,%f,%f\n",&z_z[i].tn,&z_z[i].min,&z_z[i].max);
   }

  fclose(ga);  

return;
}
