#include "global.h"


//Auxiliary program for text processing
char* eat_space(char *cptr)
{ 
  while(*cptr == ' ')
  {
    if(*cptr == '\0')
      return NULL;
    else
      *cptr++;
  }
  return cptr;
}

//Auxiliary program for text processing  
char* next_next_num(char *cptr)
{

  while(*cptr != ' ')
  {
    if(*cptr == '\0')
      return NULL;
    else
      cptr++;
  }
    
  return eat_space(cptr);
}

// Give an list of solution, return the solution with min energy (2-objective sum)
int min_cl(sol *cls)
{
  int i,j;
  int best = 0;
  double best_sum_obj;
  double curr_sum_obj;
  for(i=1;i<MAX_CLONES;i++)
    {
      best_sum_obj = 0.0;
      curr_sum_obj = 0.0;
      for (j = 0; j < objectives; j++)
	{
	  best_sum_obj += cls[best].obj[j];
	  curr_sum_obj += cls[i].obj[j];
	}
      if(curr_sum_obj < best_sum_obj)
	best = i;
    }
  return best;
}

// Give an list of solution, return the solution with max energy (2-objective sum) 
int max_cl(sol *cls)
{
  int i,j;
  int best = 0;
  double best_sum_obj;
  double curr_sum_obj;
  for(i=1;i<MAX_CLONES;i++)
    {
      best_sum_obj = 0.0;
      curr_sum_obj = 0.0;
      for (j = 0; j < objectives; j++)
	{
	  best_sum_obj += cls[best].obj[j];
	  curr_sum_obj += cls[i].obj[j];
	}
      if(curr_sum_obj > best_sum_obj)
	best = i;
    }
  return best;
} 
