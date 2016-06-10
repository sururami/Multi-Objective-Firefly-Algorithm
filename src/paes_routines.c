#include "global.h"

void evaluate(sol *s, char *problem)
{
  if(s->to_evaluate == 1)
    {
      if (!strcmp(problem, "CHARMm27"))
	{
	  CHARMm27(s);
	  num_evaluations++;
	}
      else
	{
	  printf("Invalid argument for test problem given. Exiting.\n");
	  exit(-1);
	}
    }
}

void add_to_archive(sol *s)
{
  arc[arclength] = *s;
  arclength++;
}

int compare_to_archive(sol *s)  // compares a solution to every member of the archive. Returns -1 if dominated by
{                               // any member, 1 if dominates any member, and 0 otherwise
  int i=0;
  int result=0;

  while((i<arclength)&&(result!=1)&&(result!=-1))
    {
      //MINIMIZE MAXIMIZE
      if (minmax==0)
	result = compare_min(m->obj, (&arc[i])->obj, objectives);
      else
	result = compare_max(m->obj, (&arc[i])->obj, objectives);
      i++;
    }

  return(result);
} 

int compare_min(double *first, double *second, int n)
{
  // compares two n-dimensional vectors of objective values for minimization problems
  // returns 1 if first dominates second,
  // -1 if second dominates first, and 0 otherwise

  int obj = 0;
  int deflt = 0;
  int current;

  do
    {

      if(*first < *second)
	current = 1;
      else if(*second < *first)
	current = -1;
      else
	current = 0;
      if((current)&&(current==-deflt))
	{
	  return(0);
	}
      if(current!=0)
	{
	  deflt = current;
	}
      obj++;
      *first++;
      *second++;
    }while(obj < n);

  return(deflt);
}

int compare_max(double *first, double *second, int n)
{
  // as for compare_min() but for maximization problems

  int obj = 0;
  int deflt = 0;
  int current;

  do
    {
      if(*first > *second)
	current = 1;
      else if(*second > *first)
	current = -1;
      else
	current = 0;
      if((current)&&(current==-deflt))
	{
	  return(0);
	}
      if(current!=0)
	{
	  deflt = current;
	}
      obj++;
      *first++;
      *second++;
    }while(obj < n);
  return(deflt);
}



int equal(double *first, double *second, int n)
{
  // checks to n-dimensional vectors of objectives to see if they are identical
  // returns 1 if they are, 0 otherwise

  int obj = 0;

  do
    {
      if(*first!=*second)
	return(0);
      *first++;
      *second++;
      obj++;
      // printf("%d\n",obj);
    }
  while(obj < n);
  return(1);
}


void archive_soln(sol *s)
{
  // given a solution s, add it to the archive if
  // a) the archive is empty
  // b) the archive is not full and s is not dominated or equal to anything currently in the archive
  // c) s dominates anything in the archive
  // d) the archive is full but s is nondominated and is in a no more crowded square than at least one solution
  // in addition, maintain the archive such that all solutions are nondominated.

  int i;
  int repl;
  int yes = 0;
  int most;
  int result;
  int join = 0;
  int old_arclength;
  double *evs;
  double *evli;
  int set = 0;


  int tag[MAX_ARC];
  sol *tmp;
  if (!(tmp = (sol *)malloc(MAX_ARC*sizeof(sol))))
    {
      printf("Out of memory\n");
      exit(-1);
    }


  for (i = 0; i < archive; i++)
    {
      tag[i]=0;
    }


  if (arclength == 0)
    {
      add_to_archive(s);
      return;
    }


  i = 0;
  result = 0;
  while((i < arclength)&&(result!=-1))
    {
      result = equal(s->obj, (&arc[i])->obj, objectives);
      if (result == 1)
	break;
      //MINIMIZE MAXIMIZE
      if (minmax==0)
	result = compare_min(s->obj, (&arc[i])->obj, objectives);
      else
	result = compare_max(s->obj, (&arc[i])->obj, objectives);
      //  printf("%d\n", result);

      if ((result == 1)&&(join == 0))
	{
	  arc[i] = *s;
	  join = 1;
	}
      else if (result == 1)
	{
	  tag[i]=1;
	  set = 1;
	}
      i++;
    }

  old_arclength = arclength;
  if (set==1)
    {
      for (i = 0; i < arclength; i++)
	{
	  tmp[i] = arc[i];
	}
      arclength = 0;

      for (i = 0; i < old_arclength; i++)
	{
	  if (tag[i]!=1)
	    {
	      arc[arclength]=tmp[i];
	      arclength++;
	    }
	}
    }


  if ((join==0)&&(result==0))  // ie solution is non-dominated by the list
    {
      if (arclength == archive)
	{
	  most = grid_pop[s->grid_loc];
	  for (i = 0; i < arclength; i++)
	    {
	      if (grid_pop[(&arc[i])->grid_loc] > most)
		{
		  most = grid_pop[(&arc[i])->grid_loc];
		  repl = i;
		  yes = 1;
		  //   printf("i = %d\n", i);
		}
	    }
	  if (yes)
	    {
	      arc[repl] = *s;
	    }
	}
      else
	{
	  add_to_archive(s);
	}
    }
  free(tmp);
}

int find_loc(double *eval)
{
  // find the grid location of a solution given a vector of its objective values

  int loc = 0;
  int d;
  int n = 1;

  int i;

  int inc[MAX_OBJ];
  double width[MAX_OBJ];

  // printf("obj = %d, depth = %d\n", objectives, depth);


  // if the solution is out of range on any objective, return 1 more than the maximum possible grid location number
  for (i = 0; i < objectives; i++)
    {
      if ((eval[i] < gl_offset[i])||(eval[i] > gl_offset[i] + gl_range[i]))
	return((int)pow(2,(objectives*depth)));
    }

  for (i = 0; i < objectives; i++)
    {
      inc[i] = n;
      n *=2;
      width[i] = gl_range[i];
    }

  for (d = 1; d <= depth; d++)
    {
      for (i = 0; i < objectives; i++)
	{
	  if(eval[i] < width[i]/2+gl_offset[i])
	    loc += inc[i];
	  else
	    gl_offset[i] += width[i]/2;
	}
      for (i = 0; i < objectives; i++)
	{
	  inc[i] *= (objectives *2);
	  width[i] /= 2;
	}
    }
  return(loc);
}

void update_grid(sol *s)
{
  // recalculate ranges for grid in the light of a new solution s
  static int change = 0;
  int a, b;
  int square;
  double offset[MAX_OBJ];
  double largest[MAX_OBJ];
  double sse;
  double product;
  struct solution **start;

  for (a = 0; a < objectives; a++)
    {
      offset[a] = LARGE;
      largest[a] = -LARGE;
    }

  for (b = 0; b < objectives; b++)
    {
      for (a = 0; a < arclength; a++)
	{
	  if ((&arc[a])->obj[b] < offset[b])
	    offset[b] = (&arc[a])->obj[b];
	  if ((&arc[a])->obj[b] > largest[b])
	    largest[b] = (&arc[a])->obj[b];
	}
    }
//    printf("oldCURENT:largest = %f, offset = %f\n", largest[0], offset[0]);
//    printf("oldCURENT:largest = %f, offset = %f\n", largest[1], offset[1]);


  for (b = 0; b < objectives; b++)
    {
      if (s->obj[b] < offset[b])
	offset[b] = s->obj[b];
      if (s->obj[b] > largest[b])
	largest[b] = s->obj[b];
    }

  sse = 0;
  product = 1;

  for (a = 0; a < objectives; a++)
    {

      sse += ((gl_offset[a] - offset[a])*(gl_offset[a] - offset[a]));
      sse += ((gl_largest[a] - largest[a])*(gl_largest[a] - largest[a]));
      product *= gl_range[a];
    }

  // printf("sse = %f\n", sse);
  if (sse > (0.1 * product * product))	//if the summed squared error (difference) between old and new
                                        //minima and maxima in each of the objectives
    {                                   //is bigger than 10 percent of the square of the size of the space
      change++;                         // then renormalise the space and recalculte grid locations

      for (a = 0; a < objectives; a++)
	{
	  gl_largest[a] = largest[a]+0.2*largest[a];
	  gl_offset[a] = offset[a]+0.2*offset[a];
	  gl_range[a] = gl_largest[a] - gl_offset[a];
	}

      for (a = 0; a < pow(2, (objectives*depth)); a++)
	{
	  grid_pop[a] = 0;
	}

      for (a = 0; a < arclength; a++)
	{
	  square = find_loc((&arc[a])->obj);
	  (&arc[a])->grid_loc = square;
	  grid_pop[square]++;

	}
    }
  square = find_loc(s->obj);
  s->grid_loc = square;
  grid_pop[(int)pow(2,(objectives*depth))] = -5;
  grid_pop[square]++;

}
