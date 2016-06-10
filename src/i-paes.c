// I-PAES skeleton program code for the Protein Structure Prediction problem (PSP)
/*
  Copyright (C) 2005  Vincenzo Cutello, Giuseppe Narzisi and Giuseppe Nicosia.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

   The GNU General Public License is available at:
      http://www.gnu.org/copyleft/gpl.html 

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/
//
// I-PAES is described in the paper:
//
//    Cutello V., Narzisi G. , Nicosia G.:
//    A Multi-Objective Evolutionary Approach to the Protein Structure Prediction Problem.
//    Journal of the Royal Society Interface, Royal Society Publications London.
//
// Please contact the authors if you have any comments, suggestions or questions
// about this file or the I-PAES algorithm.
// We are at:
//             {vctl,narzisi,nicosia}@dmi.unict.it
//
/*********************************************************************************************************
 * To run :                                                                                              *
 *  ./i-paes [Energy] [depth] [objectives] [residues] [archive] [iterations] [Tmax] [protein file]       *
 *           [dir] [SS]                                                                                  *
 *                                                                                                       *
 * where all parameters MUST be specified correctly (none are optional) following the instructions below:*
 *                                                                                                       *
 * [Energy] - at present there is only one Energy function included with this code: CHARMM (v. 27)       *
 * [depth] - this is the number of recursive subdivisions of the objective space carried out in order to *
 *           divide the objective space into a grid for the purposes of diversity maintenance. Values of *
 *           between 3 and 6 are useful.                                                                 *
 * [objectives] - the number of objectives to the problem (fixed to 2).                                  * 
 * [residues] - the number of residues in the protein, (only integer numbers can be accepted).           *
 * [archive] - the maximum number of solutions to be held in the nondominated solutions archive,         *
 *             (suggested values are in the range [100,1000]).                                           *
 * [iterations] - the program terminates after this number of iterations.                                *
 * [Tmax] - the program terminates after this number of fitness function evaluations.                    *
 * [minmax] - set this to 0 (zero) for minimization problems, 1 for maximization problems (for PSP 0)    *
 * [protein file] - input protein  instance file.                                                        *
 * [dir] - output directory for simulation results.                                                      *
 * [SS] - Secondary or Supersecondary Structure Constraints (1=secondary, 0=supersecondary)              *
 *                                                                                                       *
 * Example of a valid command lines for 1ZDD protein:                                                    *
 * ./i-paes CHARMm27 4 2 34 1000 250000 250000 0 instances/1ZDD.seq 1ZDD 1                                 *
 *                                                                                                       *
 ********************************************************************************************************/

#include "global.h"

FILE *fp;
int RandSeed;	//seed for the random number generator ran()

sol *FireflyArray; //All Solutions
sol *curr; // current solution
sol *cl; // clones solutions
sol *m; // mutant solution
sol *arc; // archive of solutions
sol *app;


int phi[8][2];
int psi[8][2];

int native_ss;

int seed, depth, genes, archive, objectives, iterations, Tmax, minmax; // command line parameters
int arclength=0; // current length of the archive
int num_evaluations;
int flag;
float sigma;

float gl_offset[MAX_OBJ];   // the range, offset etc of the grid
float gl_range[MAX_OBJ];
float gl_largest[MAX_OBJ];
int grid_pop[MAX_LOC];   // the array holding the population residing in each grid location


char problem[30];
char protein[200];
char dir[200];
char params[50];

double best_energy_found;

/**

MultiObjective Firefly Algorithm

*/

void MO_FA_GMJ(){
  int i, j, k, result;
  
   char M2;
  double M1;
  
  int generations=10000;
  int NumberFireflies = 1000;
  
  sol *FA, *FB, *FR;

    FA   = (sol *)malloc(MAX_POP*sizeof(sol));
    FB   = (sol *)malloc(NumberFireflies*sizeof(sol));
    FR   = (sol *)malloc(MAX_POP*sizeof(sol));
  
  for(i=1; i < generations; i++){
    for(j=0; j< NumberFireflies; j++){

  	copySolution(FA,&FireflyArray[k]);
    for(k=0; k< NumberFireflies;k++){
    	copySolution(FB,&FireflyArray[k]);
       
        if(FA->obj!=FB->obj){

          evaluate(FA, problem);
          evaluate(FB, problem);
          
          M1 = exp((-1)*((double)(2.0*num_evaluations)/(double)Tmax));
          M2 = 1 + (int)(((double)genes/(double)4) * exp(-(double)(2.0*num_evaluations)/(double)Tmax));
          
          result = compare_min(FA->obj, FB->obj, objectives);
        

          if(result==-1){ // FB Dominates FA ----------> Epsilon Dominating

              global_mutation(&FireflyArray[0], M1); // global mutation

              //local_mutation(&cl[1], M2); // local mutation
              
              evaluate(FR, problem);
              
              result = compare_min((&FB[0])->obj, (&FR[0])->obj, objectives);
              
              if(result==1) //FB dominate FR
                FireflyArray[k]=FR[0];
          }
        }       
      }      
    }    
  }  
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Initalization procedure
/*
void init()
{
  int i,j;
  char p;
  FILE *fdesc; 
	
  char name[4];
  float chi_mean[4];
  float chi_dev[4];
  int type;
	
  char buffer[160];

  num_evaluations = 0;
	 
  if (!strcmp(problem, "CHARMm27"))
	strcpy(params,"bin/params/charmm27.prm");
  else if (!strcmp(problem, "amber99"))
	strcpy(params,"bin/params/amber99.prm");
  

  // allocate memory for solutions
  curr   = (sol *)malloc(MAX_POP*sizeof(sol));
  cl  = (sol *)malloc(MAX_CLONES*sizeof(sol));
  m   = (sol *)malloc(MAX_POP*sizeof(sol));
  arc = (sol *)malloc(MAX_ARC*sizeof(sol));
  app = (sol *)malloc(MAX_ARC*sizeof(sol));

  if((!curr)||(!m)||(!arc))
    {
      printf("Out of memory. Aborting.\n");
      exit(-1);
    }
  
  if(native_ss == 1)
    {
      //Secondary structures constraints
      phi[0][0] = -67;  phi[0][1] = -47;  psi[0][0] = -57;  psi[0][1] = -37; // H
      phi[1][0] = -130; phi[1][1] = -110; psi[1][0] = 110;  psi[1][1] = 130; // B
      phi[2][0] = -130; phi[2][1] = -110; psi[2][0] = 110;  psi[2][1] = 130; // E
      phi[3][0] = -59;  phi[3][1] = -39;  psi[3][0] = -36;  psi[3][1] = -16; // G
      phi[4][0] = -67;  phi[4][1] = -47;  psi[4][0] = -80;  psi[4][1] = -60; // I
      phi[5][0] =  -180;phi[5][1] = 180;  psi[5][0] = -180;  psi[5][1] = 180;// T
      phi[6][0] =  -180;phi[6][1] = 180;  psi[6][0] = -180; psi[6][1] = 180; // S 
      phi[7][0] =  -180;phi[7][1] = 180;  psi[7][0] = -180; psi[7][1] = 180; // U (undefined)
    }
  else
    {
      //Supersecondary structures constraints
      phi[0][0] = -75;  phi[0][1] = -55;  psi[0][0] = -50;  psi[0][1] = -30; // H 
      phi[1][0] = -130; phi[1][1] = -110; psi[1][0] = 110;  psi[1][1] = 130; // E
      phi[2][0] = -150; phi[2][1] = -30;  psi[2][0] = -100; psi[2][1] = 50;  // a
      phi[3][0] = -230; phi[3][1] = -30;  psi[3][0] = 100;  psi[3][1] = 200; // b
      phi[4][0] =  30;  phi[4][1] = 130;  psi[4][0] = 130;  psi[4][1] = 260; // e
      phi[5][0] =  30;  phi[5][1] = 150;  psi[5][0] = -60;  psi[5][1] = 90;  // l
      phi[6][0] =  -160;phi[6][1] = -50;  psi[6][0] = 50;   psi[6][1] = 100; // t
      phi[7][0] =  -180;phi[7][1] = 180;  psi[7][0] = -180; psi[7][1] = 180; // U (undefined)
    }

  fdesc = fopen(protein,"r");
  if(fdesc == NULL)
  {
    printf("fopen failed to open %s\n", protein);
    exit(1);
  }
   printf("\n\n\n\n\nArchivo Abierto %s\n\n\n\n ", protein);
  for (i = 0; i < genes; i++)
    {
      fgets(buffer,159,fdesc);
      sscanf(buffer,"%s %c",curr->chrom[i].name, &p);
      if(native_ss == 1)
	curr->chrom[i].predicted = PredictedToNumber_native(p);
      else  
	curr->chrom[i].predicted = PredictedToNumber(p);
      curr->chrom[i].type = NameToType(curr->chrom[i].name);
      curr->chrom[i].num_angles = 2 + get_num_sidechain_angles(curr->chrom[i].name);
    }
  fclose(fdesc);
  
  curr->energy = 0.0;
  curr->to_evaluate = 1;

   
  // initialise curr by randomly selecting the backbone and sidechain torsion angles 
  // in the constrained regions.
  for (j = 0; j < genes; j++)
    {
      res* r = &curr->chrom[j];
      randConstAngles(r);
    }
}
*/

void init_MOFA(int N_Fireflies)
{
  int i,j;
  char p;
  FILE *fdesc; 
	
  char name[4];
  float chi_mean[4];
  float chi_dev[4];
  int type;
	


  char buffer[160];

  num_evaluations = 0;
	 
  if (!strcmp(problem, "CHARMm27"))
	strcpy(params,"bin/params/charmm27.prm");
  else if (!strcmp(problem, "amber99"))
	strcpy(params,"bin/params/amber99.prm");
  

  // allocate memory for solutions


  FireflyArray = (sol*) malloc(N_Fireflies * MAX_POP*sizeof(sol));

  curr   = (sol *)malloc(MAX_POP*sizeof(sol));
  cl  = (sol *)malloc(MAX_CLONES*sizeof(sol));
  m   = (sol *)malloc(MAX_POP*sizeof(sol));
  arc = (sol *)malloc(MAX_ARC*sizeof(sol));
  app = (sol *)malloc(MAX_ARC*sizeof(sol));

  if((!curr)||(!m)||(!arc) || (!FireflyArray))
    {
      printf("Out of memory. Aborting.\n");
      exit(-1);
    }

  if(native_ss == 1)
    {
      //Secondary structures constraints
      phi[0][0] = -67;  phi[0][1] = -47;  psi[0][0] = -57;  psi[0][1] = -37; // H
      phi[1][0] = -130; phi[1][1] = -110; psi[1][0] = 110;  psi[1][1] = 130; // B
      phi[2][0] = -130; phi[2][1] = -110; psi[2][0] = 110;  psi[2][1] = 130; // E
      phi[3][0] = -59;  phi[3][1] = -39;  psi[3][0] = -36;  psi[3][1] = -16; // G
      phi[4][0] = -67;  phi[4][1] = -47;  psi[4][0] = -80;  psi[4][1] = -60; // I
      phi[5][0] =  -180;phi[5][1] = 180;  psi[5][0] = -180;  psi[5][1] = 180;// T
      phi[6][0] =  -180;phi[6][1] = 180;  psi[6][0] = -180; psi[6][1] = 180; // S 
      phi[7][0] =  -180;phi[7][1] = 180;  psi[7][0] = -180; psi[7][1] = 180; // U (undefined)
    }
  else
    {
      //Supersecondary structures constraints
      phi[0][0] = -75;  phi[0][1] = -55;  psi[0][0] = -50;  psi[0][1] = -30; // H 
      phi[1][0] = -130; phi[1][1] = -110; psi[1][0] = 110;  psi[1][1] = 130; // E
      phi[2][0] = -150; phi[2][1] = -30;  psi[2][0] = -100; psi[2][1] = 50;  // a
      phi[3][0] = -230; phi[3][1] = -30;  psi[3][0] = 100;  psi[3][1] = 200; // b
      phi[4][0] =  30;  phi[4][1] = 130;  psi[4][0] = 130;  psi[4][1] = 260; // e
      phi[5][0] =  30;  phi[5][1] = 150;  psi[5][0] = -60;  psi[5][1] = 90;  // l
      phi[6][0] =  -160;phi[6][1] = -50;  psi[6][0] = 50;   psi[6][1] = 100; // t
      phi[7][0] =  -180;phi[7][1] = 180;  psi[7][0] = -180; psi[7][1] = 180; // U (undefined)
    }
 int k;
 for( k=0; k < N_Fireflies*MAX_POP; k+=MAX_POP){  
    fdesc = fopen(protein,"r");
    if(fdesc == NULL)
    {
      printf("fopen filed to open %s\n", protein);
      exit(1);
    }
 
    for (i = 0; i < genes; i++)
    {
      fgets(buffer,159,fdesc);
      sscanf(buffer,"%s %c",curr->chrom[i].name, &p);
      if(native_ss == 1)
	curr->chrom[i].predicted = PredictedToNumber_native(p);
      else  
	curr->chrom[i].predicted = PredictedToNumber(p);
      curr->chrom[i].type = NameToType(curr->chrom[i].name);
      curr->chrom[i].num_angles = 2 + get_num_sidechain_angles(curr->chrom[i].name);
    }
  fclose(fdesc);
  
  curr->energy = 0.0;
  curr->to_evaluate = 1;

   
  // initialise curr by randomly selecting the backbone and sidechain torsion angles 
  // in the constrained regions.
  for (j = 0; j < genes; j++)
    {
      res* r = &curr->chrom[j];
      randConstAngles(r);
    }

  copySolution(&FireflyArray[k],curr);
  
 }
 
}

void copySolution(sol* A, sol *B){

    int i,j,k;
    for(i=0; i < genes; i++){

	A->chrom[i].name[0]=B->chrom[i].name[0];
	A->chrom[i].name[1]=B->chrom[i].name[1];
	A->chrom[i].name[2]=B->chrom[i].name[2];
	A->chrom[i].name[3]=B->chrom[i].name[3];

        A->chrom[i].type=B->chrom[i].type;
        for(j=0; j < MAX_ANGLES; j++){
            A->chrom[i].angles[j]=B->chrom[i].angles[j];
        }
        A->chrom[i].num_angles=B->chrom[i].num_angles;
        A->chrom[i].predicted=B->chrom[i].predicted;
    
        for(k=0; k < MAX_OBJ; k++){
            A->obj[j]=B->obj[j];
        }
        A->grid_loc=B->grid_loc;
        A->energy=B->energy;
        A->to_evaluate=B->to_evaluate;
    }
    
}
int main(int argc, char *argv[])
{

  int i, j, k, r, w, iteration;
  int result;
  char arch[200];  
  char command[200] = "mkdir ";
  char temp[20];
  FILE *fdesc;
  double tmp_obj;
	
  int M2;
  double M1;

  if (argc!=11)
    {
      printf("Usage: %s [problem] [depth] [objectives] [residues] [archive] [iterations] [Tmax] [dir] [SS]\n", argv[0]);
      exit(-1);
    }

  //get command line parameters

  sprintf(problem, "%s", argv[1]);
  depth = atoi(argv[2]);
  objectives = atoi(argv[3]);
  genes = atoi(argv[4]);
  archive = atoi(argv[5]);
  iterations = atoi(argv[6]);
  Tmax = atoi(argv[7]);
  minmax = 0;
  sprintf(protein, "%s", argv[8]);
  sprintf(dir,"%s",argv[9]);
  native_ss = atoi(argv[10]);

  std::vector<int> first;

  srandom(time(NULL)); // seed for function random() 
  RandSeed=random();   //   "    " "  " ran()
  srand(r=random());
  //printf("Seed = %d\n",r);

  RandomInitialise(Randint(0,31328),Randint(0,30081));
   
  sigma = 1.0;

  // directory for simulation results
  //strcat(dir,"_depth");
  //strcat(dir,gcvt(depth,20,temp));
  //strcat(dir,"_arch");
  //strcat(dir,gcvt(archive,20,temp));
  //strcat(dir,"_iter");
  //strcat(dir,gcvt(iterations,20,temp));
  strcat(command,dir);
  system(command);

  print_params("parameters.dat");

  // begin (1+1)-PAES
    
  signal(SIGINT,backup);// handler for program termination by keyboard (SIGINT signal)

  int ejemplos= 100;
  init_MOFA(ejemplos);
  omp_set_num_threads(100);
  double time1=omp_get_wtime();
  #pragma omp parallel for

  for(int p=0; p < ejemplos  ; p++  ){
	evaluate(&FireflyArray[p*MAX_POP], problem);
  	printf("\n\n\n Energia Actual: %f\n",(&FireflyArray[p*MAX_POP])->energy);
  } 
  double time2=omp_get_wtime();

	printf("Tiempo: %f\n", time2-time1);

/*

  evaluate(curr, problem);  // Fitness function evaluation

  //  print_eval(c);      // Uncomment to check objective values generated

  add_to_archive(curr);

  best_energy_found = curr->energy;

  // begin main loop
  i=0; //iterations counter
  flag = 1;

*/


/*

  while(flag)
    {
      if (i%100==0)  // just print out the number of iterations done every N
	printf("Iterations: %d/%d\tFFE: %d\n", i, iterations, num_evaluations);

      //if (i%1000==0) print_arc(i); // Uncomment for snapshot data saving every N iterations

      if(curr->energy < best_energy_found) best_energy_found = curr->energy;

      print_statistics("archive_statistics.dat"); // Uncomment for data saving of archive statistics
      print_curr_sol("current_solution_statistics.dat"); // Uncomment for data saving of current solution statistics
      
      // Immune Phase: Start 

      cl[0] = *curr;   // copy the current solution
      cl[1] = *curr;   // copy the current solution
      
      //Mutation

      M1 = exp((-1)*((double)(2.0*num_evaluations)/(double)Tmax));
      M2 = 1 + (int)(((double)genes/(double)4) * exp(-(double)(2.0*num_evaluations)/(double)Tmax));

      global_mutation(&cl[0], M1); // global mutation
      local_mutation(&cl[1], M2); // local mutation
      
      //Fitness evaluation

      evaluate(&cl[0], problem); //fitness evaluation for clone cl[0]
      evaluate(&cl[1], problem); //fitness evaluation for clone cl[1]

      //Selection
      if (minmax==0)
	result = compare_min((&cl[0])->obj, (&cl[1])->obj, objectives);
      else
	result = compare_max((&cl[0])->obj, (&cl[1])->obj, objectives);

      if (result==1)  // if cl[0] dominates cl[1]
	*m = cl[0];   // replace m with cl[0]
      else if (result==-1)  // if cl[1] dominates cl[0]
	*m = cl[1];         // replace m with cl[0]
      else{           // cl[0] and cl[1] are nondominated
	*m = cl[min_cl(cl)];
	
  update_grid(&cl[max_cl(cl)]);  //calculate grid location of cl[1] solution and renormalize archive if necessary
	archive_soln(&cl[max_cl(cl)]); //update the archive by removing all dominated individuals
      }
      // Immune Phase: End

      //MINIMIZE MAXIMIZE
      if (minmax==0)
	result = compare_min(curr->obj, m->obj, objectives);
      else
	result = compare_max(curr->obj, m->obj, objectives);

      // printf("RESULT = %d\n", result);
      // printf("arclength = %d\n", arclength);

      if (result != 1)  // if mutant is not dominated by current (else discard it)
	{
	  if (result ==-1)  // if mutant dominates current
	    {
	      //printf("m dominates c\n");
	      update_grid(m);           //calculate grid location of mutant solution and renormalize archive if necessary
		archive_soln(m);          //update the archive by removing all dominated individuals

	      *curr = *m;                    // replace c with m
	    }
	  else if(result == 0)  // if mutant and current are nondominated wrt each other
	  {
	      result = compare_to_archive(m);
	      if (result != -1)  // if mutant is not dominated by archive (else discard it)
		    {
          update_grid(m);
          archive_soln(m);
          
          if((grid_pop[m->grid_loc] <= grid_pop[curr->grid_loc])||(result==1)) // if mutant dominates the archive or
          {                                                                     // is in less crowded grid loc than c
              *curr = *m; // then replace c with m
          }
        }
	    }
	}
      i++;
      //Stop after Tmax evaluations or fixed number of iteration
      if( (num_evaluations >= Tmax) || (i>iterations) ) flag=0;

    }

*/

/*
  iteration = i;
  printf("\nThe Archive is now...\n");
  for (i = 0; i < arclength; i++)
    print_eval(&arc[i]);

  printf("\nSaving of the genetic material in the archive...\n");
  print_arc(iteration);
  printf("Done.\n");*/
}
