#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/wait.h>
#include <omp.h>

using namespace std;


//#include "ran.c"

/*----------Random number generator--------------*/
#include "random/randomlib.h"
//#include "random/randomlib.c"
/*----------------------------------------------*/

# define MAX_RESIDUES 200  // change as necessary
# define MAX_OBJ 2         // change as necessary
# define MAX_ARC 1000      // change as necessary
# define MAX_LOC 32768     // number of locations in grid (set for a three-objective problem using depth 5)
# define LARGE 2000000000  // should be about the maximum size of an integer for your compiler
# define MAX_ANGLES 6      // change as necessary
# define MAX_CLONES 2      // number of fixed clones 

/*------------- Defines for random number use --------------------------*/
/* Randint(L,H) produces an integer between L and H inclusive           */
/* Rand() produces a real in [0,1]                                      */
#define Randint(L,H) (((H) > (L)) ? (rand()%((H)-(L)+1)+ (L)) : (L))
#if 0
#define Rand() (rand()*pow(RAND_MAX,-1))
#else
#define Rand() (ran())
#endif
/*----------------------------------------------------------------------*/

# define COIN(p)  (RandomUniform() * 100 < p)

extern FILE *fp;

typedef struct residue
{
  char name[4];
  int type;
  float angles[MAX_ANGLES];
  int num_angles;
  int predicted;
}res;


typedef struct solution
{
  res chrom[MAX_RESIDUES];
  double obj[MAX_OBJ];
  double objNorm[MAX_OBJ];
  int grid_loc;
  double energy;
  short to_evaluate;
}sol;

extern int RandSeed;	//seme per il generatore di numeri casuali ran()

extern sol *curr; // current solution
extern sol *cl; // clones solutions
extern sol *m; // mutant solution
extern sol *arc; // archive of solutions
extern sol *app;

extern int phi[8][2];
extern int psi[8][2];

extern int native_ss;
int compare_min(double *first, double *second, int n);
int compare_max(double *first, double *second, int n);
int equal(double *first, double *second, int n);
void init();
void print_genome(sol *s, int n, const char* dir);
void print_eval(sol *s);
void evaluate(sol *s, char *problem);
void CHARMm27(sol *s);
void add_to_archive(sol *s);
int compare_to_archive(sol *s);
void update_grid(sol *s);
void archive_soln(sol *s);
int find_loc(double *eval);

float getSidechainAngles(char *name,int n);
void protein_file(const char *file, sol *s, const char *title, const char *filexyz);
void randConstAngles(res *r);
int get_num_sidechain_angles(char *res);
void print_archive(const char *file);
void print_arc(int n);
void print_statistics(const char* name);
void print_curr_sol(const char* name);
void print_params(const char *file);
void backup(int signo);
void local_mutation(sol *s, int n_m);
void global_mutation(sol *s, double pm);
int min_cl(sol *cls);
int max_cl(sol *cls);
int PredictedToNumber(char p);
int PredictedToNumber_native(char p);
void phi_constrain(float *angle, int p, char *name);
void psi_constrain(float *angle, int p, char *name);
void chi_constrain(float *angle, char* residue, int n);
int NameToType(char* name);
double sidechain_lower_limit(char *res,int chi);
double sidechain_upper_limit(char *res,int chi);
void print_params(const char *file);


void init_MOFA(int n);
void MO_FA_GMJ();
void copySolution(sol* A,sol* B);
void Normalizar_sol(int i);
void Maximos_iniciales();
void Normalizar(int objetivo);
int comparate_min(sol* a, sol *b);

void inicializacion_MOABC();
void CopyStruct(sol *dest, sol *src);

double ran();

extern int seed, depth, genes, archive, objectives, iterations, Tmax, minmax, clones; // command line parameters
extern int arclength; // current length of the archive
extern int num_evaluations;
extern int flag;
extern float sigma;

extern float gl_offset[MAX_OBJ];   // the range, offset etc of the grid
extern float gl_range[MAX_OBJ];
extern float gl_largest[MAX_OBJ];
extern int grid_pop[MAX_LOC];   // the array holding the population residing in each grid location


extern char problem[30];
extern char protein[200];
extern char runs_state_file[200];
extern char dir[200];
extern char params[50];

extern double best_energy_found;


