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

int N_fireflies;
double max_obj1=0, max_obj2=0;


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

void MO_FA_GMJ()
{
printf("Hola desde MO_FA\n");
printf("Hola desde MO_FA\n");
printf("Hola desde MO_FA\n");
printf("Hola desde MO_FA\n");
printf("Hola desde MO_FA\n");
printf("Hola desde MO_FA\n");
    int result;
  
    int M2;
    double M1;
  
    int flag=0;
    int num_evaluations=0;
    for(int j=0; j< N_fireflies && num_evaluations <= 1000 && flag==0; j++){
        for(int k=0; k < N_fireflies && num_evaluations <= 1000 && flag==0; k++){
            if(j!=k){		

                result = comparate_min(&FireflyArray[k], &FireflyArray[j]);
                if(result==0){ // K Domina J ----------> Epsilon Dominating
                    //M1 = Distancia Euclidea de las dos soluciones dividido entre la Raiz de 2
                    M1=sqrt(pow(FireflyArray[j].objNorm[0] - FireflyArray[k].objNorm[0],2)+pow(FireflyArray[j].objNorm[1] - FireflyArray[k].objNorm[1],2))/sqrt(2);
                    //M2 = Modificador de M1 para realizarlo un 25 % de las veces. Se le suma 1 para hacerlo al menos 1 vez
                    M2= 1 + (int)(((double)(genes)/4) * M1);
                    M1*=100;

                    cl[0] = FireflyArray[j];   // Copia la Solucion en el primer Clon
                    cl[1] = FireflyArray[j];   // Copia la Solucion en el segundo Clon
    
                    global_mutation(&cl[0], M1); // Mutacion Global sobre el primer Clon
                    local_mutation(&cl[1], M2);  // Mutacion Local sobre el segundo Clon
            
                    evaluate(&cl[0], problem);
                    evaluate(&cl[1], problem);              

	            num_evaluations+=2;

                    result = comparate_min(&cl[0],&cl[1]);
                    if (result == 0) 
                        cl[0]=cl[1];
                    result = comparate_min(&cl[0],&FireflyArray[j]);
                    if(result==1){ //FB dominate FR
                            printf("Se ha producido una mutacion\n" );
                        FireflyArray[j]=cl[0];
                        Normalizar_sol(j);
                    }
                }
                //printf("\n\n-------- %d --------\n\n", num_evaluations);
                print_sol(&FireflyArray[j],j,"Energia.energia");
            }       
        }      
    }    
}

void Maximos_iniciales(){

  double a=FireflyArray[0].obj[0];

  max_obj1 = fabs(a);
  max_obj2 = fabs(FireflyArray[0].obj[1]);

  for(int i=1; i < N_fireflies; i++){
   if (fabs(FireflyArray[i].obj[0]) > max_obj1) max_obj1=fabs(FireflyArray[i].obj[0]);
   if (fabs(FireflyArray[i].obj[1]) > max_obj2) max_obj2=fabs(FireflyArray[i].obj[1]);
  }
}


void Normalizar(int objetivo){

  for (int i=0; i < N_fireflies; i++){
   FireflyArray[i].objNorm[objetivo]=(-1*(FireflyArray[i]).obj[objetivo])/max_obj1;
   FireflyArray[i].objNorm[objetivo]=((FireflyArray[i]).objNorm[objetivo]+1)/2;
  }
}


void Normalizar_sol(int i){

   if (fabs(FireflyArray[i].obj[0]) > max_obj1) { 
      max_obj1=fabs(FireflyArray[i].obj[0]);
      Normalizar(0);
   }
   else FireflyArray[i].objNorm[0]=(-1*FireflyArray[i].obj[0])/max_obj1;
   if (fabs(FireflyArray[i].obj[1]) > max_obj2) {
      max_obj2=fabs(FireflyArray[i].obj[1]);
      Normalizar(1);
   }
   else FireflyArray[i].objNorm[1]=(-1*FireflyArray[i].obj[1])/max_obj1;
}

int comparate_min(sol* s1, sol* s2){
/*Devuelve 0 si s1 domina a s2
  Devuelve 1 si s2 domina a s1
  Devuelve -1 en caso de igualdad
*/

if( (s1->obj[0] <= s2->obj[0]) && (s1->obj[1] <= s2->obj[1]) && ((s1->obj[0] < s2->obj[0]) || (s1->obj[1] < s2->obj[1] ) ) )
	return 0;

if( (s2->obj[0] <= s1->obj[0]) && (s2->obj[1] <= s1->obj[1]) && ((s2->obj[0] < s1->obj[0]) || (s2->obj[1] < s1->obj[1] ) ) )
	return 1;

return -1;


}

int init_MOFA(int N_Fireflies)
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
  

  // Reserva de memoria de los array y datos usados


  FireflyArray = (sol*) malloc(N_Fireflies *sizeof(sol));

  curr   = (sol *)malloc(sizeof(sol));
  cl  = (sol *)malloc(MAX_CLONES*sizeof(sol));
  m   = (sol *)malloc(sizeof(sol));
  arc = (sol *)malloc(MAX_ARC*sizeof(sol));
  app = (sol *)malloc(MAX_ARC*sizeof(sol));

  if((!curr)||(!m)||(!arc) || (!FireflyArray))
    {
      printf("Out of memory. Aborting.\n");
      exit(-1);
    }

  if(native_ss == 1)
    {
      //Restricciones Estructuras Secundarias
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
      //Restricciones Estructuras Supersecundarias
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
      printf("fopen filed to open %s\n", protein);
      exit(1);
    }
 int flag=0;
    for (i = 0; i < genes && flag == 0; i++)
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
    if(flag!=1){
        // Inicializa FireflyArray[k] seleccionando aleatoriamente los angulos de torsion y de la cadena principal 
        // en las regiones asociadas.

        for(int k=0; k < N_Fireflies; k++){ 
            FireflyArray[k]=*curr;
             printf(" %d\n",k);
            for (int j = 0; j < genes; j++)
            {
                res* r = &(FireflyArray[k].chrom[j]);
                randConstAngles(r);
            }
            evaluate(&FireflyArray[k], problem);
          //  printf("------>>>>num_firefly[%d] --> energy: %f\n",k, FireflyArray[k].energy);
            }
    }
printf("Fin de Init\n");
return flag;
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

  srandom(time(NULL)); // seed para la funcion ramdom() 
  RandSeed=random();   // seed para la funcion ran()
  srand(r=random());

  RandomInitialise(Randint(0,31328),Randint(0,30081));
   
  sigma = 1.0;

  strcat(command,dir);
  system(command);

  print_params("parameters.dat");

    signal(SIGINT,backup);// handler for program termination by keyboard (SIGINT signal)

  //double time1=omp_get_wtime();
  N_fireflies= 20;
  printf("\n\n\n Hola desde Antes de Init\n");
  int flag=init_MOFA(N_fireflies);
//  if(flag==0)
//    MO_FA_GMJ();
  
  free(curr);
  free(FireflyArray);
  free(cl);
  free(m);
  free(arc); 
  free(app); 
  
}

void print_sol(sol* s,int p,int i,const char* file){
    
    FILE *fd;
    char params_file[200];

    strcpy(params_file,dir);
    strcat(params_file,"/");
    strcat(params_file,file);
    fd = fopen(params_file,"aw");
 
    if(fd == NULL)
    {
        printf("fopen fallo al abrir %s\n",params_file);
        exit(-1);
    }
    fprintf(fd,"/******************************************\n");
    fprintf(fd,"Generacion %d  --- Datos de solucion %d\n",i,p );
    fprintf(fd,"Energia: %f      Bond: %lf      Non_Bond: %lf\n", s->energy,s->obj[0],s->obj[1]);
    fprintf(fd,"******************************************/\n\n");
  
    fclose(fd);
    
}

void print_sol(sol* s,int p,const char* file){
    
    FILE *fd;
    char params_file[200];

    strcpy(params_file,dir);
    strcat(params_file,"/");
    strcat(params_file,file);
    fd = fopen(params_file,"aw");
 
    if(fd == NULL)
    {
        printf("fopen fallo al abrir %s\n",params_file);
        exit(-1);
    }
    fprintf(fd,"%f,%lf,%lf\n", s->energy,s->obj[0],s->obj[1]);
  
    fclose(fd);   
}

int kbhit()
{
  struct termios oldt, newt;
  int ch;
  int oldf;
 
  tcgetattr(STDIN_FILENO, &oldt);
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);
 
  ch = getchar();
 
  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  fcntl(STDIN_FILENO, F_SETFL, oldf);
 
  if(ch != EOF)
  {
    ungetc(ch, stdin);
    return 1;
  }
 
  return 0;
}


