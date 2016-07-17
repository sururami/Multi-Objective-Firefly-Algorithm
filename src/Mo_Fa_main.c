// MO_FA_GMJ skeleton program code for the Protein Structure Prediction problem (PSP)
/*
  
*/
//

/*********************************************************************************************************
 * To run :                                                                                              *
 *  ./Mo_Fa_Main [Energy] [depth] [objectives] [residues] [archive] [iterations] [Tmax] [protein file]       *
 *           [dir] [SS]                                                                                  *
 *                                                                                                       *
 * where all parameters MUST be specified correctly (none are optional) following the instructions below:*
 *                                                                                                       *
 * [Energy] - at present there is only one Energy function included with this code: CHARMM (v. 22)       *
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
 * ./Mo_Fa_Main CHARMm22 4 2 34 1000 250000 250000 0 instances/1ZDD.seq 1ZDD 1                                 *
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
double max_obj0=0, max_obj1=0;


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

    printf("Inicio Programa\n");

    int result;
    int M2;
    double M1;
    int num_evaluations_previous;
    num_evaluations=0;
    while(num_evaluations <= 15000) {
     
    num_evaluations_previous = num_evaluations;
    for(int j=0; j< N_fireflies; j++){
        for(int k=0; k < N_fireflies; k++){
            if(j!=k){
                printf("-------------\n");
                result = dominate(&FireflyArray[k], &FireflyArray[j]);
                if(result==1){ // K Domina J ----------> Epsilon Dominating
                 //   printf("Existen dos luciernagas distintas %d y %d \n", j,k );

                    //M1 = Distancia Euclídea de las dos soluciones dividido entre la Raiz de 2
                    M1=sqrt(pow(FireflyArray[j].objNorm[0] - FireflyArray[k].objNorm[0],2)+pow(FireflyArray[j].objNorm[1] - FireflyArray[k].objNorm[1],2))/sqrt(2);
                    M1=exp(-2.0*M1);
                    //M2 = Modificador de M1 para realizarlo un 25 % de las veces. Se le suma 1 para hacerlo al menos 1 vez
                    M2= 1 + (int)(((double)(genes)/4) * M1);
                    //M1*=100;
         //           printf("Los parametros de mutacion estan activos \n" );

                    cl[0] = FireflyArray[k];   // Copia la Solucion en el primer Clon
                    cl[1] = FireflyArray[k];   // Copia la Solucion en el segundo Clon
    
                    global_mutation(&cl[0], M1); // Mutacion Global sobre el primer Clon
                    local_mutation(&cl[1], M2);  // Mutacion Local sobre el segundo Clon
                    
                    evaluate(&cl[0], problem); /*Evaluamos la primera mutacion y aumentamos el numero de evaluaciones*/
                    evaluate(&cl[1], problem); /*Evaluamos la segunda mutacion y aumentamos el numero de evaluaciones*/             

                    result = dominate(&cl[1],&cl[0]);
                    
                    printf("Comparacion\n");
                    if (result == 1){ 
                        cl[0]=cl[1];
                    }
                    
                    result = dominate(&cl[0],&FireflyArray[j]);
                    if(result==1){ //FB dominate FR
                        
                        
                        double max_local0=fabs(FireflyArray[j].obj[0]);
                        double max_local1=fabs(FireflyArray[j].obj[1]);
                        
                        FireflyArray[j]=cl[0];

                        if((max_local0 == max_obj0 && max_obj0 != fabs(cl[0].obj[0])) ||
                           (max_local1 == max_obj1 && max_obj1 != fabs(cl[0].obj[1]))){
                            Maximos_iniciales();
                            if (max_local0 == max_obj0 && max_obj0 != fabs(cl[0].obj[0])) Normalizar(0);
                            if (max_local1 == max_obj1 && max_obj1 != fabs(cl[0].obj[1])) Normalizar(1);
                        }
                        else{
                            Normalizar_sol(j);
                        }  
                    }
                }
                printf("num_evaluations realizadas  %d \n", num_evaluations );
            }       
        }      
    }    

    // Verificar estancamiento
    if (num_evaluations == num_evaluations_previous)
    {
     for(int i=0; i< N_fireflies; i++)
	print_sol(&FireflyArray[i],i,"Sol_Fin.txt");  
 
     init_MOFA(N_fireflies);
    }
 
  }
}

void Maximos_iniciales()
{

  max_obj0 = fabs(FireflyArray[0].obj[0]);
  max_obj1 = fabs(FireflyArray[0].obj[1]);

  for(int i=1; i < N_fireflies; i++){
   if (fabs(FireflyArray[i].obj[0]) > max_obj0) max_obj0=fabs(FireflyArray[i].obj[0]);
   if (fabs(FireflyArray[i].obj[1]) > max_obj1) max_obj1=fabs(FireflyArray[i].obj[1]);
  }
}
void Normalizar(int objetivo)
{

  double max_obj;
  
  if (objetivo == 0) max_obj=max_obj0;
  else max_obj=max_obj1;

  for (int i=0; i < N_fireflies; i++){
   FireflyArray[i].objNorm[objetivo]=FireflyArray[i].obj[objetivo]/max_obj;
   FireflyArray[i].objNorm[objetivo]=(FireflyArray[i].objNorm[objetivo]+1)/2;
  }
}
void Normalizar_sol(int i)
{

   if (fabs(FireflyArray[i].obj[0]) > max_obj0) { 
      max_obj0=fabs(FireflyArray[i].obj[0]);
      Normalizar(0);
   }
   else FireflyArray[i].objNorm[0]=(((FireflyArray[i].obj[0])/max_obj0)+1)/2;
   if (fabs(FireflyArray[i].obj[1]) > max_obj1) {
      max_obj1=fabs(FireflyArray[i].obj[1]);
      Normalizar(1);
   }
   else FireflyArray[i].objNorm[1]=(((FireflyArray[i].obj[1])/max_obj1)+1)/2;
}
int dominate(sol* s1, sol* s2)
{
/*Devuelve 1 si s1 domina a s2
  Devuelve -1 si s2 domina a s1
  Devuelve 0 en caso de no dominancia mutua
*/

if( (s1->obj[0] <= s2->obj[0]) && (s1->obj[1] <= s2->obj[1]) && ((s1->obj[0] < s2->obj[0]) || (s1->obj[1] < s2->obj[1] ) ) )
	return 1;

if( (s2->obj[0] <= s1->obj[0]) && (s2->obj[1] <= s1->obj[1]) && ((s2->obj[0] < s1->obj[0]) || (s2->obj[1] < s1->obj[1] ) ) )
	return -1;
return 0;
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


	 
  if (!strcmp(problem, "CHARMm22"))
	strcpy(params,"bin/params/charmm22.prm");
  

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
            Normalizar_sol(k);
          //  printf("------>>>>num_firefly[%d] --> energy: %f\n",k, FireflyArray[k].energy);
            }
            
    Normalizar(0);
    Normalizar(1);
    }
printf("Fin de Init\n");
return flag;
}

void copySolution(sol* A, sol *B)
{

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
void print_sol(sol* s,int p,int i,const char* file)
{    
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
void print_sol(sol* s,int p,const char* file)
{
    
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
    fprintf(fd,"%f,%f,%f\n", s->energy,s->obj[0],s->obj[1]);
  
    fclose(fd);   
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
  N_fireflies = atoi(argv[6]);
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

  cout << "N_Fireflies: "<< N_fireflies<<endl; 
  num_evaluations = 0;
  int flag=init_MOFA(N_fireflies);
  
  for(int i=0; i< N_fireflies; i++)
	print_sol(&FireflyArray[i],i,"Sol_Inicial.txt");  

  //printf("\n\nEstado de la Inicialiacion: %d [0 == ok --- 1 == Fallo]\n\n", flag);
  
  if(flag==0){
    MO_FA_GMJ();
   for(int i=0; i< N_fireflies; i++)
	print_sol(&FireflyArray[i],i,"Sol_Fin.txt");  
  }
  free(curr);
  free(FireflyArray);
  free(cl);
  free(m);
  free(arc); 
  free(app); 
  
}
