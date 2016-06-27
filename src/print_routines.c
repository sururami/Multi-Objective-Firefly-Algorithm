#include "global.h"

// writes the simulation parameters on file
void print_params(const char *file)
{        
  FILE *fd;
  char params_file[200];

  strcpy(params_file,dir);
  strcat(params_file,"/");
  strcat(params_file,file);

  fd = fopen(params_file,"aw");
  if(fd == NULL)
    {
      printf("fopen filed to open %s\n",params_file);
      exit(-1);
    }

  fprintf(fd,"/******************************************\n");
  fprintf(fd,"   Parameters of last execution\n");
  fprintf(fd,"   of PSP instance '%s'\n",protein);
  fprintf(fd,"******************************************/\n\n");
  fprintf(fd,"Energy     = %s # Energy function\n",problem );
  fprintf(fd,"depth      = %d # number of recursive subdivisions of the objective space\n",depth);
  fprintf(fd,"objectives = %d # number of objectives to the problem\n",objectives);
  fprintf(fd,"residues   = %d # number of residues in the protein\n",genes);
  fprintf(fd,"archive    = %d # maximum number of solutions to be held in the nondominated solutions archive\n",archive);
  fprintf(fd,"iterations = %d # max number of iterations allowed\n",iterations);
  fprintf(fd,"Tmax       = %d # max number of fitness function evaluations allowed\n",Tmax);
  fprintf(fd,"protein    = %s # input protein  instance file\n",protein);
  fprintf(fd,"dir        = %s # output directory for simulation results\n",dir);
  fprintf(fd,"SS         = %d # Structure constraints (1=Secondary,0=Supersecondary)\n",native_ss);
   
  fclose(fd);
}

// writes objectives and total energy of the 
// solutions in the archive in a file 
void print_archive(const char *file)
{
  int i,j;
  FILE *fdesc;
  float min, max;
  
  fdesc = fopen(file,"w");
  if(fdesc == NULL)
    {
      printf("fopen filed to open %s\n",file);
      exit(-1);
    }
  
  fprintf(fdesc,"#%d-objective\tTotEnergy\n",objectives);

  for (i = 0; i < arclength; i++)
    {
      for (j = 0; j < objectives; j++)
	fprintf(fdesc,"%lf\t",arc[i].obj[j]);

      fprintf(fdesc,"%lf\t",arc[i].energy);
      fprintf(fdesc,"\n");
    }
  fclose(fdesc);
}

void print_eval(sol *s)
{
  int i;
  double sum_obj = 0.0;
  for (i = 0; i < objectives; i++)
  {
	printf("%lf ", s->obj[i]);
	sum_obj += s->obj[i];
  }
  printf("%lf ", sum_obj);
  printf("\n");
}

// creates a directory with the pdb files of the 
// solutions in the archive at time step n
void print_arc(int n)
{
  int i;
  char pdb_dir[400];
  char tmp[10];
  char command[301] = "mkdir "; //command for directory creation
  char arc_file[200];

  strcpy(pdb_dir,dir);
  strcat(pdb_dir,"/pdb_files");
  strcat(pdb_dir,gcvt(n,5,tmp));
  strcat(command,pdb_dir);
  system(command);

  for (i = 0; i < arclength; i++)
    print_genome(&arc[i], i+1, pdb_dir);

  strcpy(arc_file,pdb_dir);
  strcat(arc_file,"/archive.dat");
  print_archive(arc_file);
}

// creates the pdb file of the solution s 
void print_genome(sol *s, int n, const char* dir)
{
  int i;
  char file[300];
  char xyz[300];
  char inte[300];
  char seq[300];
  char num[5];
  int fdesc;
	
  int exitcode;

  strcpy(file,dir);
  strcat(file,"/sol");
  strcat(file,gcvt(n,5,num));
  strcat(file,".dat");
  
  strcpy(xyz,dir);
  strcat(xyz,"/sol");
  strcat(xyz,gcvt(n,5,num));
  strcat(xyz,".xyz");
   
  strcpy(seq,dir);
  strcat(seq,"/sol");
  strcat(seq,gcvt(n,5,num));
  strcat(seq,".seq");
   
  strcpy(inte,dir);
  strcat(inte,"/sol");
  strcat(inte,gcvt(n,5,num));
  strcat(inte,".int");  
  
  protein_file(file, s, protein, xyz);
  
  char *args1[] = { "bin/protein", NULL };
  if (fork() == 0) {
    /* child process */

    close(0); //close stdin
    //close(1); //close stdout
    fdesc = open(file, O_RDONLY);
    dup(fdesc);
    execv(args1[0],args1);
    printf("cannot exec: %s\n",args1[0]);
    exit(-1);
  }
  else {
    /* parent process */
    if (wait(&exitcode) != -1) {
      //printf("\ndone\n");
    }
  }
  
  char *args2[] = { "bin/xyzpdb" , xyz, params, NULL };
  if (fork() == 0) {
    /* child process */
    //close(1); //close stdout
    execv(args2[0],args2);
    printf("cannot exec: %s\n",args2[0]);
    exit(-1);
  }
  else {
    /* parent process */
    if (wait(&exitcode) != -1) {
      //printf("\ndone\n");
    }
  }
  
  char *args4[] = { "/bin/rm" , seq, xyz, inte, NULL };
  if (fork() == 0) {
    /* child process */
    execv(args4[0],args4);
    printf("cannot exec: %s\n",args4[0]);
    exit(-1);
  }
  else {
    /* parent process */
    if (wait(&exitcode) != -1) {
      //printf("\ndone\n");
    }
  }
}

// print statistics results of the current archive 
void print_statistics(const char* name)
{
  int w;
  char stat_file[300];
  double sum_bond = 0.0;
  double sum_non_bond = 0.0;
  FILE *fdesc;

  strcpy(stat_file,dir);
  strcat(stat_file,"/");
  strcat(stat_file,name);
  
  fdesc = fopen(stat_file,"aw");
  if(fdesc == NULL)
    {
      printf("fopen filed to open %s\n",stat_file);
      exit(-1);
    }
 
  sum_bond = arc[0].obj[0];
  sum_non_bond = arc[0].obj[1];    

  for (w = 1; w < arclength; w++)
    {
      sum_bond += arc[w].obj[0];
      sum_non_bond += arc[w].obj[1];
    }
  fprintf(fdesc,"%lf\t%lf\t%lf\t%d\t%lf\n", 
	  sum_bond/arclength, sum_non_bond/arclength, (sum_bond+sum_non_bond)/arclength, 
	  arclength, best_energy_found);
  
  fclose(fdesc);
}

// print statistics results of the current solution
void print_curr_sol(const char* name)
{
  int w;
  char stat_file[300];
  FILE *fdesc;

  strcpy(stat_file,dir);
  strcat(stat_file,"/");
  strcat(stat_file,name);
  
  fdesc = fopen(stat_file,"aw");
  if(fdesc == NULL)
    {
      printf("fopen filed to open %s\n",stat_file);
      exit(-1);
    }

  fprintf(fdesc,"%lf\t",curr->energy);

  for(w=0;w<objectives;w++)
    fprintf(fdesc,"%lf\t",curr->obj[w]);

  fprintf(fdesc,"\n");
  
  fclose(fdesc);
}

// handler for program termination by keyboard (SIGINT signal) 
void backup(int signo)
{
  int i;
  int exitcode;
  char arch[300];
	
  printf("Signal handler started. Terminating ...\n");
	
  printf("\n#The Archive untill now is... \n");
  for (i = 0; i < arclength; i++)
    print_eval(&arc[i]);

  printf("\nSaving of the genetic material in the archive...\n");
  print_arc(iterations);

  char *args4[] = { "/bin/rm" , "prot.dat", "prot.seq", "prot.xyz", "prot.int", NULL };
  if (fork() == 0) {
    /* child process */
    execv(args4[0],args4);
    printf("cannot exec: %s\n",args4[0]);
    exit(-1);
  }
  else {
    /* parent process */
    if (wait(&exitcode) != -1) {
      //printf("\ndone\n");
    }
  }
    free(curr);
    free(FireflyArray);
    free(cl);
    free(m);
    free(arc); 
    free(app); 
    printf("Done.\n");
    exit(-1);
}
