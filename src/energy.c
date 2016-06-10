#include "global.h"


//CHARMM (v.27) energy function 
void CHARMm27(sol *s)
{
  int i;
  double bond; 
  double non_bond; 
  
  int fdesc;
  int fdesc2;
  
  double bond_stretching;
  double angle_bending;
  double urey_bradley;
  double improper_dihedral;
  double torsion_angle;
  double van_der_waals;
  double charge_charge;
  double tot;
  
  char buffer[160];
  char c[100];
  char boo[20];
  int exitcode;
  
  FILE *fd;

  // PROTEIN input file creation
  protein_file("prot.dat", s, protein, "prot");

  // .xyz file creation 
  char *args1[] = { "bin/protein", NULL };
  if (fork() == 0) {
    /* child process */
    close(0); //close stdin
    //close(1); //close stdout
    fdesc = open("prot.dat",O_RDONLY);
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

  // Energy calculation
  char *args2[] = { "bin/analyze" , "prot.xyz", params, "E", NULL };
  if (fork() == 0) {
    /* child process */
    close(1);      // close stdout (assumed open)
    fdesc = open("bin/energy.dat",O_CREAT|O_RDWR,0777);
    if(fdesc == -1)
      {
  	printf("open filed to open energy.dat\n");
  	exit(-1);
      }
    dup(fdesc);    // dups arg to min free descriptor
    execv(args2[0],args2);
    printf("cannot exec: %s\n",args2[0]);
    exit(-1);
  }
  else {
    /* parent process */
    if (wait(&exitcode) != -1) {
      close(fdesc);
      //printf("\ndone\n");
    }
  }
     
  // Bond and non-bond energy values extraction
  fd = fopen("bin/energy.dat", "r");
  if(fd == NULL)
    {
      printf("fopen filed to open energy.dat\n");
      exit(-1);
    }

  for(i=0;i<19;i++)
    fgets(buffer,159,fd);

  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %s",c, c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    bond_stretching = atof(boo);
    //printf("Bond Stretching: %lf\n", bond_stretching);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %s",c, c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    angle_bending = atof(boo);
    //printf("Angle Bending: %lf\n", angle_bending);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s",c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    urey_bradley = atof(boo);
    //printf("Urey-Bradley: %lf\n", urey_bradley);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %s",c, c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    improper_dihedral = atof(boo);
    //printf("Improper Dihedral: %lf\n", improper_dihedral);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %s",c, c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    torsion_angle = atof(boo);
    //printf("Torsional Angle: %lf\n", torsion_angle);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %s %s",c, c, c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    van_der_waals = atof(boo);
    //printf("Van der Waals: %lf\n", van_der_waals);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s",c, boo);
    for(i=0;i<20;i++)
      if(boo[i]=='D') boo[i]='E';
    charge_charge = atof(boo);
    //printf("Charge-Charge: %lf\n", charge_charge);
  }
  
  /*
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %lf",c, c, &bond_stretching);
    //printf("Bond Stretching: %lf\n", bond_stretching);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %lf",c, c, &angle_bending);
    //printf("Angle Bending: %lf\n", angle_bending);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s  %lf",c, &urey_bradley);
    //printf("Urey-Bradley: %lf\n", urey_bradley);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %lf",c, c, &improper_dihedral);
    //printf("Improper Dihedral: %lf\n", improper_dihedral);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %lf",c, c, &torsion_angle);
    //printf("Torsional Angle: %lf\n", torsion_angle);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %s %s %lf",c, c, c, &van_der_waals);
    //printf("Van der Waals: %lf\n", van_der_waals);
  }
  if (fgets(buffer,159,fd) != NULL) {
    sscanf(buffer,"%s %lf",c, &charge_charge);
    //printf("Charge-Charge: %lf\n", charge_charge);
    }
  */

  fclose(fd);
  
  bond = bond_stretching + angle_bending + urey_bradley + improper_dihedral + torsion_angle;
  non_bond = van_der_waals + charge_charge;
  
  
  //printf("Non-Bond: %f\n", non_bond);
  
  s->obj[0] = bond;
  s->obj[1] = non_bond;
  s->energy = bond + non_bond;
   
  //Elimination of superfluous files 
  char *args4[] = { "/bin/rm" ,"prot.dat", "prot.seq", "prot.xyz", "prot.int", NULL };
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

//Routine for creating the input file for the external routine PROTEIN
void protein_file(const char *file, sol *s, const char *title, const char *filexyz)
{
  int i,j;
  FILE *fdesc;
  
  fdesc = fopen(file,"w");
  if(fdesc == NULL)
    {
      printf("fopen filed to open %s\n",file);
      exit(-1);
    }
  
  fprintf(fdesc,"%s\n", filexyz);
  fprintf(fdesc,"%s\n", title);
  fprintf(fdesc,"%s\n", params);
  
	
  for(i=0;i<genes;i++)
    {
      fprintf(fdesc,"%s %.3f %.3f 180.0 ",s->chrom[i].name, s->chrom[i].angles[0], s->chrom[i].angles[1]);
      for(j=2;j<s->chrom[i].num_angles;j++)
	{
	  fprintf(fdesc,"%.3f ", s->chrom[i].angles[j]);
	}
      fprintf(fdesc,"\n");
    }
  
  fprintf(fdesc,"\nn\n");
  fprintf(fdesc,"\nEnergy: %.4f\n",s->energy);
  fclose(fdesc);
}
