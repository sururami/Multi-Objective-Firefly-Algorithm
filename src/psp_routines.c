#include "global.h"

// Return the number of sidechain torsion angles required to fix
// the 3D position of the atoms for each residue.
int get_num_sidechain_angles(char *res)
{
	if( 	strcmp(res,"GLY") == 0 ||
		strcmp(res,"ALA") == 0 ||
		strcmp(res,"PRO") == 0 )
		return 0;
	else if(strcmp(res,"SER") == 0 ||
		strcmp(res,"CYS") == 0 ||
		strcmp(res,"THR") == 0 ||
		strcmp(res,"VAL") == 0 )
		return 1;
	else if(strcmp(res,"ILE") == 0 ||
		strcmp(res,"LEU") == 0 ||
		strcmp(res,"ASP") == 0 ||
		strcmp(res,"ASN") == 0 ||
		strcmp(res,"HIS") == 0 ||
		strcmp(res,"PHE") == 0 ||
		strcmp(res,"TYR") == 0 ||
		strcmp(res,"TRP") == 0 )
		return 2;
	else if(strcmp(res,"MET") == 0 ||
		strcmp(res,"GLU") == 0 ||
		strcmp(res,"GLN") == 0 )
		return 3;
	else if(strcmp(res,"LYS") == 0 ||
		strcmp(res,"ARG") == 0 )
		return 4;
	else 
		printf("Error: residue name unknow in 'get_num_sidechain_angles'\n");
}

int PredictedToNumber(char p)
{
  switch(p){
  case 'H': 	return 0;
    break;
  case 'E': 	return 1;
    break;
  case 'a': 	return 2;
    break;
  case 'b': 	return 3;
    break;
  case 'e': 	return 4;
    break;
  case 'l': 	return 5;
    break;
  case 't': 	return 6;
    break;
  case 'U': 	return 7;
    break;
  }
}

int PredictedToNumber_native(char p)
{
  switch(p){
  case 'H': 	return 0;
    break;
  case 'B': 	return 1;
    break;
  case 'E': 	return 2;
    break;
  case 'G': 	return 3;
    break;
  case 'I': 	return 4;
    break;
  case 'T': 	return 5;
    break;
  case 'S': 	return 6;
    break;
  case 'U': 	return 7;
    break;
  }
}

int NameToType(char* name)
{
  if(strcmp(name,"ALA")==0) return 0;
  else if(strcmp(name,"ARG")==0) return 1;
  else if(strcmp(name,"ASN")==0) return 2;
  else if(strcmp(name,"ASP")==0) return 3;	
  else if(strcmp(name,"CYS")==0) return 4;	
  else if(strcmp(name,"GLN")==0) return 5;
  else if(strcmp(name,"GLU")==0) return 6;
  else if(strcmp(name,"GLY")==0) return 7;		
  else if(strcmp(name,"HIS")==0) return 8;
  else if(strcmp(name,"ILE")==0) return 9;
  else if(strcmp(name,"LEU")==0) return 10;	
  else if(strcmp(name,"LYS")==0) return 11;	
  else if(strcmp(name,"MET")==0) return 12;
  else if(strcmp(name,"PHE")==0) return 13;
  else if(strcmp(name,"PRO")==0) return 14;
  else if(strcmp(name,"SER")==0) return 15;
  else if(strcmp(name,"THR")==0) return 16;
  else if(strcmp(name,"TRP")==0) return 17;
  else if(strcmp(name,"TYR")==0) return 18;
  else if(strcmp(name,"VAL")==0) return 19;
  else {
    printf("Error:Unknow residue in 'NameToType()'\n");
    return -1;
  }
}

// phi torsion angles are {\em \red bounded} in the constraint regions
void phi_constrain(float *angle, int p, char *name)
{	

  if(*angle < phi[p][0])
    *angle = phi[p][1] - fabs(*angle - phi[p][0]);// circular constraint
  //*angle = phi[p][0];
  else if(*angle > phi[p][1])
    *angle = phi[p][0] + fabs(*angle - phi[p][1]);// circular constraint
  //*angle = phi[p][1]; 
}

// psi torsion angles are {\em \red bounded} in the constraint regions
void psi_constrain(float *angle, int p, char *name)
{	
  if(*angle < psi[p][0])
    *angle = psi[p][1] - fabs(*angle - psi[p][0]);// circular constraint
  //*angle = psi[p][0];
  else if(*angle > psi[p][1])
    *angle = psi[p][0] + fabs(*angle - psi[p][1]);// circular constraint
  //*angle = psi[p][1];
}

// chi torsion angles are {\em \red bounded} in the constraint regions
void chi_constrain(float *angle, char* residue, int n)
{	
  double low, upper;
  
  low   = sidechain_lower_limit(residue,n);
  upper = sidechain_upper_limit(residue,n);

  if(*angle < low)
    *angle = upper - fabs(*angle - low);// circular constraint
  //*angle = low;
  else if(*angle > upper)
    *angle = low + fabs(*angle - upper);// circular constraint
  //*angle = upper;
}

// Randomly reselects all the angles of the residue 
// 'res' from the constraint regions  
void randConstAngles(res *r)
{	
  int i,p;
  
  p = r->predicted;
 
  r->angles[0] = RandomInt(phi[p][0],phi[p][1]) + RandomUniform();
  r->angles[1] = RandomInt(psi[p][0],psi[p][1]) + RandomUniform();
  
  for(i=2;i<r->num_angles;i++)
    {
      r->angles[i] = getSidechainAngles(r->name,i-1);
    }
}

float getSidechainAngles(char *name,int n)
{
  return RandomInt((int)sidechain_lower_limit(name,n),(int)sidechain_upper_limit(name,n)) + RandomUniform();
}

double sidechain_lower_limit(char *res,int chi)
{
  if(strcmp(res,"ARG") == 0)
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -167.0;
	break;
      case 3: 	return -65.0;
	break;
      case 4: 	return -175.0;
	break;
      }
    }
  else if(strcmp(res,"LYS") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -68.0;
	break;
      case 3: 	return -68.0;
	break;
      case 4: 	return -65.0;
	break;
      }
    }
  else if(strcmp(res,"MET") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -65.0;
	break;
      case 3: 	return -75.0;
	break;
      }
    }
  else if(strcmp(res,"GLU") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -80.0;
	break;
      case 3: 	return -60.0;
	break;
      }
    }
  else if(strcmp(res,"GLN") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -75.0;
	break;
      case 3: 	return -100.0;
	break;
      }
    }
  else if(strcmp(res,"ASP") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -60.0;
	break;
      }
    }
  else if(strcmp(res,"ASN") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -80.0;
	break;
      }
    }
  else if(strcmp(res,"ILE") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -60.0;
	break;
      }
    }
  else if(strcmp(res,"LEU") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return 65.0;
	break;
      }
    }
  else if(strcmp(res,"HIS") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -165.0;
	break;
      }
    }
  else if(strcmp(res,"TRP") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -105.0;
	break;
      }
    }
  else if(strcmp(res,"TYR") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -85.0;
	break;
      }
    }
  else if(strcmp(res,"PHE") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      case 2: 	return -85.0;
	break;
      }
    }
  else if(strcmp(res,"PRO") == 0 )
    {
      switch(chi){
      case 1:   return -30.0;
	break;
      }
    }
  else if(strcmp(res,"THR") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      }
    }
  else if(strcmp(res,"VAL") == 0 )
    {
      switch(chi){
      case 1:   return -60.0;
	break;
      }
    }
  else if(strcmp(res,"SER") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      }
    }
  else if(strcmp(res,"CYS") == 0 )
    {
      switch(chi){
      case 1:   return -177.0;
	break;
      }
    }
  else printf("Error: residue name unknow\n");
}

double sidechain_upper_limit(char *res,int chi)
{
  if(strcmp(res,"ARG") == 0)
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 180.0;
	break;
      case 3: 	return 180.0;
	break;
      case 4: 	return 180.0;
	break;
      }
    }
  else if(strcmp(res,"LYS") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 180.0;
	break;
      case 3: 	return 180.0;
	break;
      case 4: 	return 180.0;
	break;
      }
    }
  else if(strcmp(res,"MET") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 180.0;
	break;
      case 3: 	return 180.0;
	break;
      }
    }
  else if(strcmp(res,"GLU") == 0 )
    {
      switch(chi){
      case 1:   return 70.0;
	break;
      case 2: 	return 180.0;
	break;
      case 3: 	return 60.0;
	break;
      }
    }
  else if(strcmp(res,"GLN") == 0 )
    {
      switch(chi){
      case 1:   return 70.0;
	break;
      case 2: 	return 180.0;
	break;
      case 3: 	return 100.0;
	break;
      }
    }
  else if(strcmp(res,"ASP") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 65.0;
	break;
      }
    }
  else if(strcmp(res,"ASN") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 120.0;
	break;
      }
    }
  else if(strcmp(res,"ILE") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 170.0;
	break;
      }
    }
  else if(strcmp(res,"LEU") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 175.0;
	break;
      }
    }
  else if(strcmp(res,"HIS") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 165.0;
	break;
      }
    }
  else if(strcmp(res,"TRP") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 95.0;
	break;
      }
    }
  else if(strcmp(res,"TYR") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 90.0;
	break;
      }
    }
  else if(strcmp(res,"PHE") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      case 2: 	return 90.0;
	break;
      }
    }
  else if(strcmp(res,"PRO") == 0 )
    {
      switch(chi){
      case 1:   return 30.0;
	break;
      }
    }
  else if(strcmp(res,"THR") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      }
    }
  else if(strcmp(res,"VAL") == 0 )
    {
      switch(chi){
      case 1:   return 175.0;
	break;
      }
    }
  else if(strcmp(res,"SER") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      }
    }
  else if(strcmp(res,"CYS") == 0 )
    {
      switch(chi){
      case 1:   return 62.0;
	break;
      }
    }
  else printf("Error: residue name unknow\n");
}
