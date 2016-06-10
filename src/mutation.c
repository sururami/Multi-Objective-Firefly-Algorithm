#include "global.h"

// Local mutation: it perturbs some torsion angles 
// (phi, psi, chi) of a randomly chosen residue using 
// the standard Gaussian distribution.
void local_mutation(sol *s, int n_m)
{
  int i,j,k,n,z,p;

  for(k=0;k<n_m;k++)
  {
    i = RandomInt(0,genes-1);
    n = RandomInt(0,s->chrom[i].num_angles-1);
    
    for(z=0;z<n;z++)
      {
	p = RandomInt(0,s->chrom[i].num_angles-1);

	s->chrom[i].angles[p] += RandomGaussian(0.0,sigma);
	if(p==0)
	  phi_constrain(&s->chrom[i].angles[p],s->chrom[i].predicted,s->chrom[i].name);
	else if(p==1)
	  psi_constrain(&s->chrom[i].angles[p],s->chrom[i].predicted,s->chrom[i].name);
	else if(p>1)
	  chi_constrain(&s->chrom[i].angles[p],s->chrom[i].name,p-1);
	else
	  printf("Error:selected inesistent angle");
      }
  }
  s->to_evaluate = 1;
}	

// Global mutation: all the values of the backbone and sidechain torsion angles 
// of a randomly chosen residue are re-selected from 
// their corresponding constrained regions.
void global_mutation(sol *s, double pm) 
{
  int r,k;
  
  if(COIN(pm))
    {
      r = RandomInt(0,genes-1);
      randConstAngles(&s->chrom[r]);
      s->to_evaluate = 1;
    }    
  else {
    s->to_evaluate = 0;
  }
}
