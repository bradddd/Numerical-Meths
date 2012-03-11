void interp(double ***uf, double ***uc, int nf)
/*
  Coarse-to-fine prolongation by bilinear interpolation. nf is the fine-grid dimension. The coarsegrid
  solution is input as uc[1..nc][1..nc], where nc = nf/2 + 1. The fine-grid solution is
  returned in uf[1..nf][1..nf].
*/
{
  int ic,iif,jc,jf,nc,kc,kf;
  nc=nf/2+1;


  /* Do elements that are copies.*/
	for (kc=1,kf=1;kc <= nf;kc++,kf+=2)  
  		for (jc=1,jf=1;jc<=nc;jc++,jf+=2) 
    		for (ic=1,iif=1;ic<=nc;ic++,iif+=2)
    	 	 	uf[iif][jf][kf]=uc[ic][jc][kc];
    	 	 	
    	 	 	
  /* Do odd-numbered columns, interpolating vertically.*/
  
  for (jf=1;jf<=nf;jf+=2) 
    for (iif=2;iif<nf;iif+=2)
    	for (kf=1;kf <= nf;kf+=2) 
      uf[iif][jf][kf]=0.5*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]);     
      
  /*Do even-numbered columns, interpolating horizontally.*/
  for (jf=2;jf<nf;jf+=2) 
    for (iif=1;iif <= nf;iif++)
    	for (kf=1;kf <= nf;kf+=2)  
      uf[iif][jf][kf]=0.5*(uf[iif][jf+1][kf]+uf[iif][jf-1][kf]);
      
  /*Do even-numbered columns, interpolating horizontally.*/
  for (jf=1;jf<nf;jf++) 
    for (iif=1;iif<= nf;iif++)
    	for (kf=2;kf< nf;kf+=2)  
      	uf[iif][jf][kf]=0.5*(uf[iif][jf][kf+1]+uf[iif][jf][kf-1]);
}
