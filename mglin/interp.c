#include <math.h>
#include <stdio.h>

void interp(double ***uf, double ***uc, int nf, int cubic)
/*
  Coarse-to-fine prolongation by bilinear interpolation. nf is the fine-grid dimension. The coarsegrid
  solution is input as uc[1..nc][1..nc], where nc = nf/2 + 1. The fine-grid solution is
  returned in uf[1..nf][1..nf].
*/
{
	//int cubic = 0;	
	
	if(cubic){
		
		// basic set up
		// u[i][j][k] = 		
		
		//diag ratio
		double d = sqrt(3);
		
		int ic,iif,jc,jf,nc,kc,kf;
	  nc=nf/2+1;
	
	
	  /* Do elements that are copies.*/
		for (kc=1,kf=1;kc <= nc;kc++,kf+=2)  
	  		for (jc=1,jf=1;jc<=nc;jc++,jf+=2) 
	    		for (ic=1,iif=1;ic<=nc;ic++,iif+=2)
	    	 	 	uf[iif][jf][kf]=uc[ic][jc][kc];
	    	 	 	
	    	 	 	
	  /* Do odd-numbered columns, interpolating vertically.*/
	  //printf("odd col\n");
	  for (jf=1;jf<=nf;jf+=2) 
	    for (iif=2;iif<nf;iif+=2)
	    	for (kf=1;kf <= nf;kf+=2) {
	    		if (jf ==1 || jf ==nf ){
	    			//printf("in odd if\n");	
	      		uf[iif][jf][kf]=0.5*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]);
	      	}
	      	else {
	      		//printf("j %d i %d k %d \n",jf,iif,kf);
	      		uf[iif][jf][kf]=0.25*((2*d-1)/2*d)*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]) + (0.125/d)*(uf[iif+1][jf-1][kf]+uf[iif-1][jf+1][kf]+uf[iif+1][jf+1][kf]+uf[iif-1][jf-1][kf]);
	      	}  
	      }   
	  //printf("even col\n");
	  /*Do even-numbered columns, interpolating horizontally.*/
	  for (jf=2;jf<nf;jf+=2) 
	    for (iif=1;iif <= nf;iif++)
	    	for (kf=1;kf <= nf;kf+=2) {  
	    		if (iif ==1 || iif ==nf ){	
	    			//printf("in even if\n");
	      		uf[iif][jf][kf]=0.5*(uf[iif][jf+1][kf]+uf[iif][jf-1][kf]);
	      	}
	      	else {
	      		//printf("j %d i %d k %d \n",jf,iif,kf);
	      		uf[iif][jf][kf]=0.25*((2*d-1)/2*d)*(uf[iif][jf+1][kf]+uf[iif][jf-1][kf]) + (0.125/d)*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]+uf[iif+1][jf][kf]+uf[iif-1][jf][kf]);
	      	}
	   	}
	   	
	   //printf("odd sheets\n");   
	  /*Do even-numbered columns, interpolating horizontally.*/
	  for (jf=1;jf<nf;jf++) 
	    for (iif=1;iif<= nf;iif++)
	    	for (kf=2;kf< nf;kf+=2) { 
	    		//printf("j %d i %d k %d \n",jf,iif,kf);
	    		if (kf == 2 || kf == (nf-2) || jf == 1 || jf == nf || iif == 1|| iif == nf ){
	    			uf[iif][jf][kf]=0.25*((d-1)/d)*(uf[iif][jf][kf+1]+uf[iif][jf][kf-1]);
	    		}
	    		else {
		      	uf[iif][jf][kf]=0.25*((d-1)/d)*(uf[iif][jf][kf+1]+uf[iif][jf][kf-1]) + (0.125/(2*d))*(uf[iif+1][jf+1][kf+1]
		      																							+uf[iif+1][jf+1][kf-1]+uf[iif+1][jf-1][kf+1]+uf[iif-1][jf+1][kf+1]
		      																							+uf[iif-1][jf-1][kf+1]+uf[iif-1][jf+1][kf-1]+uf[iif+1][jf-1][kf-1]
	      																								+uf[iif-1][jf-1][kf-1]);
	      	}
	      }
		}
	else {
		
	  int ic,iif,jc,jf,nc,kc,kf;
	  nc=nf/2+1;
	
	
	  /* Do elements that are copies.*/
		for (kc=1,kf=1;kc <= nc;kc++,kf+=2)  
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
	      
	  /*Do blank sheets interpolating from sheet in front and behind.*/
	  for (jf=1;jf<nf;jf++) 
	    for (iif=1;iif<= nf;iif++)
	    	for (kf=2;kf< nf;kf+=2)  
	      	uf[iif][jf][kf]=0.5*(uf[iif][jf][kf+1]+uf[iif][jf][kf-1]);
      }
}
