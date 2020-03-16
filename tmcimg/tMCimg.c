#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#include "nrutil.h"
*/

#define pi 3.1415926353
#define C_VACUUM 2.9979e11
#define TRUE 1
#define FALSE 0


double ran( int *idum, int *ncall );	/* FUNCTION DECLARATION */
int idum;	/* SEED FOR RANDOM NUMBER GENERATOR - A LARGE NEGATIVE NUMBER IS REQUIRED */
int ncall;	/* NUMBER OF CALLS TO RANDOM NUMBER GENERATOR - USED
		   BY THE GENERATOR */

int main( int argc, char *argv[] )
{
  int i,j,k,ii,jj,kk;
  int N;			/* NUMBER OF PHOTONS RUN SO FAR */
  int NT;			/* TOTAL NUMBER OF PHOTONS TO RUN */
  int Ntissue;			/* NUMBER OF TISSUE TYPES DESCRIBED IN THE IMAGE FILE */
  int VoxN;                     /* INDEX TO THE CURRENT VOXEL CONTAINING THE PHOTON */

  double foo;    		/* TEMPORARY VARIABLES */
  float ffoo;

  char ***tissueType;                           /* STORE THE IMAGE FILE */
  int tissueIndex;
  int nxstep, nystep, nzstep;                   /* DIMENSIONS OF THE IMAGE FILE */
  double xstep, ystep, zstep;                   /* VOXEL DIMENSIONS */

  double tmus[10], tmua[10], tg[10], tn[10];    /* OPTICAL PROPERTIES OF THE DIFFERENT
						   TISSUE TYPES */

  double x,y,z;			/* CURRENT PHOTON POSITION */
  double xold,yold,zold;	/* LAST POSITION OF PHOTON */
  double xi, yi, zi;            /* INITIAL POSITION OF THE PHOTON */
  
  double gg, phi,theta,sphi,cphi,stheta,ctheta;	/* SCATTERING ANGLES */
  double c1,c2,c3;		/* DIRECTION COSINES */
  double c1o, c2o, c3o;		/* OLD DIRECTION COSINES */
  double cxi, cyi, czi;  	/* INITIAL DIRECTION COSINES */

  double *II, *pII;             /* FOR STORING THE 2-PT FLUENCE */

  int Ixmin, Ixmax, Iymin, Iymax, Izmin, Izmax;   /* MIN AND MAX X,Y,Z FOR STORING THE */
  int nIxstep, nIystep, nIzstep;                  /*   2-PT FLUENCE */

  double minT, maxT;            /* MIN AND MAX TIME FOR SAMPLING THE 2-PT FLUENCE */
  double stepT, stepL;          /* TIME STEP AND CORRESPONDING LENGTH STEP FOR SAMPLING 
                                   THE 2-PT FLUENCE */
  double Lresid, Ltot, Lmin, Lmax, Lnext, step;
  int nTstep, tindex;

  int nDets, detRad;            /* SPECIFY NUMBER OF DETECTORS, DETECTOR RADIUS */
  int **detLoc;                 /*    AND X,Y,Z LOCATIONS */

  double P2pt;			/* PHOTON WEIGHT */

  float lenTiss[10];            /* THE LENGTH SPENT IN EACH TISSUE TYPE BY THE CURRENT PHOTON */

  double rnm;			/* RANDOM NUMBER */

  FILE *fp;             	/* FILE POINTERS FOR SAVING THE DATA */
  char filenm[32];	/* FILE NAME FOR DATA FILE */
  char segFile[32];     /* FILE NAME FOR IMAGE FILE */

 
  


  /* GET THE COMMAND LINE ARGUMENTS */
  if( argc!=2) {
    printf( "usage: tMCimg input_file (.inp assumed)\n" );
    exit(1);
  }


  /*********************************************************
    OPEN AND READ THE INPUT FILE 
   *********************************************************/
  sprintf( filenm, "%s.inp", argv[1] );
  if( (fp = fopen( filenm, "r" ))==NULL ) {
    printf( "usage: tMCimg input_file (.inp assumed)\n" );
    printf( "input_file = %s does not exist.\n", filenm );
    exit(1);
  }

  /* READ THE INPUT FILE */
  fscanf( fp, "%d", &NT );    /* TOTAL NUMBER OF PHOTONS */
  fscanf( fp, "%d", &idum );  /* RANDOM NUMBER SEED */
  fscanf( fp, "%lf %lf %lf", &xi, &yi, &zi );         /* INITIAL POSITION OF PHOTON */
  fscanf( fp, "%lf %lf %lf", &cxi, &cyi, &czi );      /* INITIAL DIRECTION OF PHOTON */
  fscanf( fp, "%lf %lf %lf", &minT, &maxT, &stepT );  /* MIN, MAX, STEP TIME FOR RECORDING */
  nTstep = (int)ceil((maxT-minT)/stepT);
  fscanf( fp, "%s", segFile );                        /* FILE CONTAINING TISSUE STRUCTURE */

  if( idum>=0 ) {
    printf( "ERROR: The random number seed must be a large negative number, not %d.\n", idum );
    exit(1);
  }

  /* READ IMAGE DIMENSIONS */
  fscanf( fp, "%lf %d %d %d", &xstep, &nxstep, &Ixmin, &Ixmax );   
  fscanf( fp, "%lf %d %d %d", &ystep, &nystep, &Iymin, &Iymax );
  fscanf( fp, "%lf %d %d %d", &zstep, &nzstep, &Izmin, &Izmax );
  xi--; yi--; zi--;
  Ixmin--; Ixmax--; Iymin--; Iymax--; Izmin--; Izmax--;
  nIxstep = Ixmax-Ixmin+1;
  nIystep = Iymax-Iymin+1;
  nIzstep = Izmax-Izmin+1;
  if( xstep!=1. || ystep!=1. || zstep!=1. ) {
    printf( "This code requires xstep = ystep = zstep = 1 mm.\n" );
    return 0;
  }

  /* READ NUMBER OF TISSUE TYPES AND THEIR OPTICAL PROPERTIES */
  fscanf( fp, "%d", &Ntissue );
  tmus[0] = -999.; tmua[0] = -999.; tg[0] = -999.; tn[0] = -999.;
  for( i=1; i<=Ntissue; i++ ) {
    fscanf( fp, "%lf %lf %lf %lf", &tmus[i], &tg[i], &tmua[i], &tn[i] );
    if( tn[i]!=1.0 ) {
      printf( "ERROR: The code does not yet support n!=1.0\n" );
      return(0);
    }
    if( tmus[i]==0.0 ) {
      printf( "ERROR: The code does support mus = 0.0\n" );
      return(0);
    }
  }

  /* READ NUMBER OF DETECTORS, DETECTOR RADIUS, AND DETECTOR LOCATIONS */
  fscanf( fp, "%d %d", &nDets, &detRad );
  detLoc = (int **)malloc(nDets*sizeof(int*));
  for( i=0; i<nDets; i++ )
    detLoc[i] = (int *)malloc(3*sizeof(int));
  for( i=0; i<nDets; i++ ) {
    fscanf( fp, "%d %d %d", &detLoc[i][0], &detLoc[i][1], &detLoc[i][2] );
    detLoc[i][0]--;
    detLoc[i][1]--;
    detLoc[i][2]--;
  }	

  fclose(fp);


  /* NORMALIZE THE DIRECTION COSINE OF THE SOURCE */
  foo = sqrt(cxi*cxi + cyi*cyi + czi*czi);
  cxi /= foo;
  cyi /= foo;
  czi /= foo;


  /* CALCULATE THE MIN AND MAX PHOTON LENGTH FROM THE MIN AND MAX PROPAGATION TIMES */
  Lmax = maxT * C_VACUUM / tn[1];
  Lmin = minT * C_VACUUM / tn[1];
  stepL = stepT * C_VACUUM / tn[1];


  /* READ IN THE SEGMENTED DATA FILE */
  fp = fopen( segFile, "rb" );
  if( fp==NULL ) {
    printf( "ERROR: The binary image file %s was not found!\n", segFile );
    exit(1);
  }
  tissueType = (char ***)malloc(nxstep*sizeof(char **));
  for( i=0; i<nxstep; i++ ) {
    tissueType[i] = (char **)malloc(nystep*sizeof(char *));
    for( j=0; j<nystep; j++ ) {
      tissueType[i][j] = (char *)malloc(nzstep*sizeof(char));
    }
  }
  for( k=0; k<nzstep; k++ ) {
    for( j=0; j<nystep; j++ ) {
      for( i=0; i<nxstep; i++ ) {
	fscanf( fp, "%c", &tissueType[i][j][k] );
      }
    }
    printf( "%d\n", k );
  }
  fclose(fp);


  /* ALLOCATE SPACE FOR AND INITIALIZE THE PHOTON FLUENCE TO 0 */
  II = (double *)malloc(nIzstep*nIxstep*nIystep*(nTstep+1)*sizeof(double));
  pII = &II[0];
  for( i=0; i<nIzstep; i++) {
    for( j=0; j<nIystep; j++ )
      for( k=0; k<nIxstep; k++ ) 
	for( kk=0; kk<=nTstep; kk++ ) {

	  *pII = 0.;
	  pII++;
	}
  }


  /* NUMBER OF CALLS MADE TO THE RANDOM NUMBER GENERATOR */
  ncall = 0;


  /* MAKE SURE THE SOURCE IS AT AN INTERFACE */
  i = (int)(xi);
  j = (int)(yi);
  k = (int)(zi);
  tissueIndex=tissueType[i][j][k];
  while( tissueIndex!=0 && i>0 && i<nxstep && j>0 && j<nystep &&
	 k>0 && k<nzstep ) {
    xi -= cxi;
    yi -= cyi;
    zi -= czi;
    i = (int)(xi);
    j = (int)(yi);
    k = (int)(zi);
    tissueIndex=tissueType[i][j][k];
  }
  while( tissueIndex==0 ) {
    xi += cxi;
    yi += cyi;
    zi += czi;
    i = (int)(xi);
    j = (int)(yi);
    k = (int)(zi);
    tissueIndex=tissueType[i][j][k];
  }
  

  /* NUMBER PHOTONS EXECUTED SO FAR */
  N = 0;


  /* OPEN A FILE POINTER TO SAVE THE HISTORY INFORMATION */
  sprintf( filenm, "%s.his", argv[1] );
  fp = fopen( filenm, "w" );


  /*********************************************************
     START MIGRATING THE PHOTONS 
     GENERATING PHOTONS UNTIL NUMBER OF PHOTONS EXECUTED
     (N) IS EQUAL TO THE NUMBER TO BE GENERATED (NT) 
   *********************************************************/
  while (N<NT){
    ++N;								

    /* SET THE PHOTON WEIGHT TO 1 AND INITIALIZE PHOTON LENGTH PARAMETERS */
    P2pt = 1.;
    Ltot = 0.;
    Lnext = 1.;
    Lresid = 0.;

    /* INITIALIZE THE LENGTH IN EACH TISSUE TYPE */
    for( i=0; i<=Ntissue; i++ )
      lenTiss[i] = 0.;

    
    /* INITIAL SOURCE POSITION */
    x = xi;
    y = yi;
    z = zi;	
    xold = x; yold = y; zold = z;					
		
    /* INITIAL DIRECTION OF PHOTON */
    c1 = cxi;
    c2 = cyi;
    c3 = czi;
    c1o = c1;
    c2o = c2;
    c3o = c3;


    
    /* CALCULATE SCATTERING LENGTH */
    rnm = ran( &idum, &ncall );

    /* PROPAGATE THE PHOTON */
    Lresid = -log(rnm);
    i = (int)(x);
    j = (int)(y);
    k = (int)(z);
    tissueIndex=tissueType[i][j][k];
    while( Lresid>0. && i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep && tissueIndex!=0 ) {
      
      if( Ltot>Lnext && Ltot>Lmin ) {
	if ( i>=Ixmin && i<=Ixmax && j>=Iymin && j<=Iymax && k>=Izmin && k<=Izmax) {
	  tindex = (int)floor((Ltot-Lmin)/stepL);
	  II[tindex*nIzstep*nIxstep*nIystep+(k-Izmin)*nIxstep*nIystep+(j-Iymin)*nIxstep+(i-Ixmin)] += P2pt;
	}
	Lnext += 1.;	/* THIS INSURES THAT WE DON'T DOUBLE SCORE A PHOTON IN EACH VOXEL */
      }
      
      if( (foo=tmus[tissueIndex])>Lresid ) {
	step = Lresid / foo;
	x += c1*step;
	y += c2*step;
	z += c3*step;
	Ltot += step;
	
	P2pt *= exp(-tmua[tissueIndex] * step );

	lenTiss[tissueIndex] += (float)step;
	Lresid = 0;
      } else {
	x += c1;
	y += c2;
	z += c3;
	Ltot += 1.;

	P2pt *= exp(-tmua[tissueIndex]);

	Lresid -= foo;
	lenTiss[tissueIndex] += (float)1.0;
      }


      i = (int)(x);
      j = (int)(y);
      k = (int)(z);
      if( i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep )
	tissueIndex=tissueType[i][j][k];
    }

		
    if( i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep )
      tissueIndex=tissueType[i][j][k];


    /* LOOP UNTIL TIME EXCEEDS GATE OR PHOTON ESCAPES */
    while ( Ltot<Lmax && i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep && tissueIndex!=0 ) {

      /* CALCULATE THE NEW SCATTERING ANGLE USING HENYEY-GREENSTEIN */
      gg = tg[tissueIndex];

      rnm = ran( &idum, &ncall );
      phi=2*pi*rnm;
      cphi=cos(phi);
      sphi=sin(phi);

      rnm = ran( &idum, &ncall );
      foo = (1. - gg*gg)/(1. - gg + 2.*gg*rnm);
      foo = foo * foo;
      foo = (1. + gg*gg - foo)/(2.*gg);
      theta=acos(foo);
      stheta=sin(theta);
      ctheta=foo;

      c1o = c1;
      c2o = c2;
      c3o = c3;
      if( c3<1. && c3>-1. ) {
	c1 = stheta*(c1o*c3o*cphi - c2o*sphi)/sqrt(1 - c3o*c3o) + c1o*ctheta;
	c2 = stheta*(c2o*c3o*cphi + c1o*sphi)/sqrt(1 - c3o*c3o) + c2o*ctheta;
	c3 = -stheta*cphi*sqrt(1-c3o*c3o)+c3o*ctheta;
      }
      else {
	c1 = stheta*cphi;
	c2 = stheta*sphi;
	c3 = ctheta*c3;
      }
      
      /* CALCULATE SCATTERING LENGTH */
      rnm = ran( &idum, &ncall );
      
      
      /* PROPAGATE THE PHOTON */
      Lresid = -log(rnm);
      i = (int)(x);
      j = (int)(y);
      k = (int)(z);
      tissueIndex=tissueType[i][j][k];
      while( Ltot<Lmax && Lresid>0. && i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep && (tissueIndex=tissueType[i][j][k])!=0 ) {
	
	if( Ltot>Lnext && Ltot>Lmin ) {
	  if ( i>=Ixmin && i<=Ixmax && j>=Iymin && j<=Iymax && k>=Izmin && k<=Izmax) {
	    tindex = (int)floor((Ltot-Lmin)/stepL);
	    II[tindex*nIzstep*nIxstep*nIystep+(k-Izmin)*nIxstep*nIystep+(j-Iymin)*nIxstep+(i-Ixmin)] += P2pt;
	  }
	  Lnext += 1.;
	}

	if( (foo=tmus[tissueIndex])>Lresid ) {
	  step = Lresid / foo;
	  x += c1*step;
	  y += c2*step;
	  z += c3*step;
	  Ltot += step;
	  
	  P2pt *= exp(-tmua[tissueIndex] * step );
	  
	  lenTiss[tissueIndex] += (float)step;
	  Lresid = 0;
	} else {
	  x += c1;
	  y += c2;
	  z += c3;
	  Ltot += 1.;

	  P2pt *= exp(-tmua[tissueIndex]);
	  
	  Lresid -= foo;
	  lenTiss[tissueIndex] += (float)1.0;
	}
	
	i = (int)(x);
	j = (int)(y);
	k = (int)(z);
	

      } /* PROPAGATE PHOTON */
      
    } /* LOOP UNTIL END OF SINGLE PHOTON */

    /* SCORE EXITING PHOTON AND SAVE HISTORY FILES*/
    i = (int)(x);
    j = (int)(y);
    k = (int)(z);
    
    if ( i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep ) {
      tissueIndex=tissueType[i][j][k];
      if( tissueIndex==0 ) {
	if( i>=Ixmin && i<=Ixmax && j>=Iymin && j<=Iymax && k>=Izmin && k<=Izmax ) {		    
	  tindex = (int)floor((Ltot-Lmin)/stepL);
	  VoxN = tindex*nIzstep*nIxstep*nIystep+(k-Izmin)*nIxstep*nIystep+(j-Iymin)*nIxstep+(i-Ixmin);
	  II[VoxN] -= P2pt;
	}
	
	/* LOOP THROUGH NUMBER OF DETECTORS */
	for( ii=0; ii<nDets; ii++ ) {
	  if( abs(i-detLoc[ii][0])<=detRad )
	    if( abs(j-detLoc[ii][1])<=detRad )
	      if( abs(k-detLoc[ii][2])<=detRad ) {
		
		/* WRITE TO THE HISTORY FILE	*/
		ffoo = ii;
		fwrite( &ffoo, sizeof(float), 1, fp );
		for( jj=1; jj<=Ntissue; jj++ ) {
		  fwrite( &lenTiss[jj], sizeof(float), 1, fp );
		}
	      }
	}

	/* IF NO DETECTORS THEN SAVE EXIT POSITION */
	if( nDets==0 ) {
	  ffoo=i; fwrite( &ffoo, sizeof(float), 1, fp );
	  ffoo=j; fwrite( &ffoo, sizeof(float), 1, fp );
	  ffoo=k; fwrite( &ffoo, sizeof(float), 1, fp );
	  for( jj=1; jj<=Ntissue; jj++ ) {
	    fwrite( &lenTiss[jj], sizeof(float), 1, fp );
	  }
	}
	
      }
    }

		
  } /* LOOP UNTIL ALL PHOTONS EXHAUSTED */

  /* CLOSE HISTORY FILE */
  fclose(fp); 
  

  
  /* SAVE FLUENCE DATA */
  sprintf( filenm, "%s.2pt", argv[1]);
  pII = &II[0];
  fp = fopen( filenm, "wb");
  fwrite(pII,sizeof(double),nIxstep*nIystep*nIzstep*nTstep,fp);
  fclose( fp );  

  return 1;
}





