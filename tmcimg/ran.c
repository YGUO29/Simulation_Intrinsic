#include <stdio.h>
#include <math.h>

static long idum2;

float ran( int *idum, int *ncall )
{
  float rn;
  float ran3( int *idum );
  float ran2( long *idum2 );

  *ncall = *ncall + 1;
  if( *ncall==1 ) idum2 = *idum*0.79;
  if( fmod(*ncall,1e6)==0 )  *idum=-1111-99999*ran2(&idum2);
  rn = ran3(idum);

  if( rn<1e-4 )  rn=1.e-4*ran3(idum);
  if( rn==0 )  rn=1e-12;

  return rn;
}
