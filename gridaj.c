#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main() 
{
FILE *output;
double *xi;
double *gridx,*gridy;
int i,j;
int n = 4;
double a = -1;
double b = 1;
double h = (b - a) /(n-1);

xi = (double *) malloc(n*sizeof(double));
gridx = (double *) malloc(n*n*sizeof(double));
gridy = (double *) malloc(n*n*sizeof(double));

for (i = 0; i < n; i++) {
    xi[i] = a + i*h;
    printf("xi[%d] = %f \n",i,xi[i]);
    }

for (i = 0; i < n*n; i++){
    for (j = 0; j < n; j++){
        gridx[i*n+j] = xi[j];
        gridy[j*n+i] = xi[j];  
  }
}
  
/*for (i = 0; i < n*n; i++){*/
/*    for (j = 0; j < n; j++){*/
/*        printf("xij[%d] = %f \n",i*n+j,gridx[i*n+j]);*/
/*        printf("xji[%d] = %f \n",j*n+i,gridy[i*n+j]);*/
/*    }*/
/*}*/

output = fopen ( "file.dat", "wt" );
        
if ( !output ){
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "WRITE_RESULTS - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
}

//Write the data.
for ( i = 0; i < n*n; i++ )
for ( j = 0; j < n; j++ ){
    fprintf ( output, "  %24.16g\t%24.16g\t", gridx[i*n+j], gridy[i*n+j]);
}
    return 0;
}
