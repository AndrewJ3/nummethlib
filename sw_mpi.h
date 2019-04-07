# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <mpi.h>
#define ID_2D(i,j,nx) ((i)*(nx+2)+(j))
/******************************************************************************/

void compute_flux_2d_point(double fh, double fuh, double fvh, 
                           double gh, double guh, double gvh,
                           double h, double uh, double vh, double g)
{
    fh = uh;
    fuh = uh*uh/h + 0.5*g*h*h;
    fvh = uh*vh/h; //flux for the momentum equation: u*v**h
    gh = vh; //flux for the height equation: v*h
    guh = uh*vh/h; //flux for the momentum equation: u*v**h
    gvh = vh*vh/h + 0.5*g*h*h; //flux for the momentum equation: v^2*h + 0.5*g*h^2
}
void initial_conditions ( int nx, int ny, double dx, double dy,  double x_length, double x[],double y[], double h[], double uh[] ,double vh[], int irank, int q )

/******************************************************************************/

{
  int i,j, id, id1;

  for ( i = 1; i < nx+1; i++ )
    {
      x[i-1] = -x_length/2+dx/2+(irank%q*nx+i-1)*dx;
      y[i-1] = -x_length/2+dy/2+(floor(irank/q)*ny+i-1)*dy;
    }

  for ( i = 1; i < nx+1; i++ )
    for( j = 1; j < ny+1; j++)
      {
	double xx = x[j-1];
	double yy = y[i-1];
	id=ID_2D(i,j,nx);
	h[id] = 1.0 + 0.4*exp ( -5 * ( xx*xx + yy*yy) )-0.1*cos(xx+yy);
      }
  
  for ( i = 1; i < nx+1; i++ )
    for( j = 1; j < ny+1; j++)
      {
	id=ID_2D(i,j,nx);
	uh[id] = 0.0;
	vh[id] = 0.0;
      }

  //set boundaries
  //bottom
  i=0;
  for( j = 1; j < nx+1; j++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i+1,j,nx);

      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }

  //top
  i=nx+1;
  for( j = 1; j < nx+1; j++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i-1,j,nx);

      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }

  //left
  j=0;
  for( i = 1; i < ny+1; i++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i,j+1,nx);

      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }

  //right
  j=nx+1;
  for( i = 1; i < ny+1; i++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i,j-1,nx);

      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }

  return;
}
/******************************************************************************/


void write_results ( char *output_filename, int nx, int nx_loc, int ny, int ny_loc, double x[], double y[], double h[], double uh[], double vh[], int irank, int nproc )
/******************************************************************************/

{
    int i,j, id;
    FILE *output;
    double *uh_print;
    double *vh_print;
    double *h_print;
    double *x_print;
    double *y_print;
    double *uh_send;
    double *vh_send;
    double *h_send;
    double *x_send;
    double *y_send;
    int recvcount[nproc];
    int displs[nproc];
    
    //allocate send arrays
    x_send = malloc(sizeof(double)*nx_loc*ny_loc);
    y_send = malloc(sizeof(double)*nx_loc*ny_loc);
    h_send = malloc(sizeof(double)*nx_loc*ny_loc);
    uh_send = malloc(sizeof(double)*nx_loc*ny_loc);
    vh_send = malloc(sizeof(double)*nx_loc*ny_loc);
    
    //allocate global array
    if(irank==0){
        x_print = malloc(sizeof(double)*nx*ny);
        y_print = malloc(sizeof(double)*nx*ny);
        h_print = malloc(sizeof(double)*nx*ny);
        uh_print = malloc(sizeof(double)*nx*ny);
        vh_print = malloc(sizeof(double)*nx*ny);
    }
    
    //prepare to send data
    for(i=0; i<ny_loc; i++)
    for(j=0; j<nx_loc; j++){
        id=ID_2D(i+1,j+1,nx_loc);
        x_send[i*nx_loc+j] = x[j];
        y_send[i*nx_loc+j] = y[i];
        h_send[i*nx_loc+j] = h[id];
        uh_send[i*nx_loc+j] = uh[id];
        vh_send[i*nx_loc+j] = vh[id];
    }
    
    //gather data to one process
    // mind that I am not accounting for extra points this time
    // so could have used MPI_Gather instead of MPI_Gatherv
    for(i=0; i<nproc; i++)
    recvcount[i] = nx*ny/nproc;
    
    displs[0] = 0;
    for(i=1; i<nproc; i++)
    displs[i]=displs[i-1]+recvcount[i];
    
    MPI_Gatherv(uh_send,nx_loc*ny_loc,MPI_DOUBLE,uh_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gatherv(vh_send,nx_loc*ny_loc,MPI_DOUBLE,vh_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gatherv(h_send,nx_loc*ny_loc,MPI_DOUBLE,h_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gatherv(x_send,nx_loc*ny_loc,MPI_DOUBLE,x_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gatherv(y_send,nx_loc*ny_loc,MPI_DOUBLE,y_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    if(irank==0){
        //Open the file.
        output = fopen ( output_filename, "wt" );
        
        if ( !output ){
            fprintf ( stderr, "\n" );
            fprintf ( stderr, "WRITE_RESULTS - Fatal error!\n" );
            fprintf ( stderr, "  Could not open the output file.\n" );
            exit ( 1 );
        }
        
        //Write the data.
        for ( i = 0; i < ny; i++ )
        for ( j = 0; j < nx; j++ ){
            fprintf ( output, "  %24.16g\t%24.16g\t%24.16g\t %24.16g\t %24.16g\n", x_print[i*nx+j], y_print[i*nx+j],h_print[i*nx+j], uh_print[i*nx+j], vh_print[i*nx+j]);
        }
        
        //Close the file.
        fclose ( output );
        
        //Clean-up
        free(vh_print);
        free(uh_print);
        free(h_print);
        free(y_print);
        free(x_print);
        free(vh_send);
        free(uh_send);
        free(h_send);
        free(x_send);
        free(y_send);
    }
    return;
}
/******************************************************************************/

void getArgs_mpi(int *nx, double *dt, double *x_length, double *t_final, int argc, char *argv[], int irank, MPI_Comm comm)
{
    
    if(irank==0){
        /*
         Get the quadrature file root name:
         */
        if ( argc <= 1 ){
            *nx = 401;
        }else{
            *nx = atoi ( argv[1] );
        }
        
        if ( argc <= 2 ){
            *dt = 0.002;
        }else{
            *dt = atof ( argv[2] );
        }
        
        if ( argc <= 3 ){
            *x_length = 10.0;
        }else{
            *x_length = atof ( argv[3] );
        }
        
        if ( argc <= 4 ){
            *t_final = 0.5;
        }else{
            *t_final = atof ( argv[4] );
        }
    }
    
    MPI_Bcast(nx,1,MPI_INT,0,comm);
    MPI_Bcast(dt,1,MPI_DOUBLE,0,comm);
    MPI_Bcast(x_length,1,MPI_DOUBLE,0,comm);
    MPI_Bcast(t_final,1,MPI_DOUBLE,0,comm);
}
