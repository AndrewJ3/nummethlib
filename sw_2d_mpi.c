#include "sw_mpi.h"

/******************************************************************************/
/*
 Purpose:
 MAIN is the main program for SHALLOW_WATER_2D.
 
 Discussion:
 SHALLOW_WATER_2D approximates the 2D shallow water equations.
 The version of the shallow water equations being solved here is in
 conservative form, and omits the Coriolis force.  The state variables
 are H (the height), UH and VH (the mass velocity).

 The initial conditions are used to specify ( H, UH, VH ) at an equally
 spaced set of finite volumes, and then the Lax-Friedrichs method is used to advance
 the solution until a final time t_final.
 Some input values will result in an unstable calculation that
 quickly blows up.  This is related to the Courant-Friedrichs-Levy
 condition, which requires that DT be small enough, relative to DX and
 the velocity, that information cannot cross an entire cell.
 
 A "reasonable" set of input quantities is
 shallow_water_1d 400 0.002 10.0 0.5
 
 Licensing:
 This code is distributed under the GNU LGPL license.
 
 Modified:
 26 March 2019 by Michal A. Kopera - complete overhaul of the code
 07 April 2019 by Andrew M. Jones -parallelization
 
 Reference:
 John Burkardt,
 Cleve Moler,
 "The Shallow Water Equations",
 Experiments with MATLAB.
 
 */
int main ( int argc, char *argv[] )
{
    double dx;
    double dy;
    double dt;
    double g = 9.81; //[m/s^2] gravitational acceleration
    double *h;
    double *fh;
    double *gh;
    double *hm;
    int i,j, id, id_left, id_right, id_bottom, id_top,id1,id2;
    int nx, ny, nx_loc, ny_loc;
    double t_final;
    double *uh;
    double *fuh;
    double *guh;
    double *uhm;
    double *vh;
    double *fvh;
    double *gvh;
    double *vhm;
    double *x;
    double *y;
    double x_length, time;
	double *buf_left_recv; 
	double *buf_right_recv; 
	double *buf_left_send;
	double *buf_right_send;
    double *buf_up_recv;
	double *buf_down_recv; 
	double *buf_up_send; 
	double *buf_down_send;
    int irank, nproc;
	double t1,t2,tmpi;

	//initialize MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

	//get command line arguments
    getArgs_mpi(&nx, &dt, &x_length, &t_final, argc, argv, irank, MPI_COMM_WORLD);

    //divide data among processors
    int q = (int) sqrt((double) nproc);
    ny=nx;
    nx_loc = nx/q;
    ny_loc = ny/q;
	int n_buf = 3*(ny_loc);
    
	//Allocate space (nx+2)((nx+2) long, to accound for ghosts
    //height array
    h = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    hm = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    fh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    gh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    //x momentum array
    uh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    uhm = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    fuh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    guh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    //y momentum array
    vh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    vhm = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    fvh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
    gvh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
	/*LEFT/RIGHT BUFS*/
    buf_left_recv = ( double * ) malloc ( n_buf * sizeof ( double ) );
    buf_right_recv = ( double * ) malloc ( n_buf * sizeof ( double ) );
    buf_left_send = ( double * ) malloc ( n_buf * sizeof ( double ) );
    buf_right_send = ( double * ) malloc ( n_buf * sizeof ( double ) );
    /*TOP/BOTTOM BUFS*/
    buf_up_recv = ( double * ) malloc ( n_buf * sizeof ( double ) );
    buf_down_recv = ( double * ) malloc ( n_buf * sizeof ( double ) );
    buf_up_send = ( double * ) malloc ( n_buf * sizeof ( double ) );
    buf_down_send = ( double * ) malloc ( n_buf * sizeof ( double ) );
    // location arrays
    x = ( double * ) malloc ( nx_loc * sizeof ( double ) );
    y = ( double * ) malloc ( nx_loc * sizeof ( double ) );
    
    //Define the locations of the nodes and time steps and the spacing.
    dx = x_length / ( double ) ( nx );
    dy = x_length / ( double ) ( nx );
    
    //Apply the initial conditions.
    initial_conditions ( nx_loc, ny_loc, dx, dy, x_length,  x, y, h, uh, vh, irank ,q);
    
    //Write initial condition to a file
    write_results("sw2d_mpi_init.dat",nx,nx_loc,nx,ny_loc,x,y,h,uh,vh, irank,nproc);
    
    double lambda_x = 0.5*dt/dx;
    double lambda_y = 0.5*dt/dy;
    
// --------------------------------------------------------------------------------------
    /* REQUESTING & INITIALIZATION OF BUF_FIELDS (TOP/BOTTOM & LEFT/RIGHT) */
//     MPI_Status status;
    MPI_Request request_h_recv_left;//, request_uh_recv_left,request_vh_recv_left;
    MPI_Request request_h_recv_right;//, request_uh_recv_right,request_vh_recv_right;
    MPI_Request request_h_send_left;//, request_uh_send_left,request_vh_send_left;
    MPI_Request request_h_send_right;//, request_uh_send_right,request_vh_send_right;
        
    MPI_Request request_h_recv_up;//, request_uh_recv_up,request_vh_recv_up;
    MPI_Request request_h_recv_down;//, request_uh_recv_down,request_vh_recv_down;
    MPI_Request request_h_send_up;//, request_uh_send_up,request_vh_send_up;
    MPI_Request request_h_send_down;//, request_uh_send_down,request_vh_send_down;
    
    MPI_Status status;
    
    /*SET-UP CARTESIAN COMMUNICATOR*/
    MPI_Comm Comm_cart;
    int ndims=2;
    int rank_left, rank_right,rank_up,rank_down;
    int dims[ndims];
    int wrap_around[2];
    
    //no periodicity
    wrap_around[0] = 0;
    wrap_around[1] = 0;

    dims[0] = q;
    dims[1] = q;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, 0, &Comm_cart);
    MPI_Cart_shift(Comm_cart, 1 ,1, &rank_left, &rank_right);
    MPI_Cart_shift(Comm_cart, 0 ,1, &rank_down, &rank_up);

//-------------------------------------------------------------------------------------

    time=0.0;
	t1=MPI_Wtime();
    while (time<t_final)
    {
        //  Take a time step
        time=time+dt;
//         printf("time = %f\n",time);

        /*PARALLELIZATION BEGINS*/
//--------------------------------------------------------------
        /*SEND RECVS FOR BUF FIELDS ON TOP/BOTTOM & LEFT/RIGHT*/

        /***********LEFT/RIGHT************/
        if(rank_left>=0){
            //set up receiving from the left
            MPI_Irecv(buf_left_recv, n_buf, MPI_DOUBLE,rank_left,
            101, MPI_COMM_WORLD, &request_h_recv_left);

            j=1;
            for ( i = 0; i < ny_loc; i++ ){
                id=ID_2D(i+1,j,ny_loc);		
				//set up sending to the left
				buf_left_send[i] = h[id]; 
				buf_left_send[i+ny_loc] = uh[id];
				buf_left_send[i+2*ny_loc] = vh[id];
            }
            MPI_Isend(buf_left_send,n_buf,MPI_DOUBLE,rank_left,
            101,MPI_COMM_WORLD,&request_h_send_left); 
        
        }

        if(rank_right>=0){
            MPI_Irecv(buf_right_recv, n_buf, MPI_DOUBLE,rank_right,
            101, MPI_COMM_WORLD, &request_h_recv_right);
            j=nx_loc;
            for ( i = 0; i < ny_loc; i++ ){
                id=ID_2D(i+1,j,ny_loc);
				buf_right_send[i] = h[id]; 
				buf_right_send[i+ny_loc] = uh[id];
				buf_right_send[i+2*ny_loc] = vh[id];
			}
            MPI_Isend(buf_right_send,n_buf,MPI_DOUBLE,rank_right,
            101,MPI_COMM_WORLD,&request_h_send_right); 		
        }

		/**********TOP/BOTTOM**********/

        if(rank_up>=0){           
            MPI_Irecv(buf_up_recv, n_buf, MPI_DOUBLE,rank_up,
            101, MPI_COMM_WORLD, &request_h_recv_up);

            i=ny_loc;
            for ( j = 0; j < nx_loc; j++ ){
                id=ID_2D(i,j+1,nx_loc);

				//set up sending to the top
				buf_up_send[j] = h[id]; 
				buf_up_send[j+ny_loc] = uh[id];
				buf_up_send[j+2*ny_loc] = vh[id];
			}
			MPI_Isend(buf_up_send,n_buf,MPI_DOUBLE,rank_up,
            101,MPI_COMM_WORLD,&request_h_send_up); 
        }

        if(rank_down>=0){
			 //set up receiving from the top
            MPI_Irecv(buf_down_recv, n_buf, MPI_DOUBLE,rank_down,
            101, MPI_COMM_WORLD, &request_h_recv_down);
			
			i=1;
            for ( j = 0; j < nx_loc; j++ ){
                id=ID_2D(i,j+1,nx_loc);

				//set up sending to the bottom
				buf_down_send[j] = h[id]; 
				buf_down_send[j+nx_loc] = uh[id]; 
				buf_down_send[j+2*nx_loc] = vh[id];
            }

			 MPI_Isend(buf_down_send,n_buf,MPI_DOUBLE,rank_down,
            101,MPI_COMM_WORLD,&request_h_send_down); 
		}
//-----------------------------------------------------------       
        //Compute interior fluxes
        for ( i = 1; i < ny_loc+1; i++ )
        for ( j = 1; j < nx_loc+1; j++){
            id=ID_2D(i,j,nx_loc);
            
            fh[id] = uh[id]; //flux for the height equation: u*h
            fuh[id] = uh[id]*uh[id]/h[id] + 0.5*g*h[id]*h[id]; //flux for the momentum equation: u^2*h + 0.5*g*h^2
            fvh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gh[id] = vh[id]; //flux for the height equation: v*h
            guh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gvh[id] = vh[id]*vh[id]/h[id] + 0.5*g*h[id]*h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
        }
        
//---------------------------------------------------------------
        /*********WAIT FOR COMMUNICATION TO COMPLETE***********/
        /**********LEFT/RIGHT***********/
        /****WAIT & RECV****/
        if(rank_left>=0){
	        MPI_Wait(&request_h_recv_left,&status);
			j=0;
			for(i = 0; i < ny_loc ; i++){
                id=ID_2D(i+1,j,ny_loc);
                h[id]=buf_left_recv[i];
                uh[id]=buf_left_recv[i+ny_loc];
                vh[id]=buf_left_recv[i+2*ny_loc];
            }

        }
		/****BOUNDARY_CONDITION****/
        else{ 
            j=0;
   	   		for( i = 1; i < ny_loc+1 ; i++){
                id1=ID_2D(i,j,ny_loc);
                id2=ID_2D(i,1,ny_loc);
                h[id1]=h[id2];
                uh[id1]=-uh[id2];
                vh[id1]=vh[id2];
            }
        }
        /****WAIT & RECV****/
        if(rank_right>=0){
		    MPI_Wait(&request_h_recv_right,&status);
            j=nx_loc+1;
	   		for(i = 0; i < ny_loc ; i++){
	            id=ID_2D(i+1,j,ny_loc);
	            h[id]=buf_right_recv[i];
	            uh[id]=buf_right_recv[i+ny_loc];
	            vh[id]=buf_right_recv[i+2*ny_loc];
	        }
        }
        /****BOUNDARY_CONDITION****/
        else{ 
            j=ny_loc+1;
   	   		for( i = 1; i < ny_loc+1 ; i++){
                id1=ID_2D(i,j,ny_loc);
                id2=ID_2D(i,ny_loc,ny_loc);
                h[id1]=h[id2];
                uh[id1]=-uh[id2];
                vh[id1]=vh[id2];
            }
        }

        /*************TOP/BOTTOM**************/
        /****WAIT & RECV****/
        if(rank_up>=0){
            MPI_Wait(&request_h_recv_up,&status);
            i=nx_loc+1;
	   		for( j = 0; j < nx_loc ; j++){
	            id=ID_2D(i,j+1,nx_loc);
	            h[id]=buf_up_recv[j];
    	        uh[id]=buf_up_recv[j+nx_loc];
        	    vh[id]=buf_up_recv[j+2*nx_loc];
    	    }
        }
        /****BOUNDARY_CONDITION****/
		else{
            i=nx_loc+1;
   	   		for( j= 1; j < nx_loc+1 ; j++){
				id1=ID_2D(i,j,nx_loc);
				id2=ID_2D(nx_loc,j,nx_loc);
                h[id1]=h[id2];
                uh[id1]=uh[id2];
                vh[id1]=-vh[id2];
            }
        }
        /****WAIT & RECV****/
        if(rank_down>=0){
            MPI_Wait(&request_h_recv_down,&status);
			i=0;
   	   		for( j = 0; j < nx_loc ; j++){
                id=ID_2D(i,j+1,nx_loc);
	            h[id]=buf_down_recv[j];
    	        uh[id]=buf_down_recv[j+nx_loc];
        	    vh[id]=buf_down_recv[j+2*ny_loc];
            }
        }
        /****BOUNDARY_CONDITION****/
		else{ 
            i=0;
   	   		for( j = 1; j < nx_loc+1 ; j++){
				id1=ID_2D(i,j,nx_loc);
				id2=ID_2D(1,j,nx_loc);
                h[id1]=h[id2];
                uh[id1]=uh[id2];
                vh[id1]=-vh[id2];
            }
        }
//---------------------------------------------------------------
        
        //Compute ghost fluxes (need ghost values)
        //left ghost
        j=0;
        for ( i = 1; i < ny_loc+1; i++ ){
            id=ID_2D(i,j,nx_loc);
            fh[id] = uh[id];
            fuh[id] = uh[id]*uh[id]/h[id] + 0.5*g*h[id]*h[id];
            fvh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gh[id] = vh[id]; //flux for the height equation: v*h
            guh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gvh[id] = vh[id]*vh[id]/h[id] + 0.5*g*h[id]*h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
        }
        
        //right ghost
        j = nx_loc+1;
        for ( i = 1; i < ny_loc+1; i++ ){
            id=ID_2D(i,j,nx_loc);
            fh[id] = uh[id];
            fuh[id] = uh[id]*uh[id]/h[id] + 0.5*g*h[id]*h[id];
            fvh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gh[id] = vh[id]; //flux for the height equation: v*h
            guh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gvh[id] = vh[id]*vh[id]/h[id] + 0.5*g*h[id]*h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
        }
        
        //bottom ghost
        i = 0;
        for ( j = 1; j < nx_loc+1; j++ ){
            id=ID_2D(i,j,nx_loc);
            fh[id] = uh[id];
            fuh[id] = uh[id]*uh[id]/h[id] + 0.5*g*h[id]*h[id];
            fvh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gh[id] = vh[id]; //flux for the height equation: v*h
            guh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gvh[id] = vh[id]*vh[id]/h[id] + 0.5*g*h[id]*h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
        }
        
        //top ghost
        i = ny_loc+1;
        for ( j = 1; j < nx_loc+1; j++ ){
            id=ID_2D(i,j,nx_loc);
            fh[id] = uh[id];
            fuh[id] = uh[id]*uh[id]/h[id] + 0.5*g*h[id]*h[id];
            fvh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gh[id] = vh[id]; //flux for the height equation: v*h
            guh[id] = uh[id]*vh[id]/h[id]; //flux for the momentum equation: u*v**h
            gvh[id] = vh[id]*vh[id]/h[id] + 0.5*g*h[id]*h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
        }
        
        //Compute updated variables
        for ( i = 2; i < ny_loc; i++ )
        for ( j = 2; j < nx_loc; j++ )
        {
            id=ID_2D(i,j,nx_loc);
            id_left=ID_2D(i,j-1,nx_loc);
            id_right=ID_2D(i,j+1,nx_loc);
            id_bottom=ID_2D(i-1,j,nx_loc);
            id_top=ID_2D(i+1,j,nx_loc);
            hm[id] = 0.25*(h[id_left]+h[id_right]+h[id_bottom]+h[id_top])
            - lambda_x * ( fh[id_right] - fh[id_left] )
            - lambda_y * ( gh[id_top] - gh[id_bottom] );
            uhm[id] = 0.25*(uh[id_left]+uh[id_right]+uh[id_bottom]+uh[id_top])
            - lambda_x * ( fuh[id_right] - fuh[id_left] )
            - lambda_y * ( guh[id_top] - guh[id_bottom] );
            vhm[id] = 0.25*(vh[id_left]+vh[id_right]+vh[id_bottom]+vh[id_top])
            - lambda_x * ( fvh[id_right] - fvh[id_left] )
            - lambda_y * ( gvh[id_top] - gvh[id_bottom] );
        }
      
        //Compute updated variables at the edges (need ghost values)
        //left boundary
        j=1;
        for ( i = 2; i < ny_loc; i++ ){
            id=ID_2D(i,j,nx_loc);
            id_left=ID_2D(i,j-1,nx_loc);
            id_right=ID_2D(i,j+1,nx_loc);
            id_bottom=ID_2D(i-1,j,nx_loc);
            id_top=ID_2D(i+1,j,nx_loc);
            hm[id] = 0.25*(h[id_left]+h[id_right]+h[id_bottom]+h[id_top])
            - lambda_x * ( fh[id_right] - fh[id_left] )
            - lambda_y * ( gh[id_top] - gh[id_bottom] );
            uhm[id] = 0.25*(uh[id_left]+uh[id_right]+uh[id_bottom]+uh[id_top])
            - lambda_x * ( fuh[id_right] - fuh[id_left] )
            - lambda_y * ( guh[id_top] - guh[id_bottom] );
            vhm[id] = 0.25*(vh[id_left]+vh[id_right]+vh[id_bottom]+vh[id_top])
            - lambda_x * ( fvh[id_right] - fvh[id_left] )
            - lambda_y * ( gvh[id_top] - gvh[id_bottom] );
        }
        //right boundary
        j=nx_loc;
        for ( i = 2; i < ny_loc; i++ ){
            id=ID_2D(i,j,nx_loc);
            id_left=ID_2D(i,j-1,nx_loc);
            id_right=ID_2D(i,j+1,nx_loc);
            id_bottom=ID_2D(i-1,j,nx_loc);
            id_top=ID_2D(i+1,j,nx_loc);
            hm[id] = 0.25*(h[id_left]+h[id_right]+h[id_bottom]+h[id_top])
            - lambda_x * ( fh[id_right] - fh[id_left] )
            - lambda_y * ( gh[id_top] - gh[id_bottom] );
            uhm[id] = 0.25*(uh[id_left]+uh[id_right]+uh[id_bottom]+uh[id_top])
            - lambda_x * ( fuh[id_right] - fuh[id_left] )
            - lambda_y * ( guh[id_top] - guh[id_bottom] );
            vhm[id] = 0.25*(vh[id_left]+vh[id_right]+vh[id_bottom]+vh[id_top])
            - lambda_x * ( fvh[id_right] - fvh[id_left] )
            - lambda_y * ( gvh[id_top] - gvh[id_bottom] );
        }
        
        //bottom boundary
        i=1;
        for ( j = 1; j < nx_loc+1; j++ ){
            id=ID_2D(i,j,nx_loc);
            id_left=ID_2D(i,j-1,nx_loc);
            id_right=ID_2D(i,j+1,nx_loc);
            id_bottom=ID_2D(i-1,j,nx_loc);
            id_top=ID_2D(i+1,j,nx_loc);
            hm[id] = 0.25*(h[id_left]+h[id_right]+h[id_bottom]+h[id_top])
            - lambda_x * ( fh[id_right] - fh[id_left] )
            - lambda_y * ( gh[id_top] - gh[id_bottom] );
            uhm[id] = 0.25*(uh[id_left]+uh[id_right]+uh[id_bottom]+uh[id_top])
            - lambda_x * ( fuh[id_right] - fuh[id_left] )
            - lambda_y * ( guh[id_top] - guh[id_bottom] );
            vhm[id] = 0.25*(vh[id_left]+vh[id_right]+vh[id_bottom]+vh[id_top])
            - lambda_x * ( fvh[id_right] - fvh[id_left] )
            - lambda_y * ( gvh[id_top] - gvh[id_bottom] );
        }
        
        //top boundary
        i=ny_loc;
        for ( j = 1; j < nx_loc+1; j++ ){
            id=ID_2D(i,j,nx_loc);
            id_left=ID_2D(i,j-1,nx_loc);
            id_right=ID_2D(i,j+1,nx_loc);
            id_bottom=ID_2D(i-1,j,nx_loc);
            id_top=ID_2D(i+1,j,nx_loc);
            hm[id] = 0.25*(h[id_left]+h[id_right]+h[id_bottom]+h[id_top])
            - lambda_x * ( fh[id_right] - fh[id_left] )
            - lambda_y * ( gh[id_top] - gh[id_bottom] );
            uhm[id] = 0.25*(uh[id_left]+uh[id_right]+uh[id_bottom]+uh[id_top])
            - lambda_x * ( fuh[id_right] - fuh[id_left] )
            - lambda_y * ( guh[id_top] - guh[id_bottom] );
            vhm[id] = 0.25*(vh[id_left]+vh[id_right]+vh[id_bottom]+vh[id_top])
            - lambda_x * ( fvh[id_right] - fvh[id_left] )
            - lambda_y * ( gvh[id_top] - gvh[id_bottom] );
        }
        
        //update interior state variables
        for (i = 1; i < ny_loc+1; i++)
        for (j = 1; j < nx_loc+1; j++){
            id=ID_2D(i,j,nx_loc);
            h[id] = hm[id];
            uh[id] = uhm[id];
            vh[id] = vhm[id];
        }
    }
	/*************PARALLELIZATION ENDS HERE*************/
	/**TIMING**/
	t2=MPI_Wtime();
	tmpi=t2-t1;
    if(irank==0){
		printf("%d %d %f \n",nx,nproc,tmpi);
	}
	
	// Write data to file
    write_results("sw2d_mpi_final.dat",nx,nx_loc,nx,ny_loc,x,y,h,uh,vh, irank,nproc);
    
    //Free memory.
    free(h);
    free(uh);
    free(vh);
    free(x);
    free(y);
    free(uhm);
    free(vhm);
    free(fuh);
    free(guh);
    free(fvh);
    free(gvh);
    free(fh);
    free(gh);
    free(buf_left_recv);
    free(buf_right_recv);
    free(buf_left_send);
    free(buf_right_send);  
    free(buf_up_recv);
    free(buf_down_recv);
    free(buf_up_send); 
    free(buf_down_send);
    MPI_Finalize();
    return 0;
}
