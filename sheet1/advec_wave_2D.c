#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

static const int halo_width= 2;    // fill in halo_width
static int rank;

const float pi=3.14159;
const float u_x= 1, u_y= 0, c_amp=2;  // choose velocity components and amplitude of initial condition
const float cdt=.3;               // safety factor for timestep (experiment!)
static float dx, dy;              // grid spacings

float ugrad_upw(int i, int j, int ny, float data[][ny]){

    // u.grad operator with upwinding acting on field in data at point i,j.
    //
    const float coeff[]={-3./2.,4./2.,-1./2.};
    float sum_x=0., sum_y=0.;

    int inc = -copysign(1.0, u_x);
    for (int k=0; k<=halo_width; k++){
        sum_x += coeff[k]*data[i+inc*k][j];
    }
    sum_x *= fabs(u_x)/dx;

    inc = -copysign(1.0, u_y);
    for (int k=0; k<=halo_width; k++){
        sum_y += coeff[k]*data[i][j+inc*k];
    }
    sum_y *= fabs(u_y)/dy;

    return sum_x + sum_y;
}

int* find_proc(int rank,int ipx, int ipy,int nprocx,int nprocy,int nprocs)
{
    // Implement finding process rank from coordinates ipx, ipy in process grid!
    // up,down,left and right in order for rank_neighbours
    int *rank_neighbours=(int*)malloc(4*sizeof(int));

    rank_neighbours[0] = rank - nprocx;
    if (rank_neighbours[0]<0) rank_neighbours[0] = ipx + (nprocx*(nprocy-1));

    rank_neighbours[1] = rank + nprocx;
    if (rank_neighbours[1]>nprocs-1) rank_neighbours[1] = ipx;

    rank_neighbours[2] = rank - 1;
    if (rank_neighbours[2]<0 || rank_neighbours[2]<(ipy*nprocx)) rank_neighbours[2] = (ipy*nprocx-1) + nprocx;

    rank_neighbours[3] = rank + 1;
    if (rank_neighbours[3]>nprocx*(ipy+1)-1) rank_neighbours[3] = rank_neighbours[3]-nprocx;
    
    return rank_neighbours;
}

int* find_proc_coords(int rank, int nprocx, int nprocy)
{
    // Implement finding process coordinates ipx, ipy in process grid from process rank!
    // printf("The rank is %d\n",rank);
    // printf("The nprocx is %d\n",nprocx);
    // printf("The nprocy is %d\n",nprocy);

    // The template passes 3 arguments, but either the second or the last one can be optional depending on how to arrage the ranks 
    
    int *rank_coords = (int*)malloc(2*sizeof(int));

    rank_coords[0] = (int)rank % nprocx;
    rank_coords[1] = (int)rank / nprocx;

    return rank_coords;
}

void initcond(int nx, int ny, float x[], float y[], float data[][ny+2*halo_width])
{
    // Initialisation of field in data: harmonic function in x (can be modified to a harmonic in y or x and y):
    for (int ix = halo_width; ix < halo_width+nx; ++ix)
    {
        for (int iy = halo_width; iy < halo_width+ny; ++iy)
        {
            data[ix][iy] = c_amp*sin((double) x[ix]);
            // other choices:
            // data[ix][iy] = c_amp*sin((double) y[iy]);
            // data[ix][iy] = c_amp*sin((double) x[ix])*sin((double) y[iy]);
        }
    }
}

void rhs(const int xrange[2], const int yrange[2], int ny, float data[][ny], float d_data[][ny])
{
    //Right-hand side d_data of pde for field in data for a subdomain defined by xrange, yrange:
    int ix,iy;

    for (ix = xrange[0]; ix < xrange[1]; ++ix)
        for (iy = yrange[0]; iy < yrange[1]; ++iy)
        {
            d_data[ix][iy] = ugrad_upw(ix, iy, ny, data);
        }
}

FILE*
get_file_ptr(const char* prefix, const int pid)
{
    char name[4098];
    sprintf(name,"%s%d",prefix,pid);
    return fopen(name,"w");
}
int main(int argc, char** argv)
{
    int nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Check compatibility of argv parameters!
    int nprocx = atoi(argv[1]),nprocy = atoi(argv[2]); 

    int *proc_coords = find_proc_coords(rank,nprocx,nprocy);
    int ipx=proc_coords[0], ipy=proc_coords[1];

    if (nprocs != nprocx*nprocy) 
    {
        if (rank==0) printf("nprocs != nprcox*nprocy - no meaningful simulation!");
        exit(1);
    }

    // Find neighboring processes!
    int *rank_neighbours = find_proc(rank,ipx, ipy, nprocx,nprocy,nprocs);
    int up=rank_neighbours[0],down=rank_neighbours[1],left=rank_neighbours[2],right=rank_neighbours[3];

    // if (rank == 0)
    // {
    //     printf("For rank %d:\n",rank);
    //     printf("The rank of up is %d\n",up);
    //     printf("The rank of down is %d\n",down);
    //     printf("The rank of left is %d\n",left);
    //     printf("The rank of right is %d\n",right);
    // }

    int domain_nx = atoi(argv[3]),                 // number of gridpoints in x direction
        subdomain_nx = domain_nx / nprocx,                            // subdomain x-size w/o halos
        subdomain_mx = subdomain_nx + 2*halo_width;                               //                  with halos

    int domain_ny = atoi(argv[4]),                 // number of gridpoints in y direction
        subdomain_ny = domain_ny / nprocy,                            // subdomain y-size w/o halos
        subdomain_my = subdomain_ny + 2*halo_width;                               //                  with halos

    // if (rank==0)
    // {
    //     printf("The domain_nx is %d\n",domain_nx);
    //     printf("The domain_ny is %d\n",domain_ny);
    //     printf("The subdomain_nx is %d\n",subdomain_nx);
    //     printf("The subdomain_ny is %d\n",subdomain_ny);
    //     printf("The subdomain_mx is %d\n",subdomain_mx);
    //     printf("The subdomain_my is %d\n",subdomain_my);
    // }

    float data[subdomain_mx][subdomain_my], d_data[subdomain_mx][subdomain_my];

    float xextent=2.*pi, yextent=2.*pi;            // domain has extents 2 pi x 2 pi

    // Set grid spacings dx, dy:
    dx=xextent/domain_nx, dy=yextent/domain_ny;

    float x[subdomain_mx], y[subdomain_my];
    int ix, iy;

    // Populate grid coordinate arrays x,y (equidistant): 
    for (ix=0;ix<subdomain_mx; ix++) x[ix] = (ipx*subdomain_nx - halo_width + ix + 0.5)*dx;
    for (iy=0;iy<subdomain_my; iy++) y[iy] = (ipy*subdomain_ny - halo_width + iy + 0.5)*dy;
    
    // Initialisation of data.
    
    initcond(subdomain_nx, subdomain_ny, x, y, data);

    // Think about convenient data types to access non-contiguous portions of array data!
    float buffer_up[subdomain_mx][subdomain_my];
    float buffer_down[subdomain_mx][subdomain_my];
    float buffer_left[subdomain_mx][subdomain_my];
    float buffer_right[subdomain_mx][subdomain_my];

    for (int i=0;i<subdomain_mx;i++)
    {
        for (int j=0;j<subdomain_my;j++)
        {
            buffer_up[i][j]=0;
            buffer_down[i][j]=0;
            buffer_left[i][j]=0;
            buffer_right[i][j]=0;
        }
    }

    //Create MPI Windows for one-sided communication 
    MPI_Win win;
    MPI_Win_create(&data,subdomain_mx*subdomain_my*sizeof(float),sizeof(float),MPI_INFO_NULL,MPI_COMM_WORLD,&win);

    int sizes_up_down[2] = {subdomain_mx, subdomain_my};      
    int subsizes_up_down[2] = {halo_width, subdomain_nx};                 

    MPI_Datatype subarray_type_up;
    int start_lower_part_of_upper_rank[2] = {subdomain_ny-halo_width, halo_width};            
    MPI_Type_create_subarray(2, sizes_up_down, subsizes_up_down, start_lower_part_of_upper_rank, MPI_ORDER_C, MPI_FLOAT, &subarray_type_up);
    MPI_Type_commit(&subarray_type_up);

    MPI_Datatype subarray_type_down;
    int start_upper_part_of_down_rank[2] = {halo_width, halo_width};                            
    MPI_Type_create_subarray(2, sizes_up_down, subsizes_up_down, start_upper_part_of_down_rank, MPI_ORDER_C, MPI_FLOAT, &subarray_type_down);
    MPI_Type_commit(&subarray_type_down);

    int sizes_left_right[2] = {subdomain_mx, subdomain_my};      
    int subsizes_left_right[2] = {subdomain_ny, halo_width};                 

    MPI_Datatype subarray_type_left;
    int start_right_parts_of_left_rank[2] = {halo_width, subdomain_nx-halo_width};          
    MPI_Type_create_subarray(2, sizes_left_right, subsizes_left_right, start_right_parts_of_left_rank, MPI_ORDER_C, MPI_FLOAT, &subarray_type_left);
    MPI_Type_commit(&subarray_type_left);

    MPI_Datatype subarray_type_right;
    int start_left_parts_of_right_rank[2] = {halo_width, halo_width};                           
    MPI_Type_create_subarray(2, sizes_left_right, subsizes_left_right, start_left_parts_of_right_rank, MPI_ORDER_C, MPI_FLOAT, &subarray_type_right);
    MPI_Type_commit(&subarray_type_right);

    unsigned int iterations = atoi(argv[5]);       // number of iterations=timesteps

    if (u_x==0 && u_y==0) {
      if (rank==0) printf("velocity=0 - no meaningful simulation!");
      exit(1);
    }

    // CFL condition for timestep:
    float dt = cdt*(u_x==0 ? (u_y==0 ? 0 : dy/u_y) : (u_y==0 ? dx/u_x : fmin(dx/u_x,dy/u_y)));

    float t=0.;

    // Consider proper synchronization measures!
    MPI_Win_fence(0,win);

    // Initialize timing!
    double start = MPI_Wtime();

    // Setup domain bounds
    int ixstart=halo_width, ixstop=subdomain_mx-halo_width, iystart=halo_width, iystop=subdomain_my-halo_width;

    // Construct file name for data chunk of process.
    FILE* fptr_computational = get_file_ptr("field_chunk_computational_",rank);

    int range_x[2], range_y[2];

    for (unsigned int iter = 0; iter < iterations; ++iter)
    {
        // Get the data from neighbors!
    
        // if (rank!=3)
        // {
        //     printf("Data of rank %d is:\n",rank);
        //     for (int i = 0; i < subdomain_mx;i++) {
        //         for (int j = 0; j < subdomain_my;j++) {
        //             printf("%f ", data[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        //     printf("\n");
        // }
        
        MPI_Win_fence(0,win); 
        
        MPI_Get(&data[0][subdomain_nx], 1, subarray_type_right, right, 0, 1, subarray_type_right, win);
        MPI_Get(&data[subdomain_ny][0], 1, subarray_type_down, down, 0, 1, subarray_type_down, win);
        
        MPI_Get(&data[0][subdomain_nx-subdomain_mx+halo_width], 1, subarray_type_left, left, 0, 1, subarray_type_left, win);
        MPI_Get(&data[subdomain_ny-subdomain_my+halo_width][0], 1, subarray_type_up, up, 0, 1, subarray_type_up, win);

        // MPI_Get(&buffer_right[0][subdomain_nx], 1, subarray_type_right, right, 0, 1, subarray_type_right, win);
        // MPI_Get(&buffer_down[subdomain_ny][0], 1, subarray_type_down, down, 0, 1, subarray_type_down, win);
        
        // MPI_Get(&buffer_left[0][subdomain_nx-subdomain_mx+halo_width], 1, subarray_type_left, left, 0, 1, subarray_type_left, win);
        // MPI_Get(&buffer_up[subdomain_ny-subdomain_my+halo_width][0], 1, subarray_type_up, up, 0, 1, subarray_type_up, win);


        // if (rank == 0 )
        // {
        //     printf("The content of buffer_right is :\n");
        //     for (int i=0;i<subdomain_mx;i++)
        //     {
        //         for (int j=0;j<subdomain_my;j++)
        //         {
        //             printf("%f\t",buffer_right[i][j]);
        //         }
        //         printf("\n\n");
        //     }

        //     printf("The content of buffer_down is :\n");
        //     for (int i=0;i<subdomain_mx;i++)
        //     {
        //         for (int j=0;j<subdomain_my;j++)
        //         {
        //             printf("%f\t",buffer_down[i][j]);
        //         }
        //         printf("\n\n");
        //     }

        //     printf("The content of buffer_left is :\n");
        //     for (int i=0;i<subdomain_mx;i++)
        //     {
        //         for (int j=0;j<subdomain_my;j++)
        //         {
        //             printf("%f\t",buffer_left[i][j]);
        //         }
        //         printf("\n\n");
        //     }

        //     printf("The content of buffer_up is :\n");
        //     for (int i=0;i<subdomain_mx;i++)
        //     {
        //         for (int j=0;j<subdomain_my;j++)
        //         {
        //             printf("%f\t",buffer_up[i][j]);
        //         }
        //         printf("\n\n");
        //     }


        //     printf("Data of rank %d is:\n",rank);
        //     for (int i = 0; i < subdomain_mx;i++) {
        //         for (int j = 0; j < subdomain_my;j++) {
        //             printf("%f ", data[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        //     printf("\n");
        // }

        // Compute rhs. Think about concurrency of computation and data fetching by MPI_Get!
        range_x[0] = 2*halo_width;
        range_x[1] = subdomain_nx-halo_width;
        range_y[0] = 2*halo_width;
        range_y[1] = subdomain_ny-halo_width;

        rhs(range_x,range_y,subdomain_my,data,d_data);

        // Data arrived -> compute stencils in all points that *are* affected by halo points.
        MPI_Win_fence(0,win);

        range_x[0] = halo_width;
        range_x[1] = 2*halo_width;
        range_y[0] = halo_width;
        range_y[1] = 2*halo_width;

        rhs(range_x,range_y,subdomain_my,data,d_data);

        range_x[0] = subdomain_nx-halo_width;
        range_x[1] = subdomain_nx;
        range_y[0] = subdomain_ny-halo_width;
        range_y[1] = subdomain_ny;

        rhs(range_x,range_y,subdomain_my,data,d_data);

        // Update field in data using rhs in d_data (Euler's method):
        for (ix = ixstart; ix < ixstop; ++ix)
            for (iy = iystart; iy < iystop; ++iy)
            {
                data[ix][iy] += dt*d_data[ix][iy];
            }
        t = t+dt;

        // Output solution for checking/visualisation with choosable cadence!
        if (u_y==0.)
        {
            for (int ix=ixstart; ix < ixstop; ++ix) fprintf(fptr_computational,"%f ",data[ix][iystart]);
        }
        
        else if(u_x == 0.)
        {
            for (int iy=iystart; iy < iystop; ++iy) fprintf(fptr_computational,"%f ",data[ixstart][iy]);
        }
        
        else
        {
            for (int ix=ixstart; ix < ixstop; ++ix) 
                    for (int iy=iystart; iy < iystop; ++iy)
                        fprintf(fptr_computational,"%f ",data[ix][iy]);
        }
        fprintf(fptr_computational,"\n");

    }

    MPI_Type_free(&subarray_type_up);
    MPI_Type_free(&subarray_type_down);
    MPI_Type_free(&subarray_type_left);
    MPI_Type_free(&subarray_type_right);

    // // Finalize timing!
    double end = MPI_Wtime();

    // analytic solution in array data_an:
    float data_an[iterations][subdomain_mx][subdomain_my], xshift[subdomain_mx], yshift[subdomain_my];

    // Construct file name for data chunk of process.
    FILE* fptr_analytical = get_file_ptr("field_chunk_analytical_",rank);

    t=0.;
    for (int iter=0; iter<iterations; iter++) {
        for (ix=0;ix<subdomain_mx;ix++) xshift[ix] = x[ix] - u_x*t;
        for (iy=0;iy<subdomain_my;iy++) yshift[iy] = y[iy] - u_y*t;

        initcond(subdomain_nx, subdomain_ny, xshift, yshift, (float (*)[subdomain_my]) &data_an[iter][0][0]);

        if (u_y==0.)
        {
            for (int ix=ixstart; ix < ixstop; ++ix) fprintf(fptr_analytical,"%f ",data_an[iter][ix][iystart]);
        }
        
        else if(u_x == 0.)
        {
            for (int iy=iystart; iy < iystop; ++iy) fprintf(fptr_analytical,"%f ",data_an[iter][ixstart][iy]);
        }
        
        else
        {
            for (int ix=ixstart; ix < ixstop; ++ix) 
                    for (int iy=iystart; iy < iystop; ++iy)
                        fprintf(fptr_analytical,"%f ",data_an[iter][ix][iy]);
        }
        fprintf(fptr_analytical,"\n");
        t += dt;
    }
    fclose(fptr_analytical);

    MPI_Win_free(&win);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

