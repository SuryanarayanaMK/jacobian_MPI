/** check test1.cpp program for demonstration 

Things yet to be done is 
1.  import the matrix
2.  check for the index of x_old since it's written from 0

**/

#include<iostream>
#include<time.h>
#include<stdio.h>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include "mpi.h"
#define N 32
#define NU (N-2)*(N-2)
#define MAX 1000
#define MASTER 0
void jacobian(int offset,int chunksize,double *a,double *f,double *temp,double *x_old);

using namespace std;

ifstream a_ptr;


int main(){

 int rank,size;
 double start,finish;
 int offset,chunksize;
 int i,j,n;
 int dest;
 int tag1,tag2,tag3,tag4; 
 double *a;
 double *f;
 double *temp;
 double *x_old;
 double t1;
 double t2;    
 
     a = (double *) calloc((NU+1)*(NU+1),sizeof(double));
 x_old = (double *) calloc(NU+1,sizeof(double));
     f = (double *) calloc(NU+1,sizeof(double));


 MPI_Init(NULL,NULL);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);

 tag1 = 1;
 tag2 = 2;
 tag3 = 3;
 tag4 = 4;

 MPI_Status status; 
 chunksize = NU/size;

 temp = (double*) calloc(chunksize,sizeof(double));
 
 start  = clock();

 for(n=0;n<MAX;n++){

    if(n==0){
          
        if(rank == MASTER){
             
             a_ptr.open("matrix_disp.dat",ios::in);
             for(j=0;j<NU;j++){
                for(i=1;i<=NU;i++){
                   a_ptr >>  a[i+j*NU] ;
                }
                   a_ptr >>  f[j+1];
                   x_old[j+1] = 0.5;
             }      

             a_ptr.close();
 
            offset = chunksize;
            for(dest=1;dest<size;dest++){       
              MPI_Send(&offset,1,MPI_INT,dest,tag1,MPI_COMM_WORLD);
              MPI_Send(&f[1],NU,MPI_DOUBLE,dest,tag2,MPI_COMM_WORLD); 
              MPI_Send(&x_old[1],NU,MPI_DOUBLE,dest,tag3,MPI_COMM_WORLD);
              MPI_Send(&a[1+offset*NU],chunksize*NU,MPI_DOUBLE,dest,tag4,MPI_COMM_WORLD);
              offset += chunksize;
            }
        } // end of MASTER
        else{
              MPI_Recv(&offset,1,MPI_INT,0,tag1,MPI_COMM_WORLD,&status);
              MPI_Recv(&f[1],NU,MPI_DOUBLE,0,tag2,MPI_COMM_WORLD,&status); 
              MPI_Recv(&x_old[1],NU,MPI_DOUBLE,0,tag3,MPI_COMM_WORLD,&status);
              MPI_Recv(&a[1+offset*NU],chunksize*NU,MPI_DOUBLE,0,tag4,MPI_COMM_WORLD,&status);
        }
   } // end of n==0
  

                MPI_Bcast(&x_old[1],NU,MPI_DOUBLE,0,MPI_COMM_WORLD);
                
                MPI_Scatter(&x_old[1],chunksize,MPI_DOUBLE,&temp[0],chunksize,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
                if(rank==0){
                     offset = 0;
                }        

                jacobian(offset,chunksize,a,f,temp,x_old);

                MPI_Allgather(&temp[0],chunksize,MPI_DOUBLE,&x_old[1],chunksize,MPI_DOUBLE,MPI_COMM_WORLD);                          

               } // end of n i.e iterations

		  if(rank == 1){
                    for(i=1;i<=NU;i++)
                      cout << "iteration  " << n << "   " << i << "     "  << x_old[i] << endl; 
                   }

// MPI_Barrier(MPI_COMM_WORLD);
  finish = clock();

// if(rank==0){
      cout << "total time is " << (finish-start) << " for " << rank << "  processor " << endl;
// }
 MPI_Finalize(); 

 return 0;
}

void jacobian(int offset,int chunksize,double *a,double *f,double *temp,double *x_old){

    int i,j;
    double sum = 0;
    int count = 0;
    for(j=offset;j<offset+chunksize;j++){
        sum = 0;
        for(i=1;i<=NU;i++){
            if((j+1)!=i)
            sum += a[i+j*NU]*x_old[i];
        }
            temp[count] = (f[j+1]-sum)/a[(j+1)+j*(N-2)*(N-2)];
            count++;
    }        
    
}   // end of jacobian
/*
             if(rank==1){
                 for(j=0;j<NU;j++){
                    for(i=1;i<=NU;i++){
                      cout <<  setw(5) << a[i+j*NU] << setw(5);
                    }
                      cout <<  setw(10) << f[j+1] << endl;
                 }
             }
*/

