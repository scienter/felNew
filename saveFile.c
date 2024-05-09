#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include <gsl/gsl_sf_bessel.h>

void saveTotalBFactor(Domain *D)
{
   int maxStep,step;
	double dz,minZ,z;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status; 

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   maxStep=D->maxStep;
   dz=D->dz;
   minZ=D->minZ;

   if(myrank==0) {
      out=fopen("totalB","w");
      for(step=0; step<maxStep; step++) {
         z=step*dz+minZ;
         fprintf(out,"%.15g %g",z,D->totalBunch[step]);
      }
      fclose(out);
      printf("totalB is made.\n");
	} else ;
}




void updateBFactor(Domain *D,int iteration)
{
   int i,j,s,startI,endI,minI,sliceI,sliceN,rank;
	double bucketZ,z,theta,send[3],recv[3],cnt;
	double complex sum;
   LoadList *LL;
	ptclList *p;
   char name[100];
   FILE *out;

   int myrank, nTasks;    
   MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;  endI=D->subSliceN+1;
	minI=D->minI;

	sliceN=D->sliceN;

   sum=0.0+I*0.0;
	cnt=0.0;
   for(i=startI; i<endI; i++)
   {
     for(s=0; s<D->nSpecies; s++)  {
       p=D->particle[i].head[s]->pt;
       while(p) {
         theta=p->theta;
			sum+=cexp(I*theta);
			cnt+=1.0;
			p=p->next;
		 }
	  }
	}
   send[0]=creal(sum);
   send[1]=cimag(sum);
   send[2]=cnt;

	recv[0]=0.0;
	recv[1]=0.0;
	recv[2]=0.0;

   if(myrank==0)  {
     for(rank=1; rank<nTasks; rank++) {
	    MPI_Recv(recv,3,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
       sum+=recv[0]+I*recv[1];
       cnt+=recv[2];
     }
	} else {
	  for(rank=1; rank<nTasks; rank++) {
	    if(myrank==rank)
		   MPI_Send(send,3,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
   }

   D->totalBunch[iteration]=cabs(sum)/cnt;

}


/*
void saveBFactor(Domain *D,int iteration)
{
   int i,j,s,startI,endI,minI,sliceI,sliceN,rank,cnt;
	double bucketZ,z,theta,*data,*recvData;
	double complex sum;
   LoadList *LL;
	ptclList *p;
   char name[100];
   FILE *out;

   int myrank, nTasks;    
   MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;  endI=D->subSliceN+1;
	minI=D->minI;

	sliceN=D->sliceN;
	recvData=(double *)malloc((sliceN+2)*sizeof(double ));
	data=(double *)malloc((sliceN+2)*sizeof(double ));
	for(i=0; i<(sliceN+2); i++) data[i]=0.0;

	
   for(i=startI; i<endI; i++)
   {
	  sliceI=i+minI;
	  cnt=0;
	  sum=0.0+I*0.0;
     for(s=0; s<D->nSpecies; s++)  {
       p=D->particle[i].head[s]->pt;
       while(p) {
         theta=p->theta;
			sum+=cexp(I*theta);
			cnt++;
			p=p->next;
		 }
	  }
     if(cnt>0) data[sliceI]=cabs(sum)/(cnt*1.0);
	}

   if(myrank==0)  {
     for(rank=1; rank<nTasks; rank++) {
	    MPI_Recv(recvData,sliceN+2,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
       for(i=0; i<sliceN+2; i++)
         data[i]+=recvData[i];
     }
	} else {
	  for(rank=1; rank<nTasks; rank++) {
	    if(myrank==rank)
		   MPI_Send(data,sliceN+2,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
   }

   bucketZ=D->numSlice*D->lambda0;
   if(myrank==0) {
     sprintf(name,"bFact%d",iteration);
     out = fopen(name,"w");

     for(i=1; i<sliceN+1; i++) {
       z=(i-1)*bucketZ;
       fprintf(out,"%g %g\n",z,data[i]);
     }
	  fclose(out);
     printf("%s is made.\n",name);
   } else ;

}
*/
