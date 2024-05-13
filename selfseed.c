#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>

double complex ***complexMemory3Asign(int harmony,int nz,int nx,int ny);
void complexDeleteField3(double complex ***field,int harmony,int subNz);


void washingOut(Domain *D,int iteration)
{
   int i,startI,endI,s,n,numInBeamlet;
	double dg,aveTh,minTh;
   Particle *particle;
   particle=D->particle;
   LoadList *LL;
   ptclList *p;

   int nTasks,myrank;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;  endI=D->subSliceN+1;

	LL=D->loadList;
   s=0;
   while(LL->next) {
      numInBeamlet=LL->numInBeamlet;
		dg=2*M_PI/(numInBeamlet*1.0);

      for(i=startI; i<endI; i++)
      {
	      p=particle[i].head[s]->pt;
		   while(p)  {
			   aveTh=0.0;
			   for(n=0; n<numInBeamlet; n++) aveTh+=p->theta[n];
				aveTh/=1.0*numInBeamlet;

				minTh=aveTh-M_PI;
			   for(n=0; n<numInBeamlet; n++) p->theta[n]=minTh+n*dg;

				p=p->next;
			}
		}

      LL=LL->next;
		s++;
	}
}

void selfSeed_Field(Domain *D,int iteration)
{
   int h,i,j,n,rank,N,startI,endI,numHarmony,minI,maxI,nn,ii;
	int dataNum,start,sliceN,subN,*minmax,numSlice,indexI;
	double delayT,dt,val,tmp,tau,ctau,arg,coef,sinTh;
	double k0,shiftT,d,extincL,*sendData,*recvData,realV,imagV,*U;
	double complex chi0,tmpComp,compVal,J,first,result,*data;
   int myrank, nTasks;
	FILE *out;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
	minI=D->minI; maxI=D->maxI;
	sliceN=D->sliceN;
   N=D->nx*D->ny;

   dt = D->lambda0*D->numLambdaU/velocityC;
	delayT = D->chi_delay;
	sinTh = sin(D->bragTh);
	d=D->chi_d;
	extincL=D->extincL;
	chi0=D->chi0;
	k0=D->ks;
   coef=M_PI*M_PI*d*sinTh/(extincL*extincL);	


	dataNum=D->numHarmony*N*(sliceN+2)*2;
   sendData=(double *)malloc(dataNum*sizeof(double ));
   recvData=(double *)malloc(dataNum*sizeof(double ));
   U=(double *)malloc(dataNum*sizeof(double ));
	for(i=0; i<dataNum; i++) U[i]=0.0;
   numSlice=sliceN+2;
   //U=complexMemory3Asign(D->numHarmony,sliceN+2,D->nx,D->ny);

   //double maxT,dtt;
	//int numT;
   //maxT=100e-15;
	//numT=200;
	//dtt=maxT/(numT*1.0);
	//chi0=-2.12e-5;
   //sinTh=sin(67.2966*M_PI/180.0);
	//extincL=22.8e-6;
	//d=100e-6;
   //coef=M_PI*M_PI*d*sinTh/(extincL*extincL);
   data=(double complex *)malloc(numSlice*sizeof(double complex));
	for(i=0; i<numSlice; i++) data[i]=0.0+I*0.0;

   for(h=0; h<1; h++) 
     for(i=startI; i<endI; i++) 
	  {
       shiftT=(sliceN - (minI+i-startI))*dt + delayT;

       for(j=0; j<N; j++) 
		 {
         compVal=D->U[h][i][j];
		  
         for(nn=1; nn<numSlice; nn++) {
           tau = shiftT-(nn-1)*dt;
			  if(tau>=0) {
             ctau=velocityC*tau;

             tmp=ctau*(2.0*d/sinTh+ctau/sinTh/sinTh);
	          if(tmp==0.0) J=0.5;
	          else if(tmp<0) {
		         arg=M_PI/extincL*sqrt(fabs(tmp));
		         J=gsl_sf_bessel_I1(arg); J/=arg;
				   J*=-I;
			    } else {
			      arg=M_PI/extincL*sqrt(tmp);
			      J=gsl_sf_bessel_J1(arg); J/=arg;
			    }

             tmpComp=chi0*k0*(d+ctau/sinTh)/2.0/sinTh;
             first=cexp(I*tmpComp);
             result=coef*first*J*compVal;
//				 data[nn]+=result*dt*velocityC;
             U[h*numSlice*N*2 + (nn-1)*(N*2) + j*2 + 0]+=creal(result)*dt*velocityC;
             U[h*numSlice*N*2 + (nn-1)*(N*2) + j*2 + 1]+=cimag(result)*dt*velocityC;
			  } else ;
         }
		 }
	  }
	 

   for(i=0; i<dataNum; i++) sendData[i]=U[i];

   if(myrank==0)  {
     for(rank=1; rank<nTasks; rank++) {
       MPI_Recv(recvData,dataNum,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
		 for(i=0; i<dataNum; i++)
		   U[i]+=recvData[i];
     }
	} else {
     for(rank=1; rank<nTasks; rank++) {
	    if(myrank==rank) 
         MPI_Send(sendData,dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
	}
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(U,dataNum,MPI_DOUBLE,0,MPI_COMM_WORLD);

/*
   out=fopen("conField","w");
	for(n=1; n<sliceN+1; n++) {
	  tau=(n-1)*dt;
	  fprintf(out,"%g %g\n",tau,cabs(data[n]));
	}
	fclose(out);
	printf("'conField' is made.\n");
*/
   for(h=0; h<numHarmony; h++) 
     for(i=startI; i<endI; i++) {
       nn=minI+i;
       for(j=0; j<N; j++) { 
         realV=U[h*numSlice*N*2 + nn*(N*2) + j*2 + 0];
         imagV=U[h*numSlice*N*2 + nn*(N*2) + j*2 + 1];
         D->U[h][i][j]=realV+I*imagV;
  //       D->U[h][i][j]=data[nn];
		 }
     }

   free(U);
	free(sendData);
	free(recvData);
	free(data);
   //complexDeleteField3(U,D->numHarmony,sliceN+2);
}


void seed_Field_test(Domain *D,int iteration)
{
   int nn;
	double delayT,dt,tmp,tau,ctau,arg,coef,sinTh;
	double k0,shiftT,d,extincL,maxT;
	double complex chi0,tmpComp,compVal,first,result,*U,*listJ,J;
   FILE *out;
   char fileName[100];
   int myrank, nTasks;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   dt = D->lambda0*D->numLambdaU/velocityC;
	delayT = D->chi_delay;
	sinTh = sin(D->bragTh);
	d=D->chi_d;
	extincL=D->extincL;
	chi0=D->chi0;
	k0=D->ks;
   coef=M_PI*M_PI*D->chi_d*sinTh/(D->extincL*D->extincL);	
   //printf("bragg=%g, extincL=%g\n",D->bragTh,extincL);

   U=(double complex *)malloc(200*sizeof(double complex));
   listJ=(double complex *)malloc(200*sizeof(double complex));
   maxT=100e-15;
   dt=maxT/200.0;

   for(nn=0; nn<200; nn++) {
      tau = nn*dt;
      ctau=velocityC*tau;

      tmp=ctau*(2.0*d/sinTh+ctau/sinTh/sinTh);
		printf("nn=%d, tmp=%g, d=%g, sinTh=%g\n",nn,tmp,d,sinTh);
      if(tmp==0.0) J=0.5;
      else if(tmp<0) {
        arg=M_PI/extincL*sqrt(fabs(tmp));
        J=gsl_sf_bessel_I1(arg); J/=arg;
        J*=-I;
		} else {
		  arg=M_PI/extincL*sqrt(tmp);
		  J=gsl_sf_bessel_J1(arg); J/=arg;
		}

      listJ[nn]=J;

      tmpComp=chi0*k0*(d+ctau/sinTh)/2.0/sinTh;
      first=cexp(I*tmpComp);
      result=coef*first*J;
      
      U[nn]=result;
   }

   if(myrank==0) {
      sprintf(fileName,"sswake");
      out=fopen(fileName,"w");
      for(nn=0; nn<200; nn++) {
        tau = nn*dt;
        fprintf(out,"%g %g %g\n",tau,cabs(U[nn]),cabs(listJ[nn]));
      }
      fclose(out);
      printf("%s is made.\n",fileName);
    }

   free(U);
   free(listJ);
}


void whatCrystal(double ks,ChiList *Chi,char *str)
{
   double d,chi0R,chi0I,energy;
	double a0,a1,a2,a3,a4,a5,a6;
	double b0,b1,b2,b3,b4,b5,b6;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   energy=ks*velocityC*hbar;   //eV
	if(energy<4000) {
	   if(myrank==0) 
		   printf("photon energy is %g eV, which is not allowed for now!\n",energy); 
		else ;
		exit(0);
	} else ;

   if(strstr(str,"Si_100"))  {
     d=0.543096892e-9;     //grating constant

	  a0=-0.291154;
	  a1=-1037.78;
	  b0=4.90666e-5;
	  b1=-3.54453e-9;
	  b2=9.15755e-14;
	  chi0R=a0/(energy+a1) + b0+b1*energy+b2*pow(energy,2);
	  a0=-0.0298072;
	  a1=-1799.71;
	  b0=1.61241e-5;
	  b1=-2.53438e-9;
	  b2=1.65622e-13;
	  b3=-3.88831e-18;
	  chi0I=a0/(energy+a1) + b0+b1*energy+b2*pow(energy,2)+b3*pow(energy,3);
	  Chi->chi0=chi0R+I*chi0I;
	  a0=5.03719e7;
	  a1=497.096;
	  a2=-0.0225941;
	  a3=2.2854e-8;
	  a4=1.63355e-11;
     Chi->extincL=(a0+a1*energy+a2*pow(energy,2)+a3*pow(energy,3)+a4*pow(energy,4))*1e-6;
	  Chi->bragTh=asin(M_PI/(d*ks));
   }	
   else   {
     printf("No crystall! Define the crystall\n"); 
     exit(0);
   }
}




