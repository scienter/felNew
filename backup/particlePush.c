#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>

void transversePush_3D(Domain *D,int iteration);
void push_theta_gamma_3D(Domain *D,int iteration);
void push_theta_gamma_1D(Domain *D,int iteration);
void drift_theta_gamma_3D(Domain *D,int iteration);
void drift_theta_gamma_1D(Domain *D,int iteration);

void transversePush(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
	  ;
//    particlePush1D(D);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
//    transversePush_3D(D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}


void push_theta_gamma(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    push_theta_gamma_1D(D,iteration);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
//    push_theta_gamma_3D(D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}

void drift_theta_gamma(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    drift_theta_gamma_1D(D,iteration);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
//    drift_theta_gamma_3D(D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}

/*
void transversePush_3D(Domain *D,int iteration)
{
    int i,j,s,sliceI,startI,endI;
    LoadList *LL;
    double dz,ks,kx,ky,K0,g;
    double x,y,z,px,py,gamma,invGam,v0[2],v1[2],M[2][2];
    double qx,qy,sqrtQx,sqrtQy,cnt,avePx,avePy,aveGam,send[4],recv[4];
    ptclList *p;

    int myrank, nTasks;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    dz=D->dz;
    ks=D->ks;
    kx=D->kx;    ky=D->ky;
    K0=D->K0;
    g=D->g;
    startI=1;  endI=D->subSliceN+1;

    avePx=avePy=aveGam=cnt=0.0;
    for(sliceI=startI; sliceI<endI; sliceI++)
    {    
      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[sliceI].head[s]->pt;
        while(p) {
       	  gamma=p->gamma;  invGam=1.0/gamma;
	
          qx=(K0*K0*invGam*kx*kx+eCharge*g/eMass/velocityC)*invGam;
          qy=(K0*K0*invGam*ky*ky-eCharge*g/eMass/velocityC)*invGam;
          sqrtQx=sqrt(fabs(qx));
          sqrtQy=sqrt(fabs(qy));

	  // calculate for x direction
  	  if(qx>0) {
            M[0][0]=cos(sqrtQx*dz*0.5);
            M[0][1]=sin(sqrtQx*dz*0.5)/sqrtQx*invGam;
            M[1][0]=-gamma*sqrtQx*sin(sqrtQx*dz*0.5);
            M[1][1]=M[0][0];
  	  } else if (qx<0) {
            M[0][0]=cosh(sqrtQx*dz*0.5);
            M[0][1]=sinh(sqrtQx*dz*0.5)/sqrtQx*invGam;
            M[1][0]=gamma*sqrtQx*sinh(sqrtQx*dz*0.5);
            M[1][1]=M[0][0];
  	  } else {
            M[0][0]=1.0;
            M[0][1]=dz*0.5*invGam;
            M[1][0]=0.0;
            M[1][1]=1.0;
  	  }
          v0[0]=p->x; v0[1]=p->px;
  	  v1[0]=0.0;  v1[1]=0.0;
	  for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
              v1[i]+=M[i][j]*v0[j];
          p->x=v1[0];
          p->px=v1[1];
		
          // calculate for y direction
  	  if(qy>0) {
            M[0][0]=cos(sqrtQy*dz*0.5);
            M[0][1]=sin(sqrtQy*dz*0.5)/sqrtQy*invGam;
            M[1][0]=-gamma*sqrtQy*sin(sqrtQy*dz*0.5);
            M[1][1]=M[0][0];
  	  } else if (qy<0) {
            M[0][0]=cosh(sqrtQy*dz*0.5);
            M[0][1]=sinh(sqrtQy*dz*0.5)/sqrtQy*invGam;
            M[1][0]=gamma*sqrtQy*sinh(sqrtQy*dz*0.5);
            M[1][1]=M[0][0];
  	  } else {
            M[0][0]=1.0;
            M[0][1]=dz*0.5*invGam;
            M[1][0]=0.0;
            M[1][1]=1.0;
  	  }
          v0[0]=p->y; v0[1]=p->py;
  	  v1[0]=0.0;  v1[1]=0.0;
	  for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
              v1[i]+=M[i][j]*v0[j];
          p->y=v1[0];
          p->py=v1[1];

	  avePx+=p->px*p->weight;
	  avePy+=p->py*p->weight;
	  aveGam+=p->gamma*p->weight;
	  cnt+=p->weight;

	  p=p->next;
        }
      }		//End of for(s)
    }		//End of for(i)
         
}
*/
/*
void push_theta_gamma_3D(Domain *D,int iteration)
{
    int i,j,N,ii,jj,s,h,H,numHarmony,nx,ny,order,ll,L,f,F;
    int startI,endI,minI,maxI,sliceI,indexJ,idx;
    LoadList *LL;
    double complex U[D->numHarmony],compVal,Em[D->SCLmode];
    double dz,dx,dy,ku,ks,kx,ky,K0,K,xi,e_mc2,r,dr,dBessel;
    double x,y,z,px,py,gamma,theta,invGam,minX,minY,wakeE,invBeta,phi;
    double coef,JJ,J1,J2,wx[2],wy[2],sumTh,sumGam,sumEzPart,w[2];
    double k1,k2,k3,k4,l1,l2,l3,l4,tmp,dPhi,amp,arg,chi;
    ptclList *p;

    dz=D->dz;    K0=D->K0;
    ku=D->ku;    ks=D->ks;
    numHarmony=D->numHarmony;
    kx=D->kx; ky=D->ky;
    dx=D->dx; dy=D->dy; dr=D->dr;
    nx=D->nx; ny=D->ny;
	 dBessel = D->dBessel;
    minX=D->minX;  minY=D->minY;
    dPhi=2*M_PI*D->numSlice;
    e_mc2 = eCharge/eMass/velocityC/velocityC;	 
    
    L = D->SCLmode; F = D->SCFmode;	 
    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;
    N=D->nx*D->ny;

    for(sliceI=startI; sliceI<endI; sliceI++)
    {
      if(D->wakeONOFF==ON) wakeE=D->wakeE[sliceI-startI+minI]/mc2*1e-6;
      else                wakeE=0.0;      

      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[sliceI].head[s]->pt;
        while(p) {
          x=p->x;    y=p->y;  r=sqrt(x*x+y*y);
          px=p->px;  py=p->py;
		    if(x==0) phi = 0;
          else     phi = atan2(y,x);


  	       i=(int)((x-minX)/dx);
	       j=(int)((y-minY)/dy);
	       if(i>=0 && i<nx-1 && j>=0 && j<ny-1)  {
			   indexJ = (int)(r/dr);
				wy[1]=(r/dr-indexJ); wy[0]=1.0-wy[1];
     			for(ll=0; ll<L; ll++) {
				  Em[ll]=0.0+I*0.0;
              for(f=0; f<F; f++) {
//					 amp=0.0; arg=0.0;
				    for(jj=0; jj<2; jj++) {
//				      amp+=cabs(D->Ez[sliceI][indexJ][ll][f])*wy[jj];
//				      arg+=carg(D->Ez[sliceI][indexJ][ll][f])*wy[jj];
				      Em[ll]+=D->Ez[sliceI][indexJ][ll][f]*wy[jj];
					 }
//				    Em[ll] += amp*cexp(I*(f*phi+arg));
				  }
				}


       	   wx[1]=(x-minX)/dx-i; wx[0]=1.0-wx[1];
	         wy[1]=(y-minY)/dy-j; wy[0]=1.0-wy[1];
	         for(h=0; h<numHarmony; h++)  {
				
				  U[h]=0.0+I*0.0;
	           for(ii=0; ii<2; ii++) 
                for(jj=0; jj<2; jj++)  
      			   U[h]+=D->U[h][sliceI][(j+jj)*nx+(i+ii)]*wx[ii]*wy[jj];
					
//      		  U[h]=D->U[h][sliceI][(j)*nx+(i)];
	         }

	         K=K0*(1.0+kx*kx*x*x+ky*ky*y*y);
            
				l2=l3=l4=k2=k3=k4=0.0;
				// Step 1
	         theta=p->theta;
	         gamma=p->gamma; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
            sumTh=sumGam=0.0;
	         //xi=ks/ku*0.25*K*K*invGam*invGam;
	         xi=K*K*0.5/(1+K*K);
            for(h=0; h<numHarmony; h++)  {
				  H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                //chi = H*sqrt(2.0)*K*ks/ku*px*invGam*invGam;
                chi = H*sqrt(2.0)*K*px*2/(1+K*K);
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k1=dz*(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l1=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);
            
				// Step 2 2 
	         theta=p->theta+0.5*k1;
	         gamma=p->gamma+0.5*l1; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
	         xi=K*K*0.5/(1+K*K);
            sumTh=sumGam=0.0;
            for(h=0; h<numHarmony; h++)  {
	           H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*px*2/(1+K*K);
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);					 
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k2=dz*(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l2=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);

				// Step 3 
	         theta=p->theta+0.5*k2;
	         gamma=p->gamma+0.5*l2; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
	         xi=K*K*0.5/(1+K*K);
            sumTh=sumGam=0.0;
            for(h=0; h<numHarmony; h++)  {
	           H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*px*2/(1+K*K);
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k3=dz*(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l3=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);

				// Step 4 
	         theta=p->theta+k3;
	         gamma=p->gamma+l3; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K + px*px + py*py)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;			 
	         xi=K*K*0.5/(1+K*K);
            sumTh=sumGam=0.0;
            for(h=0; h<numHarmony; h++)  {
	           H = D->harmony[h];
              if(H%2==1)  {  //odd harmony
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
              } else {
                chi = H*sqrt(2.0)*K*px*2/(1+K*K);
                coef=pow(-1.0,(H-2)*0.5);
                idx=(int)(H*xi/dBessel);
					 w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-2)*0.5;
					 J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+2)*0.5;
					 J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=chi*0.5*coef*(J1-J2);
              }
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k4=dz*(ku-ks*(1.0+K*K+px*px+py*py+sumTh)*0.5*invGam*invGam)*invBeta;
	         l4=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);

            tmp=1.0/6.0*(k1+2*k2+2*k3+k4);
	         if(tmp>=dPhi || tmp<=-dPhi) {
              printf("iteration=%d, dTheta=%g, r=%g,sumEzPart=%g\n",iteration,tmp,r,sumEzPart); 
              exit(0);
	         } else;
            p->theta+=tmp;
            p->gamma+=1.0/6.0*(l1+2*l2+2*l3+l4)-dz*wakeE;

	       } else ;

	       p=p->next;
        }
      }		//End of for(s)
    }		//End of for(sliceI)

}
*/


void push_theta_gamma_1D(Domain *D,int iteration)
{
    int n,i,s,h,H,numHarmony,order,startI,endI,minI,maxI,ll,L,idx,intThe,bn,*N;
    LoadList *LL;
    double complex U[D->numHarmony],Em[D->SCLmode],compVal;
    double dz,ku,ks,K0,K,xi,e_mc2,dBessel,w[2];
    double z,gamma,theta,invGam,invBeta,tmp,dPhi,sumGam,sumTh,sumEzPart,prevThe;
    double coef,JJ,J1,J2,wakeE,absU,absU2;
    double k1,k2,k3,k4,l1,l2,l3,l4;
	 double beta1,beta2,beta3,beta4;
	 double sumGam1,sumGam2,sumGam3,sumGam4;
    ptclList *p;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    L = D->SCLmode;
    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;
	 dz = D->dz;     K0=D->K0;
    ku=D->ku;       ks=D->ks;
    numHarmony = D->numHarmony;
    dBessel = D->dBessel;
    dPhi=2*M_PI*D->numSlice;
    e_mc2 = eCharge/eMass/velocityC/velocityC;	 
	 bn=D->bn;

    N=(int *)malloc(D->nSpecies*sizeof(int ));
    LL=D->loadList;
    s=0;
    while(LL->next) {
       N[s]=LL->numInBeamlet;
       LL=LL->next;
       s++;
    }

    for(i=startI; i<endI; i++)
    {
      if(D->wakeONOFF==ON) wakeE=D->wakeE[i-startI+minI]/mc2*1e-6;
      else                 wakeE=0.0;      

      for(ll=0; ll<L; ll++) 
		  Em[ll]=D->Ez[i][0][ll][0];

      for(h=0; h<numHarmony; h++)  
		  U[h]=D->U[h][i][0];
      
      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[i].head[s]->pt;

        while(p) {
          for(n=0; n<N[s]; n++) {
            K=K0;

            theta=p->theta[n];
            gamma=p->gamma[n]; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;
            sumTh=sumGam=0.0;
            //xi=ks/ku*0.25*K*K*invGam*invGam;
	         xi=K*K*0.5/(1+K*K);
            for(h=0; h<numHarmony; h++)  {
			     H = D->harmony[h];
			     if(H%2==1) {
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
                if(idx>bn-1) { printf("idx=%d\n",idx); idx=bn-2; }
                w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
                J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
                J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
				  } else JJ=0.0;
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k1=dz*(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;
            l1=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);
			 
	         theta=p->theta[n]+0.5*k1;
            gamma=p->gamma[n]+0.5*l1; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;
            sumTh=sumGam=0.0;
	         xi=K*K*0.5/(1+K*K);
            for(h=0; h<numHarmony; h++)  {
			     H = D->harmony[h];
			     if(H%2==1) {
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
				    if(idx>bn-1) { printf("idx=%d\n",idx); idx=bn-2; }
                w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
                J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
                J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
				  } else JJ=0.0;
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k2=dz*(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;
            l2=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);

	         theta=p->theta[n]+0.5*k2;
            gamma=p->gamma[n]+0.5*l2; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;
            sumTh=sumGam=0.0;
	         xi=K*K*0.5/(1+K*K);
            for(h=0; h<numHarmony; h++)  {
			     H = D->harmony[h];
			     if(H%2==1) {
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
				    if(idx>bn-1) { printf("idx=%d\n",idx); idx=bn-2; }
                w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
                J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
                J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
				  } else JJ=0.0;
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k3=dz*(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;
            l3=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);

	         theta=p->theta[n]+k3;
            gamma=p->gamma[n]+l3; invGam=1.0/gamma;
            invBeta = 1.0-(1.0 + K*K)*0.5*invGam*invGam;
            invBeta = 1.0/invBeta;
            sumTh=sumGam=0.0;
	         xi=K*K*0.5/(1+K*K);
            for(h=0; h<numHarmony; h++)  {
			     H = D->harmony[h];
			     if(H%2==1) {
                coef=pow(-1.0,(H-1)*0.5);
                idx=(int)(H*xi/dBessel);
				    if(idx>bn-1) { printf("idx=%d\n",idx); idx=bn-2; }
                w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                order=(H-1)*0.5;
                J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                order=(H+1)*0.5;
                J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                JJ=coef*(J1-J2);
				  } else JJ=0.0;
              compVal=U[h]*cexp(I*H*theta);
              sumTh -=2*JJ*K*cimag(compVal);
              sumGam-=2*JJ*K*creal(compVal);
            }
            sumEzPart = 0.0;
            for(ll=0; ll<L; ll++)  {
              tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
              sumEzPart += 2.0*tmp;
            }
            k4=dz*(ku-ks*(1.0+K*K+sumTh)*0.5*invGam*invGam)*invBeta;
            l4=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);

            tmp=dz/6.0*(k1+2*k2+2*k3+k4);
	         if(tmp>dPhi || tmp<-dPhi) { 
              printf("iteration=%d, dTheta=%g, sumEzPart=%g, U[0]=%g\n",iteration,tmp,sumEzPart,cabs(D->U[0][i][0]));  //lala
              exit(0);
            }
            p->theta[n]-=tmp;
            p->gamma[n]+=dz/6.0*(l1+2*l2+2*l3+l4) - wakeE*dz;
          }   //End of for(n)

	       p=p->next;
        }
      }		//End of for(s)
    }		//End of for(i)

   free(N);
}


/*
void drift_theta_gamma_3D(Domain *D,int iteration)
{
   int i,j,N,ii,jj,s,h,H,numHarmony,nx,ny,order,ll,L,f,F;
   int startI,endI,minI,maxI,sliceI,indexJ,idx;
   LoadList *LL;
   double complex U[D->numHarmony];
   double dz,dx,dy,ks,kx,ky,K0,K,xi,e_mc2,r,dr;
   double x,y,z,px,py,gamma,invGam,minX,minY,wakeE;
   double coef,JJ,J1,J2,wx[2],wy[2],w[2],dPhi;
   ptclList *p;

   dz=D->dz;    K0=D->K0;   ks=D->ks;
   numHarmony=D->numHarmony;
   dx=D->dx; dy=D->dy;
   nx=D->nx; ny=D->ny;
   kx=D->kx; ky=D->ky;
   minX=D->minX;  minY=D->minY;
   dPhi=2*M_PI*D->numSlice;
   e_mc2 = eCharge/eMass/velocityC/velocityC;	 
    
   startI=1;       endI=D->subSliceN+1;
   minI=D->minI;   maxI=D->maxI;
   
   for(sliceI=startI; sliceI<endI; sliceI++)
   {
     if(D->wakeONOFF==ON) wakeE=D->wakeE[sliceI-startI+minI]/mc2*1e-6;
     else                wakeE=0.0;      

     for(s=0; s<D->nSpecies; s++)  {
       p=D->particle[sliceI].head[s]->pt;
       while(p) {
         x=p->x;    y=p->y;  r=sqrt(x*x+y*y);
         px=p->px;  py=p->py;

	      K=K0*(1.0+kx*kx*x*x+ky*ky*y*y);

  	      i=(int)((x-minX)/dx);
	      j=(int)((y-minY)/dy);
	      if(i>=0 && i<nx-1 && j>=0 && j<ny-1)  {

       	  wx[1]=(x-minX)/dx-i; wx[0]=1.0-wx[1];
	        wy[1]=(y-minY)/dy-j; wy[0]=1.0-wy[1];
	        for(h=0; h<numHarmony; h++)  {				
			    U[h]=0.0+I*0.0;
	          for(ii=0; ii<2; ii++) 
               for(jj=0; jj<2; jj++)  
      	        U[h]+=D->U[h][sliceI][(j+jj)*nx+(i+ii)]*wx[ii]*wy[jj];					
	        }

	        gamma=p->gamma; invGam=1.0/gamma;
			  p->theta-=ks*0.5*invGam*invGam*(1+px*px+py*py)*dz;
           p->gamma-=wakeE*dz;

	      } else ;

	      p=p->next;
       }
     }		//End of for(s)
   }		//End of for(sliceI)

}
*/

void drift_theta_gamma_1D(Domain *D,int iteration)
{
   int n,i,s,h,H,numHarmony,startI,endI,minI,maxI,*N;
   LoadList *LL;
   double complex U[D->numHarmony];
   double dz,ks,K0,K,e_mc2,gamma,invGam,dPhi,px,py,wakeE,invGam0,invBeta0,tmp;
   ptclList *p;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;       endI=D->subSliceN+1;
   minI=D->minI;   maxI=D->maxI;
	dz = D->dz;     K0=D->K0;
   ks=D->ks;
   numHarmony = D->numHarmony;
   dPhi=2*M_PI*D->numSlice;
   e_mc2 = eCharge/eMass/velocityC/velocityC;	 
	invGam0=1.0/D->gamma0;
	invBeta0=1.0/sqrt(1-invGam0*invGam0);

   N=(int *)malloc(D->nSpecies*sizeof(int ));
   LL=D->loadList;
   s=0;
   while(LL->next) {
      N[s]=LL->numInBeamlet;
      LL=LL->next;
      s++;
   }

   if(myrank==0) printf("iteration=%d, drift ON\n",iteration);

   for(i=startI; i<endI; i++)
   {
     if(D->wakeONOFF==ON) wakeE=D->wakeE[i-startI+minI]/mc2*1e-6;
     else                 wakeE=0.0;      

     for(h=0; h<numHarmony; h++) U[h]=D->U[h][i][0];
      
     for(s=0; s<D->nSpecies; s++)  
	  {
       p=D->particle[i].head[s]->pt;

       while(p) {
         for(n=0; n<N[s]; n++) {
           px=p->px[n];  py=p->py[n];

           K=K0;
           gamma=p->gamma[n]; invGam=1.0/gamma;
           tmp = ks*dz*0.5*invBeta0*(invGam0*invGam0-invGam*invGam);
			  p->theta[n]+=tmp;
           p->gamma[n]-=wakeE*dz;
         }

         p=p->next;
       }
     }		//End of for(s)
   }		//End of for(i)

   free(N);
}

	
double Runge_Kutta_gamma(double complex *U,double theta,double harmony,double ks,double K,double xi,double gamma,double invBeta0,double dz)
{
  double sum,sinX,cosX,coef,J1,J2,JJ,k1,k2,k3,k4;
  double complex compVal;
  int h,order;

  sum=0.0;
  sinX=sin(theta); cosX=cos(theta);
  for(h=1; h<=harmony; h++)  {
    if(h%2==1)  {  //odd harmony
      coef=pow(-1.0,h-1);
      order=(h-1)/2;
      J1=gsl_sf_bessel_Jn(order,h*xi);
      order=(h+1)/2;
      J2=gsl_sf_bessel_Jn(order,h*xi);
      JJ=coef*(J1-J2);
      compVal=U[h-1]*cexp(-I*theta);
      sum=JJ*2*creal(compVal)*h;
    } else {
      JJ=0.0;
      compVal=U[h-1]*cexp(-I*theta);
      sum=JJ*2*creal(compVal)*h;
    }
  }

  k1=-1.0*K*ks*invBeta0*sum/(gamma);
  
  k2=-1.0*K*ks*invBeta0*sum/(gamma+0.5*k1);
  
  k3=-1.0*K*ks*invBeta0*sum/(gamma+0.5*k2);
  
  k4=-1.0*K*ks*invBeta0*sum/(gamma+k3);

  return (k1/6.0+k2/3.0+k3/3.0+k4/6.0)*dz;
}

/*
void phaseShift(Domain *D,int iteration)
{
	 int s,i,sliceI,startI,endI;
	 double shiftValue,theta;
    LoadList *LL;
    ptclList *p;
	 PhaseShifter *PS;
    int myrank, nTasks;
    MPI_Status status;
	 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	 MPI_Comm_size(MPI_COMM_WORLD, &nTasks); 

    startI=1;       endI=D->subSliceN+1;

    PS=D->psList;
    while(PS->next) {
      for(i=0; i<PS->num; i++) {
        if(iteration==PS->step[i]) {
          shiftValue=PS->phase;
          if(myrank==0) printf("phase shift with %g is done at step%d.\n",shiftValue,iteration);  else ;

          for(sliceI=startI; sliceI<endI; sliceI++)  {
            for(s=0; s<D->nSpecies; s++)  {
              p=D->particle[sliceI].head[s]->pt;
              while(p) {
			       theta=p->theta;
                p->theta=theta-shiftValue;

	             p=p->next;
              }
            }		//End of for(s)
          }		//End of for(sliceI)

		  } else ;
	   }
		PS=PS->next;
	 }

}
*/
