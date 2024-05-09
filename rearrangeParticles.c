#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

/*
void rearrangeParticles(Domain *D,int iteration)
{
    Particle *particle;
    particle=D->particle;

    int i,s,intZ,cnt,deleteFlag=0;
    int startI,endI,nSpecies;
    double dPhi,theta;
    ptclList *p,*New,*prev,*tmp;

    startI=1;  endI=D->subSliceN+1;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

    for(i=startI; i<endI; i++)
    {
      for(s=0; s<nSpecies; s++) {
        cnt=1;
        p=particle[i].head[s]->pt;
        while(p)  {
          if(cnt==1)
            prev=p;
          deleteFlag=0;
              
          theta=p->theta;
          if(theta>=dPhi)  {
            intZ=1;
	    theta-=dPhi;
            deleteFlag=1;
          }
          else if(theta<0) {              
            intZ=-1;
	    theta+=dPhi;
            deleteFlag=1;
          } 
          else   intZ=0;

          if(deleteFlag==1)  {
            if(cnt==1)  {
              p->theta=theta;    
              particle[i].head[s]->pt = p->next;
              p->next = particle[i+intZ].head[s]->pt;
              particle[i+intZ].head[s]->pt = p;
              p=particle[i].head[s]->pt;
              cnt=1;
            } else {
              prev->next = p->next;
              p->theta=theta;    
              p->next = particle[i+intZ].head[s]->pt;
              particle[i+intZ].head[s]->pt = p;
              p=prev->next;
            }
          }		//End of if(deleteFlag==1)
          else {
            prev=p;
            p=p->next;
            cnt++;
          }              
        }	//End of while(p)
      }		//End of for(s)
    }		//End of for(i)
}
*/

void periodicParticles(Domain *D,int iteration)
{
    Particle *particle;
    particle=D->particle;

    int i,s,intThe,n,*N,calFlag;
    int startI,endI,nSpecies;
    double dPhi,aveTh,delTh;
    LoadList *LL;    
    ptclList *p;

    startI=1;  endI=D->subSliceN+1;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

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
       for(s=0; s<nSpecies; s++) {
          p=particle[i].head[s]->pt;
          while(p)  {
             calFlag=OFF;
             aveTh=0.0;
             for(n=0; n<N[s]; n++) aveTh+=p->theta[n];
             aveTh/=1.0*N[s];

             if(aveTh>=dPhi)  {
                intThe=(int)(aveTh/dPhi);
                delTh=dPhi*intThe;
                calFlag=ON;
	             //theta-=dPhi;
             } else if(aveTh<0) {              
                intThe=(int)(aveTh/dPhi);
                delTh=dPhi*(intThe-1);
                calFlag=ON;
	             //theta+=dPhi;
             } else ;

             if(calFlag==ON) {
                for(n=0; n<N[s]; n++) p->theta[n]-=delTh;    
                aveTh-=delTh;
             } else ;

             if(aveTh>=dPhi || aveTh<0) { printf("In rearrange, intThe%d,delTh=%g,newAve=%g,iteration=%d,dPhi=%g\n",intThe,delTh,aveTh,iteration,dPhi);  }
          
	          p=p->next;
	       }
       }
    }		//End of for(i)

    free(N);
}

