EXEC = show
CC = /opt/ompi/3.1.6/bin/mpicc
#CC = mpicc
OBJS = main.o parameterSetting.o findparam.o boundary.o loadBeam.o updateTotalEnergy.o saveFile.o saveParticleHDF.o saveFieldHDF.o particlePush.o rearrangeParticles.o updateK_quadG.o solveField.o fieldShareZ.o chicane.o 

#OBJS = loadSeed.o clean.o shiftField.o selfseed.o 


#-------- for Beam server ----------# 
CFLAGS = -I/opt/gsl/2.6/include -I/opt/hdf5/1.10.3/include 
LDFLAGS = -L/opt/gsl/2.6/lib -L/opt/hdf5/1.10.3/lib

#-------- for PAL ----------#
#CFLAGS = -I/home/scienter/gsl2.6/include
#LDFLAGS = -L/home/scienter/gsl2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi

#-------- for home computer ----------#
#CFLAGS= -I/usr/include/hdf5/openmpi
#LDFLAGS= -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi


INCL = constants.h mesh.h particle.h
LIBS = -lm -lgsl -lgslcblas -lhdf5 -lm -lz

$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
#	$(CC)           $(OBJS)            $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
