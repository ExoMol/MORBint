goal:   morbint

tarball:
	tar cf renner.tar makefile *.f
        
checkin:
	ci -l Makefile *.f

FOR = ifort 
FFLAGS = -ip   -openmp -O3  -static


LAPACK = -mkl=parallel

LIB     =   $(LAPACK) 

%.o : %.f
	$(FOR) -c $(FFLAGS) $<

OBJ    = yfunct.o strtch.o mutens.o mues.o  movibr.o morrea.o mormai.o morbid.o  \
         moparm.o mocori.o intone.o genstr.o diagro.o delks.o  clock.o

morbint:	$(OBJ)
	$(FOR) -o morbint.x $(FFLAGS) $(OBJ) $(LIB)
clean:
	rm $(OBJ)

yfunct.o:	yfunct.f;	$(FOR) -c $(FFLAGS) yfunct.f
strtch.o:	strtch.f;	$(FOR) -c $(FFLAGS) strtch.f
mutens.o:	mutens.f;	$(FOR) -c $(FFLAGS) mutens.f
mues.o  :	mues.f  ;	$(FOR) -c $(FFLAGS) mues.f  
movibr.o:	movibr.f;	$(FOR) -c $(FFLAGS) movibr.f
morrea.o:	morrea.f;	$(FOR) -c $(FFLAGS) morrea.f
mormai.o:	mormai.f;	$(FOR) -c $(FFLAGS) mormai.f
morbid.o:	morbid.f;	$(FOR) -c $(FFLAGS) morbid.f
moparm.o:	moparm.f;	$(FOR) -c $(FFLAGS) moparm.f
mocori.o:	mocori.f;	$(FOR) -c $(FFLAGS) mocori.f
intone.o:	intone.f;	$(FOR) -c $(FFLAGS) intone.f
genstr.o:	genstr.f;	$(FOR) -c $(FFLAGS) genstr.f
diagro.o:	diagro.f;	$(FOR) -c $(FFLAGS) diagro.f
delks.o :	delks.f ;	$(FOR) -c $(FFLAGS) delks.f 
clock.o :	clock.f ;	$(FOR) -c $(FFLAGS) clock.f 
