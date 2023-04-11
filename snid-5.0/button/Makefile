#
FC= f77
FSOURCE = buttmis.f button.f buttqbr.f buttqcf.f buttqch.f buttqex.f \
          buttqit.f buttqpr.f buttqxb.f buttqyb.f buttqytext.f buttsbr.f \
          buttscf.f buttsch.f buttsex.f buttsit.f buttspr.f buttsxb.f \
          buttsyb.f buttsytext.f ifbutton.f rpgband.f rpgbegin.f rpgbegok.f \
          rpgenv.f rpgeras.f rpgerasb.f rpgerasw.f
OBJECTS = buttmis.o button.o buttqbr.o buttqcf.o buttqch.o buttqex.o \
          buttqit.o buttqpr.o buttqxb.o buttqyb.o buttqytext.o buttsbr.o \
          buttscf.o buttsch.o buttsex.o buttsit.o buttspr.o buttsxb.o \
          buttsyb.o buttsytext.o ifbutton.o rpgband.o rpgbegin.o rpgbegok.o \
          rpgenv.o rpgeras.o rpgerasb.o rpgerasw.o
# Default rule to make object modules
.f.o:
	$(FC) -c $*.f
# Generate Library
libbutton.a: $(OBJECTS)
	ar cr libbutton.a $(OBJECTS)
	ranlib libbutton.a
# Target line to make object modules
buttmis.o: buttmis.f button.inc
button.o: button.f button.inc
buttqbr.o: buttqbr.f button.inc
buttqcf.o: buttqcf.f button.inc
buttqch.o: buttqch.f button.inc
buttqex.o: buttqex.f button.inc
buttqit.o: buttqit.f button.inc
buttqpr.o: buttqpr.f button.inc
buttqxb.o: buttqxb.f button.inc
buttqyb.o: buttqyb.f button.inc
buttqytext.o: buttqytext.f button.inc
buttsbr.o: buttsbr.f button.inc
buttscf.o: buttscf.f button.inc
buttsch.o: buttsch.f button.inc
buttsex.o: buttsex.f button.inc
buttsit.o: buttsit.f button.inc
buttspr.o: buttspr.f button.inc
buttsxb.o: buttsxb.f button.inc
buttsyb.o: buttsyb.f button.inc
buttsytext.o: buttsytext.f button.inc
ifbutton.o: ifbutton.f button.inc
rpgband.o: rpgband.f button.inc
rpgbegin.o: rpgbegin.f button.inc
rpgbegok.o: rpgbegok.f button.inc
rpgenv.o: rpgenv.f button.inc
rpgeras.o: rpgeras.f button.inc
rpgerasb.o: rpgerasb.f button.inc
rpgerasw.o: rpgerasw.f button.inc
# Target to clean object modules
clean:    $(OBJECTS)
	rm -f $(OBJECTS)
# Target to touch source modules
touch:
	touch $(FSOURCE)
