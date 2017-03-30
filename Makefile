FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check -mcmodel=medium
#FFLAGS = -O2 -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace
#FFLAGS = -O3 -march=native -fimplicit-none -Wall -Wline-truncation -fwhole-file -std=f2008
SOURCES = bsort.f calc_kinf.f canon_min.f canon_rates.f canon_writeout.f check_input.f chemformula.f chkdens.f covsrt.f dnt.f dsort.f eis.f element.f finalprint.f find_jmax.f find_ming.f find_tmins.f funcs.f gaussj.f get_viblo.f hsort.f indexx.f jaread_input.f javerage.f jthermavg.f lfit.f lowerc.f main.f micro_rates.f micro_writeout.f multiwell_write.f nduhrlev.f nghrlev.f nkinfint.f nkrotlev.f nmorlev.f nrgconvert.f nshrlev.f parabfit.f japrewrite.f qele.f qhinf.f qroq.f qvibf.f rank.f read_input.f sestamp.f sort_input.f stamp.f sterab.f sterabj.f super_mol.f tghrlev.f tqstop.f tshrlev.f veff.f write_input.f write_mat.f zpe_calc.f
OBJECTS = $(SOURCES: .f=.o)
EXECUTABLE = ktools

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $@
	cp $(EXECUTABLE) ../../bin/ 
.f.o:
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)
	rm -f ../../bin/$(EXECUTABLE)
