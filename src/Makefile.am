AUTOMAKE_OPTIONS = no-dependencies dist-zip 

HFILES= structs.h defs.h messages.h xdr.h  list.h utlsup.h ReadDCD.h defconf.h specdata.h sginfo.h

UTILC = algorith.c alloc.c auxil.c ansi.c convert.c eigens.c kernel.c input.c matrix.c quaterns.c restart.c startup.c values.c xdr.c ReadDCD.c utlsup.c molout.c

UTILS =  moldyext dumpconv dumpext  manalyze bdist mdshak msd mdavpos mdbond ransub mdvaf syswrite

bin_PROGRAMS = moldy $(UTILS)

bin_SCRIPTS = plotrdf moldyext2

moldy_SOURCES = accel.c algorith.c alloc.c ansi.c auxil.c convert.c dump.c ewald.c force.c input.c eigens.c kernel.c leapfrog.c main.c matrix.c output.c quaterns.c rdf.c restart.c startup.c values.c xdr.c parallel.c  $(HFILES)

msd_SOURCES = msd.c $(UTILC)

mdshak_SOURCES = mdshak.c  $(UTILC) 

mdavpos_SOURCES = mdavpos.c  $(UTILC) 

mdbond_SOURCES = mdbond.c  $(UTILC) 

mdvaf_SOURCES = mdvaf.c  $(UTILC) 

moldyext_SOURCES = moldyext.c 

dumpconv_SOURCES = dumpconv.c xdr.c 

dumpext_SOURCES = dumpext.c xdr.c utlsup.c matrix.c

manalyze_SOURCES = manalyze.c 

bdist_SOURCES = bdist.c

ransub_SOURCES = ransub.c sgclib.c sginfo.h sgio.c spgrexp.c transparse.c readers.c $(UTILC) 

syswrite_SOURCES = syswrite.c sgclib.c sginfo.h sgio.c spgrexp.c transparse.c readers.c $(UTILC)

EXTRA_msd_SOURCES = getopt.c

FIGURES = fig_arralloc fig_dostep-a fig_dostep-b \
          fig_ewald fig_link-cell fig_main \
          fig_skewstart fig_startup-a fig_startup-b

FIGDIST  = $(FIGURES:%=figures/%.ftx) $(FIGURES:%=figures/%.ps) \
	   $(FIGURES:%=figures/%.pdf) $(FIGURES:%=figures/%-eepic.ftx) 

INFILES = 	control.water control.tip4p control.mgclh2o control.clay \
                control.argon control.tips2 control.quartz \
		tips2.in tip4p.in mgclh2o.in \
		argon.in quartz-vbst.in methane.in mcy.in \
		water-example.out tip4p-example.out \
		mgclh2o-example.out clay-example.out \
		argon-example.out tips2-example.out \
		quartz-example.out

EXTRA_DIST = figures moldy.tex moldy.sty moldy.perl moldy.bbl plotrdf moldyext2 $(INFILES) ewald-RIL.c LICENSE-SGINFO

dist-local:	
		bibtex moldy

CLEANFILES = moldy.aux moldy.lof moldy.log moldy.lot moldy.out moldy.toc   

DISTCLEANFILES = moldy.dvi moldy.pdf moldy.ps

COMPILE_NOOPT = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) 

COMPILE = $(COMPILE_NOOPT) $(OPT)

COMPILE_OPT2 = $(COMPILE_NOOPT) $(OPT2) $(INLINE)

#
#  Special options for performance-critical modules
#
auxil.o:	auxil.c defs.h
	$(COMPILE_NOOPT) $(OPTAUX) $(INLINE) -c  $<

ewald.o:	ewald.c structs.h defs.h
	$(COMPILE_OPT2) -c  $<

rdf.o:		rdf.c structs.h defs.h
	$(COMPILE_OPT2) -c  $<

kernel.o:	kernel.c structs.h defs.h
	$(COMPILE_OPT2) -c  $<

force.o:	force.c structs.h defs.h
	$(COMPILE_OPT2) -c  $<
#
# Lint targets
#
LCHECK =$(LINT) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) $(LFLAGS) 

lint:		
		$(LCHECK) $(moldy_SOURCES) $(LIBS)

lint-utils:		
		$(LCHECK) $(mdshak_SOURCES) 
		$(LCHECK) $(msd_SOURCES) 
		$(LCHECK) $(mdbond_SOURCES) 
		$(LCHECK) $(mdavpos_SOURCES) 
		$(LCHECK) $(ransub_SOURCES) 
		$(LCHECK) $(mdvaf_SOURCES) 
		$(LCHECK) $(moldyext_SOURCES) 
		$(LCHECK) $(dumpconv_SOURCES) 
		$(LCHECK) $(dumpext_SOURCES) 
		$(LCHECK) $(manalyze_SOURCES) 
		$(LCHECK) $(bdist_SOURCES) 
		$(LCHECK) $(syswrite_SOURCES) 
#
# Figures and Manual
#
moldy.dvi:	moldy.tex moldy.bbl moldy.aux $(FIGDIST)
		TEXINPUTS=figures: latex moldy
		TEXINPUTS=figures: latex moldy
		@touch moldy.bbl

moldy.pdf:	moldy.tex moldy.bbl moldy.aux $(FIGDIST)
		TEXINPUTS=figures: TEXPSHEADERS=figures: pdflatex moldy
		TEXINPUTS=figures: TEXPSHEADERS=figures: pdflatex moldy
		@touch moldy.bbl

moldy.aux:	moldy.tex
		TEXINPUTS=figures:: latex moldy

#moldy.bbl:	moldy.aux
#		-bibtex moldy

ewald-RIL.c:
	 cvs co -p -r RIL moldy/src/ewald.c > ewald-RIL.c


#
# Dependencies -- put in by hand since we don't use gcc.
#
ReadDCD.o: ReadDCD.c ReadDCD.h
accel.o: accel.c defs.h config.h structs.h
algorith.o: algorith.c defs.h config.h messages.h
alloc.o: alloc.c defs.h config.h messages.h
ansi.o: ansi.c defs.h config.h
bdist.o: bdist.c defs.h config.h list.h
convert.o: convert.c defs.h config.h structs.h messages.h
dump.o: dump.c defs.h config.h structs.h messages.h xdr.h
dumpconv.o: dumpconv.c defs.h config.h messages.h structs.h xdr.h
dumpext.o: dumpext.c defs.h config.h structs.h xdr.h utlsup.h messages.h
eigens.o: eigens.c defs.h config.h
input.o: input.c defs.h config.h structs.h messages.h
leapfrog.o: leapfrog.c defs.h config.h structs.h messages.h
main.o: main.c defs.h config.h structs.h messages.h
manalyze.o: manalyze.c defs.h config.h structs.h
matrix.o: matrix.c defs.h config.h messages.h
mdavpos.o: mdavpos.c defs.h config.h structs.h messages.h utlsup.h
mdbond.o: mdbond.c defs.h config.h structs.h messages.h list.h utlsup.h
mdshak.o: mdshak.c defs.h config.h structs.h utlsup.h
mdvaf.o: mdvaf.c defs.h config.h structs.h messages.h utlsup.h
moldyext.o: moldyext.c defs.h config.h
molout.o: molout.c defs.h config.h structs.h ReadDCD.h
msd.o: msd.c defs.h config.h structs.h messages.h utlsup.h
output.o: output.c defs.h config.h structs.h messages.h
parallel.o: parallel.c defs.h config.h structs.h messages.h
quaterns.o: quaterns.c defs.h config.h
ransub.o: ransub.c readers.c defs.h config.h structs.h messages.h utlsup.h \
        specdata.h
restart.o: restart.c defs.h config.h structs.h messages.h xdr.h
sgclib.o: sgclib.c sginfo.h
sgio.o: sgio.c sginfo.h
spgrexp.o: spgrexp.c sginfo.h specdata.h messages.h defs.h config.h
startup.o: startup.c defs.h config.h structs.h messages.h
syswrite.o: syswrite.c readers.c defs.h config.h structs.h messages.h utlsup.h \
        sginfo.h specdata.h
transparse.o: transparse.c sginfo.h
utlsup.o: utlsup.c defs.h config.h structs.h messages.h specdata.h
values.o: values.c defs.h config.h structs.h messages.h
xdr.o: xdr.c structs.h defs.h config.h xdr.h
