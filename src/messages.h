/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1988, 1992, 1993 Keith Refson
 
This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 
In other words, you are welcome to use, share and improve this program.
You are forbidden to forbid anyone else to use, share and improve
what you give them.   Help stamp out software-hoarding!  */
/*
 * $Header$
 */

#ifndef SYSRD   /* Skip if already defined				      */

		/* Severities to be passed to 'message'.		      */
#define		INFO		0
#define		WARNING		1
#define		ERROR		2
#define		FATAL		3
#define	NULLP		(char*)0
#define	NULLI		(int*)0

#define	BACKUP	"Simulation restarted from backup file %s"
#define RNFAIL	"Rename to \"%s\" failed -- file left in \"%s\"\n    -  %s"
#define SYSRD	"reading system specification file"
#define TOMANY	"Too many species in system specification file"
#define NOSPEC	"failed to read any molecular species info"
#define NONUM	"number of molecules of %s not specified"
#define NOMOLS	"invalid number (%d) of molecules of %s"
#define NOSITE	"invalid number (%d) of sites on molecule %s"
#define INVSID	"invalid site id %d"
#define NCONF	"site %d - name conflicts with previously set value %s"
#define CCONF	"site %d - charge conflicts with previously set value %f"
#define MCONF	"site %d - mass conflicts with previously set value %f"
#define	INVMAS	"site %d - negative value (%g) supplied for mass"
#define MISSCO	"only %d site co-ordinates specified for site %d"
#define BROKEN	"internal error - sscanf returned %d items"
#define NOMASS	"no mass specified for site id %d"
#define NOCGRG	"no charge specified for site id %d"
#define NONAME	"no name specified for site id %d"
#define NOTUSD	"unused site identifier %d"
#define	UNKPOT	"unknown type of potential %s (lennard-jones, buckingham, mcy, generic, hiw, morse)"
#define NOPAIR	"site id pair is required"
#define NOPOTP	"insufficient potential parameters supplied (%d needed)"
#define IDOUTR	"site id %d out of range"
#define EXTPOT	"potentials specified for unused sites - ignored"
#define DUPPOT	"potentials already set for this site pair"
#define NOPOT	"no potential parameters given between sites %d and %d"
#define ERRS	"system specification file contains %d error%c"
#define SUCCES	"system specification file successfully read in"
#define LATTIC	"system successfully initialised from lattice start"
#define	REFORM	"wrong length record in restart file at byte %ld\
 - found %lu, expected %lu"
#define	REREAD	"read error on restart file at byte %ld\n    -  %s"
#define	REEOF	"unexpected end of file reading restart file"
#define	REWRT	"write error on restart file\n    -  %s"
#define OLDRST  "restart file version %s was written by wrong Moldy version - can not read"
#define	NOVAL	"no value associated with name \"%s\""
#define	NOTFND	"keyword \"%s\" not found"
#define	BADVAL	"value \"%s\" is wrong type for keyword \"%s\""
#define	INVVAL	"value \"%d\" for keyword \"%s\" is invalid or out of range"
#define	INVVLF	"value \"%g\" for keyword \"%s\" is invalid or out of range"
#define	ERRCON	"control file contains %d error%c"
#define SUCCON  "control file read in successfully"
#define	BADUNI	"Overflow during conversion of units - ln(scale) = %f"
#define	ZMASS	"Mass of %s molecule (%f) is less than 1 amu"
#define	OCFAIL	"Failed to open file \"%s\" for reading control info"
#define ODFAIL	"Failed to open file \"%s\" for reading system specification\n    -  %s"
#define ORFAIL	"Failed to open file \"%s\" for reading restart configuration\n    -  %s"
#define OSFAIL	"Failed to open file \"%s\" for writing restart configuration\n    -  %s"
#define	OOFAIL	"Failed to open file \"%s\" as main output file\n    -  %s"
#define SEFAIL  "Rewind or seek failed on file  \"%s\"\n    -  %s"
#define	RESUCC	"restart file \"%s\" successfully read in"
#define	NEWTS	"Changing from old timestep %g to new %g"
#define	NSPCON	"Number of species in sysdef file (%d) and restart file (%d)\
must be the same"
#define	NMLCON	"Number of molecules of %s in sysdef and restart files must\
be the same (%d vs %d)"
#define	NDFCON	"Molecules of %s in sysdef and restart files must have same \
rotational degrees of freedom (%d vs %d)"
#define CROWDED "Cells too large - multiple occupation (mol=%d, cell=%d)"
#define NABORS  "Neighbour list contains %d cells"
#define TONAB	"Too many sites (%d) for neighbour list arrays (%d slots) \n\
  -  Increase NMULT in \"force.c\""
#define CUTOFF  "Cutoff radius > %d * cell dimension.\n\
  -  Increase table size IMCELL_XTRA in \"force.c\""
#define CUTRDF  "RDF limit > %d * cell dimension.\n\
  -  Increase NSH in \"force.c\""
#define TOOCLS	"Sites %d and %d closer than %fA."
#define TOODIM  "Arralloc request for %d dimensions - max %d"
#define INSIDE  "Array bounds [%d...%d] inside out (from ARRALLOC)"
#define UNKPTY  "KERNEL called with unknown potential type %d"
#define NOCELL  "Not enough values to specify unit cell\n\
     expected 3 lengths, 3 angles and number of unit cells in MD cell."
#define INVCEL  "Invalid unit cell parameters - must be +ve and angles < 180"
#define UNKSPE  "\"%s\" is not recognised as a molecular species"
#define	FEWCOO  "Too few co-ordinates for species \"%s\" - 3 needed"
#define	FEWQUA  "Too few quaternions for species \"%s\" - 4 needed"
#define FRACCO  "Fractional co-ordinates must be in range [0,1) - (%f,%f,%f)"
#define QNORM   "Quaternion (%f,%f,%f,%f) is not normalised"
#define QNORM2  "Quaternion %d (%f,%f,%f,%f) - normalisation error"
#define Qconst  "Quaternion %d - constraint error (%g)"
#define NIMOLS  "Wrong number of molecules of %s - given %d, expected %d"
#define INITER  "Initialisation file contains %d error%c"
#define ROTLEN  "Length (%d) not a multiple of quaternion number (%d) in \
\"rotate\""
#define OVRLAP  "%s - Result matrix overlaps input"
#define OVRLP1  "mat_vec_mul - Input and output vectors overlap \
(in=%x, out=%x, len=%x)"
#define SNGMAT  "%s - matrix is singular"
#define AVNOC   "%s - init_averages has not been called"
#define AVBNDS  "add_average - offset (%d) out of bounds for type %s"
#define NEGVAR  "%s - variance < 0 (%f) for type %s"
#define NOMEM   "Memory allocation fails at line %d in \"%s\" \
(%d item%sof %lu bytes)"
#define	DUMPST	"Started dumping data to file \"%s\" at timestep %ld"
#define CONTIG	"Dump file \"%s\" and restart file do not match"
#define DUMPTS  "Dump file \"%s\"'s timestep (%ld) does not match current timestep (%ld)."
#define DRESET  "Problems prevent continuing existing dump.  Starting new sequence."
#define	SHTDMP	"Records missing from dump file \"%s\" - found %d, expected %d"
#define LNGDMP	"Extra records in dump file \"%s\" - found %d, expected %d"
#define	CORUPT	"Dump file \"%s\" corrupt - expected %d bytes, found %d"
#define DOERRR	"Failed to open dump file \"%s\"\n    -  %s"
#define DRERR	"Read from dump file \"%s\" failed\n    -  %s"
#define DHDERR	"Failed to read header from dump file \"%s\""
#define RDHERR  "Internal error in read_sump_header (%lu,%lu)"
#define DMPALT	"Dump-level altered on restart: new dump run started"
#define DMPAL2	"Dump-interval or maxdumps altered on restart: new dump run started"
#define DMPEXS	"File \"%s\" exists - Dumps will be written to \"%s\""
#define DOERRW	"Failed to open dump file \"%s\" for writing\n    -  %s"
#define DWERR	"Write to dump file \"%s\" failed\n    -  %s"
#define MUFAIL  "Unable to mutate dump file name \"%s\" (in %d attempts)\n"
#define UNKEY   "Unknown keyword \"%s\""
#define	SYSCHG	"System has net electric charge of %.2g - \
correction of %g kJmol(-1) added to self energy"
#define	FRACHG	"Framework has net electric charge of %.2g - \
correction of %g kJmol(-1) added to self energy"
#define WDPTR	"Arralloc called for object with size (%lu) not an integral \
number of words"
/* #define NCNVRG  "Do_step: velocities failed to converge after %d iterations \
- dist = %f"*/
#define INRVSN	"Init_averages: Invalid RCS version \"%s\""
#define INDVSN  "Read_dump_header:  Invalid RCS version \"%s\""
#define OLDVSN	"Attempt to append to incompatible version of dump files (RCS version \"%s\")"
#define CPOTFL  "Not enough space to store potential parameters from restart\
 file.  Increase NPOTP in \"defs.h\" from %d to at least %d and recompile"
#define LOCKED   "Another instance of MOLDY is accessing the same %s files.\n\
     Delete lockfile \"%s\" if the previous run crashed."
#define LOCKFL  "Failed to create lockfile \"%s\""
#define XDRFL	"XDR call failed"
#define GPFAIL  "Get file position failed on file  \"%s\"\n    -  %s"
#endif
#define SYSEOF	"Unexpected end of file encountered while reading %s"
#define MAXCUT  "Optimum cutoff of %2.2f A is larger than allowed maximum \
of %2.2f.\n Maximum will be used but accuracy of Ewald sum is compromised."
#define GANDP   "Application of Parinello-Rahman algorithm along with \
Gaussian thermostat would interfere with each other."
#define	NOKEY	"no keyword found on line"
#define NOCUT   "No real-space cutoff specified or cutoff invalid"
#define BODGUP  "Reading bad 2.10 restart file - thermostat parameters invalid.\
\n     Setting ttmass and rtmass to %f and const-temp=0"
#define DESYNC  "Trajectories on parallel threads are diverging.\n     Thread %d: %20.17g != %20.17g"
#define RDFALT  "One of the RDF sampling parameters changed on restart\n     stored RDF data has been discarded"
#define RDFDIS  "RDF output frequency parameters changed on restart\n     Excess stored RDF data may be discarded"
#define CNVRDF "Converting stored RDF data from old integer format (%s) to new format"
#define RDFFIS "Can not convert partial RDF data from restart file version %s on this machine because sizeof(int)(%d) != sizeof(float)\n -- This might be possible on a different platform (eg Intel x86)"
#define DPSCHG  "Surface dipole term switched ON but some species are charged\n     Expect non-conservation of energy and large, impulsive forces!"
#define FIXROL "Changing roll-interval when converting pre 2.19 restart file is not supported"
#define SAVINC "Save file would overwrite input restart file.  Filename changed to %s"
#define IVEWAC "Ewald Sum accuracy parameter (%g) invalid.  Should be > 0"
#define ERFACC "Error function approximation for real-space Ewald Sum can not meet requested level of accuracy (%g)\n     Algorithm changes in kernel.c would be needed.  See Abramowitz and Stegun p299."
#define EWINAC "Explicitly specified Ewald sum parameters give accuracy parameter of %g, not %g as requested"
