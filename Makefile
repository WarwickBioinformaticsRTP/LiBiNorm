################################################################################
#
# Make all or Make debug makes the debug version
# Make release to make the release version
#
################################################################################
# Standard rules
CCC = g++

# Flags required buy all stages of C++ compiler
# CCCALLFLAGS= -std=c++0x
CCCALLFLAGS= -std=gnu++11 -Wno-format-security

# Directory information

BAMTOOLSDIR = ../bamtools/
BIOINFORMATICSLIBDIR = ../bioinformaticsLib/
BIOINFORMATICSLIBDESTDIR = bioinformaticsLib/
MCMCLIBDIR = mcmcLib/
LIBINORMSRCDIR = LiBiNormSrc/
RELDIR = Release
DEBUGDIR = Debug

ifeq "$(findstring release, $(MAKECMDGOALS))" ""
BUILD=$(DEBUGDIR)
CCCALLFLAGS += -g3 -O0 -D_DEBUG -Wall -Wno-unknown-pragmas \
	 -Wno-reorder -Wno-format-security
else
BUILD=$(RELDIR)
CCCALLFLAGS += -O3
endif

INCLUDES= -I$(BAMTOOLSDIR) -I$(BIOINFORMATICSLIBDIR) -I$(MCMCLIBDIR)

##################################################################
#
#	Instructions for building the pre compiled header.   This is not 'cleaned', but is rebuilt if
#	the pch.h file or the Makefile is changed.   A different one is created for each build type in the .h.pch 
#	directory.   g++ works out which one to use based on its contents not the name


$(PCHOUT): $(addprefix ../bamtools/, shared/bamtools_global.h \
	$(addprefix api/, BamWriter.h  api_global.h BamAux.h BamReader.h BamAlignment.h \
	BamConstants.h BamIndex.h SamHeader.h SamProgramChain.h SamProgram.h \
	SamReadGroupDictionary.h SamReadGroup.h SamSequenceDictionary.h SamSequence.h ) )


##################################################################
#
#	Instructions for building release and debug object files  These are dependant on the Makefile so Makefile changes
#	force a rebuild.  ALso dependant on the precompiled header so the pch gets built first and a change to the pch
#	forces a rebuild

$(BUILD)/%.o : %.cpp $(PCHOUT)
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<

$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)%.o : $(BIOINFORMATICSLIBDIR)%.cpp
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<


################################################################################
# Outputs of this Makefile

LIBINORM = LiBiNorm

ifeq ($(OS),Windows_NT)
EXE = .exe
endif

LIBINORMEXE = $(BUILD)/$(LIBINORM)$(EXE) 

TARGS =  $(LIBINORMEXE) 

################################################################################
# Libraries to be linked.   Assumes that there is a RElease and Debug version of Bamtools

LIBS        = -lbamtools -lz 
LIBPATH     = -L$(BAMTOOLSDIR)$(BUILD)

################################################################################
# Source files and the build specific outputs
#  All of the .cpp files in LIBDIR are included in the build and the objects are in <build>/libdir
#  The files that are from external directories are handled slightly differently so that there object files
#  are still within the $(BUILD) directory

LIBINORMSRC = LiBiNormSrc/LiBiNorm.cpp

BIOLIBSRC = $(shell find $(BIOINFORMATICSLIBDIR) -name *.cpp)
MCMCLIBSRC =  $(shell find $(MCMCLIBDIR) -name *.cpp)

LIBINORMSRCEX = $(addprefix $(LIBINORMSRCDIR), LiBiCount.cpp LogLiklihoods.cpp GtfFileEx.cpp transcriptData.cpp Regions.cpp) 


SOURCES = $(LIBINORMSRC) $(MCMCLIBSRC)
        
COREOBJS :=  $(LIBINORMSRCEX:%.cpp=$(BUILD)/%.o) $(MCMCLIBSRC:%.cpp=$(BUILD)/%.o)  $(BIOLIBSRC:%.cpp=$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)%.o)

################################################################################
# For building necessary outout directories.  DIRMARKERS are a set of empty files called .z.   If they are not present then the
# mkdir, touch code will make the directory and insert the file.

DIRECTORIES = $(addprefix $(BUILD)/,  $(BIOINFORMATICSLIBDESTDIR) $(dir $(SOURCES)))

DIRMARKERS = $(addsuffix .z , $(DIRECTORIES) )

%.z:
	mkdir -p $(@D)
	touch $@


################################################################################
# The main builds



debug : all
	
release : all    

all:  $(DIRMARKERS) $(TARGS)
	@echo "%% $(BUILD) sysmedibd code built"

#	This is the format for a manual final make rule
$(LIBINORMEXE) :$(BUILD)/$(LIBINORMSRC:%.cpp=%.o) $(COREOBJS)
	$(CCC)  $^ -o $@ $(CCCALLFLAGS) $(INCLUDES) $(LIBPATH) $(LIBS)

# and this could be used as a common rule to build them all.   Note that the $< and $@ dont seem to work
#	requires a fairly recent version of make to work though
#define MAINRULE =
#Release/$(FILEROOT) : $(BUILD)/$(FILEROOT)Src/$(FILEROOT).o $(COREOBJS)
#	$(CCC)  $(BUILD)/$(FILEROOT)Src/$(FILEROOT).o $(COREOBJS) -o Release/$(FILEROOT).exe $(CCCALLFLAGS) $(LIBPATH) $(LIBS)
#endef

#$(foreach FILEROOT,$(FILETEST),$(eval  $(MAINRULE)))
	
clean : 
	rm -r -f $(RELDIR)
	rm -r -f $(DEBUGDIR)

#	do make depend to update dependancies.  This makes a dependancy list that is dynamically dependant 
#	on the build type.  There are two make depends, one for all of the sources within this directory (SOURCES)
#	and one for the sources that are in another library directory (FASTASRC)

depend :
	makedepend  -Y $(CCCAALLFLAGS) $(INCLUDES) $(SOURCES) -p'$$(BUILD)/'
	makedepend  -Y -a $(CCCAALLFLAGS) $(INCLUDES) $(BIOLIBSRC) -p'$$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)'	
		
# DO NOT DELETE THIS LINE -- make depend depends on it.

$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: mcmcLib/mcmc.h mcmcLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: mcmcLib/params.h LiBiNormSrc/logLiklihoods.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiNorm.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/transcriptData.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiCount.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/GtfFileEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/api/SamSequence.h
$(BUILD)/mcmcLib/dataVec.o: mcmcLib/dataVec.h
$(BUILD)/mcmcLib/mcmc.o: mcmcLib/mcmc.h mcmcLib/dataVec.h mcmcLib/params.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/mcmcLib/params.o: mcmcLib/params.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/mcmcLib/params.o: mcmcLib/dataVec.h

$(BUILD)/bioinformaticsLib/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/codFile.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/fastaFile.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/parser.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/parser.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/printEx.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/printEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/printEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/smithWaterman.o: ../bioinformaticsLib/smithWaterman.h
