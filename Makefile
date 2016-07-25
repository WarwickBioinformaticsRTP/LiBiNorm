################################################################################
#
# Make all or Make debug makes the debug version
# Make release to make the release version
#
################################################################################
# Standard rules
CCC = g++

# Flags required buy all stages of C++ compiler
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
#	Instructions for building release and debug object files  These are dependant on the Makefile so Makefile changes
#	force a rebuild.  

$(BUILD)/%.o : %.cpp $(PCHOUT) Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<

$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)%.o : $(BIOINFORMATICSLIBDIR)%.cpp Makefile
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
#  All of the .cpp files in LIBDIR are included in the build and the objects are in <build>/mcmcLib
#  The files that are from external directories are handled slightly differently so that their object files
#  are still within the $(BUILD) directory and are deleted with a make clean

LIBINORMSRC = LiBiNormSrc/LiBiNorm.cpp

LIBINORMSRCEX = $(addprefix $(LIBINORMSRCDIR), LiBiCount.cpp LogLiklihoods.cpp GtfFileEx.cpp transcriptData.cpp Regions.cpp) 

BIOLIBSRC = $(shell find $(BIOINFORMATICSLIBDIR) -name *.cpp)
MCMCLIBSRC =  $(shell find $(MCMCLIBDIR) -name *.cpp)


SOURCES = $(LIBINORMSRC) $(MCMCLIBSRC) $(LIBINORMSRCEX)
        
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
$(BUILD)/LiBiNormSrc/LiBiNorm.o: mcmcLib/params.h LiBiNormSrc/LogLiklihoods.h
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
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/LiBiCount.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/GtfFileEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: LiBiNormSrc/LogLiklihoods.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: mcmcLib/mcmc.h mcmcLib/dataVec.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: LiBiNormSrc/GtfFileEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/BamReader.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/api_global.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/BamAux.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/BamConstants.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/BamIndex.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamHeader.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamProgram.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bamtools/api/SamSequence.h
$(BUILD)/LiBiNormSrc/transcriptData.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/transcriptData.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/transcriptData.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/transcriptData.o: LiBiNormSrc/transcriptData.h
$(BUILD)/LiBiNormSrc/transcriptData.o: mcmcLib/dataVec.h
$(BUILD)/LiBiNormSrc/Regions.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/BamReader.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/api_global.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/BamAux.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/BamConstants.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/BamIndex.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamHeader.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamProgram.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/api/SamSequence.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/codFile.h
$(BUILD)/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/../bioinformaticsLib/codFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/fastaFile.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/../bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/parser.h
$(BUILD)/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/printEx.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/parser.h
$(BUILD)/../bioinformaticsLib/gtfFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/parser.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/../bioinformaticsLib/parser.o: ../bioinformaticsLib/parser.h
$(BUILD)/../bioinformaticsLib/printEx.o: ../bioinformaticsLib/printEx.h
$(BUILD)/../bioinformaticsLib/printEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/../bioinformaticsLib/printEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/../bioinformaticsLib/smithWaterman.o: ../bioinformaticsLib/smithWaterman.h

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
