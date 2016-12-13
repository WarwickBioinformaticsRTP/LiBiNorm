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
BAMTOOLSDESTDIR = bamtools/
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

$(BUILD)/%.o : %.cpp Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<

$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)%.o : $(BIOINFORMATICSLIBDIR)%.cpp Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<

$(BUILD)/$(BAMTOOLSDESTDIR)%.o : $(BAMTOOLSDIR)%.cpp Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -Wno-sign-compare -o $@ $<

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

LIBS        = -lz -pthread
LIBPATH     = 

################################################################################
# Source files and the build specific outputs
#  All of the .cpp files in MCMCLIBDIR are included in the build and the objects are in <build>/mcmcLib
#  The files that are from the external directories are handled slightly differently so that their object files
#  are also placed within the $(BUILD) directory and so are deleted with a make clean

LIBINORMSRC = LiBiNormSrc/LiBiNorm.cpp

LIBINORMSRCEX = $(addprefix $(LIBINORMSRCDIR), LiBiDedup.cpp LiBiCount.cpp LogLiklihoods.cpp GtfFileEx.cpp transcriptData.cpp Regions.cpp MakeFastq.cpp) 

BIOLIBSRC = $(shell find $(BIOINFORMATICSLIBDIR) -name *.cpp)
MCMCLIBSRC =  $(shell find $(MCMCLIBDIR) -name *.cpp)
BAMTOOLSFILES = api/BamAlignment.cpp api/BamReader.cpp api/BamWriter.cpp api/SamHeader.cpp api/SamProgram.cpp api/SamProgramChain.cpp \
    api/SamReadGroup.cpp api/SamReadGroupDictionary.cpp api/SamSequence.cpp api/SamSequenceDictionary.cpp \
    api/BamMultiReader.cpp \
    api/internal/bam/BamHeader_p.cpp api/internal/bam/BamMultiReader_p.cpp api/internal/bam/BamRandomAccessController_p.cpp \
    api/internal/bam/BamReader_p.cpp api/internal/bam/BamWriter_p.cpp api/internal/index/BamIndexFactory_p.cpp \
    api/internal/index/BamStandardIndex_p.cpp api/internal/index/BamToolsIndex_p.cpp api/internal/sam/SamFormatParser_p.cpp \
    api/internal/sam/SamFormatPrinter_p.cpp api/internal/sam/SamHeaderValidator_p.cpp \
    api/internal/utils/BamException_p.cpp api/internal/io/BamDeviceFactory_p.cpp \
    api/internal/io/BamFile_p.cpp api/internal/io/BamFtp_p.cpp \
    api/internal/io/BamHttp_p.cpp api/internal/io/BamPipe_p.cpp api/internal/io/BgzfStream_p.cpp \
    api/internal/io/ByteArray_p.cpp api/internal/io/HostAddress_p.cpp \
    api/internal/io/HostInfo_p.cpp api/internal/io/HttpHeader_p.cpp \
    api/internal/io/ILocalIODevice_p.cpp api/internal/io/RollingBuffer_p.cpp \
    api/internal/io/TcpSocket_p.cpp api/internal/io/TcpSocketEngine_p.cpp \
    api/internal/io/TcpSocketEngine_unix_p.cpp \
    toolkit/bamtools_sort.cpp utils/bamtools_options.cpp
    
BAMTOOLSSRC = $(addprefix $(BAMTOOLSDIR), $(BAMTOOLSFILES) )    

SOURCES = $(LIBINORMSRC) $(MCMCLIBSRC) $(LIBINORMSRCEX)
        
COREOBJS :=  $(LIBINORMSRCEX:%.cpp=$(BUILD)/%.o) \
	$(MCMCLIBSRC:%.cpp=$(BUILD)/%.o)  \
	$(BAMTOOLSSRC:%.cpp=$(BUILD)/$(BAMTOOLSDESTDIR)%.o) \
	$(BIOLIBSRC:%.cpp=$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)%.o)


	
################################################################################
# For building necessary outout directories.  DIRMARKERS are a set of empty files called .z.   If they are not present then the
# mkdir, touch code will make the directory and insert the file.

BAMDIRS =  $(addprefix $(BAMTOOLSDESTDIR), $(dir $(BAMTOOLSSRC) ) )

DIRECTORIES = $(addprefix $(BUILD)/,  $(BAMDIRS) $(BIOINFORMATICSLIBDESTDIR) $(dir $(SOURCES)))

DIRMARKERS = $(addsuffix .z , $(DIRECTORIES) )

%.z:
	mkdir -p $(@D)
	touch $@


################################################################################
# The main builds

debug : all
	
release : all    

all:  $(DIRMARKERS) $(TARGS)
	@echo "%% $(BUILD) LiBiNorm code built"

#	The final make rule
$(LIBINORMEXE) :$(BUILD)/$(LIBINORMSRC:%.cpp=%.o) $(COREOBJS)
	$(CCC)  $^ -o $@ $(CCCALLFLAGS) $(INCLUDES) $(LIBPATH) $(LIBS)

#	All intermediate and final build products go in one directory and its sub directories 
#	for each of the build types, making a make clean very simple
	
clean : 
	rm -r -f $(RELDIR)
	rm -r -f $(DEBUGDIR)

#	do make depend to update dependancies.  This makes a dependancy list that is dynamically dependant 
#	on the build type.  There are three make depends, one for all of the sources within this directory (SOURCES)
#	and one for the sources that are in other library directories (BIOLIBSRC & BAMTOOLSSRC)

depend :
	makedepend  -Y $(CCCAALLFLAGS) $(INCLUDES) $(SOURCES) -p'$$(BUILD)/'
	makedepend  -Y -a $(CCCAALLFLAGS) $(INCLUDES) $(BIOLIBSRC) -p'$$(BUILD)/$(BIOINFORMATICSLIBDESTDIR)'	
	makedepend  -Y -a $(CCCAALLFLAGS) $(INCLUDES) $(BAMTOOLSSRC) -p'$$(BUILD)/$(BAMTOOLSDESTDIR)'		

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: mcmcLib/mcmc.h mcmcLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: mcmcLib/params.h LiBiNormSrc/LogLiklihoods.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiNorm.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/transcriptData.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiCount.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/GtfFileEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/genbankFile.h
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
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiDedup.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/MakeFastq.h
$(BUILD)/mcmcLib/dataVec.o: mcmcLib/dataVec.h
$(BUILD)/mcmcLib/dataVec.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/mcmcLib/dataVec.o: ../bioinformaticsLib/printEx.h
$(BUILD)/mcmcLib/dataVec.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/mcmcLib/mcmc.o: mcmcLib/mcmc.h mcmcLib/dataVec.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/printEx.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/inQuotes.h mcmcLib/params.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/mcmcLib/params.o: mcmcLib/params.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/mcmcLib/params.o: mcmcLib/dataVec.h ../bioinformaticsLib/libCommon.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/printEx.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: LiBiNormSrc/LiBiDedup.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/api/BamWriter.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/LiBiCount.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/GtfFileEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/printEx.h
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
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/transcriptData.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: mcmcLib/dataVec.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: LiBiNormSrc/LogLiklihoods.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: mcmcLib/mcmc.h mcmcLib/dataVec.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LogLiklihoods.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: LiBiNormSrc/GtfFileEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/gtfFile.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/GtfFileEx.o: ../bioinformaticsLib/printEx.h
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
$(BUILD)/LiBiNormSrc/transcriptData.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/transcriptData.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/Regions.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/inQuotes.h
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
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/parser.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: LiBiNormSrc/MakeFastq.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/BamReader.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/api_global.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/BamAux.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/BamConstants.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/BamIndex.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamHeader.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamProgram.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/api/SamSequence.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/fastaFile.h

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
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/parser.h
$(BUILD)/bioinformaticsLib/../bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/libCommon.h
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

$(BUILD)/bamtools/../bamtools/api/BamAlignment.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/BamAlignment.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/BamAlignment.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/BamAlignment.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/BamAlignment.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/BamReader.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/BamWriter.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/internal/bam/BamWriter_p.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/BamWriter.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamConstants.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/internal/sam/SamFormatParser_p.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/internal/sam/SamFormatPrinter_p.h
$(BUILD)/bamtools/../bamtools/api/SamHeader.o: ../bamtools/api/internal/sam/SamHeaderValidator_p.h
$(BUILD)/bamtools/../bamtools/api/SamProgram.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/SamProgram.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamProgram.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamProgramChain.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/SamProgramChain.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamProgramChain.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamProgramChain.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroup.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroup.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroup.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroupDictionary.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroupDictionary.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroupDictionary.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamReadGroupDictionary.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/SamSequence.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/SamSequence.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamSequence.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamSequenceDictionary.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/SamSequenceDictionary.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/SamSequenceDictionary.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/SamSequenceDictionary.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/BamMultiReader.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/internal/bam/BamMultiReader_p.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/internal/bam/BamMultiMerger_p.h
$(BUILD)/bamtools/../bamtools/api/BamMultiReader.o: ../bamtools/api/algorithms/Sort.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamHeader_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/BamMultiReader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/SamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/algorithms/Sort.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/internal/bam/BamMultiReader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamMultiReader_p.o: ../bamtools/api/internal/bam/BamMultiMerger_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/internal/index/BamIndexFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/index/BamStandardIndex_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/index/BamToolsIndex_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamReader_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/internal/bam/BamWriter_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/bam/BamWriter_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/internal/index/BamIndexFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/internal/index/BamStandardIndex_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamIndexFactory_p.o: ../bamtools/api/internal/index/BamToolsIndex_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/index/BamStandardIndex_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamStandardIndex_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/index/BamToolsIndex_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/index/BamToolsIndex_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/internal/sam/SamFormatParser_p.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatParser_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/api/internal/sam/SamFormatPrinter_p.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/internal/sam/SamHeaderValidator_p.h
$(BUILD)/bamtools/../bamtools/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/api/internal/sam/SamHeaderVersion_p.h
$(BUILD)/bamtools/../bamtools/api/internal/utils/BamException_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/internal/io/BamFile_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/internal/io/BamFtp_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/internal/io/BamHttp_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamDeviceFactory_p.o: ../bamtools/api/internal/io/BamPipe_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFile_p.o: ../bamtools/api/internal/io/BamFile_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFile_p.o: ../bamtools/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFile_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFile_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFile_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/internal/io/BamFtp_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamFtp_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/BamHttp_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/HttpHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamHttp_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamPipe_p.o: ../bamtools/api/internal/io/BamPipe_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamPipe_p.o: ../bamtools/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamPipe_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamPipe_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BamPipe_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/BgzfStream_p.o: ../bamtools/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ByteArray_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ByteArray_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ByteArray_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostAddress_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostAddress_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostAddress_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostInfo_p.o: ../bamtools/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostInfo_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostInfo_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostInfo_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HostInfo_p.o: ../bamtools/api/internal/io/NetUnix_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HttpHeader_p.o: ../bamtools/api/internal/io/HttpHeader_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HttpHeader_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/HttpHeader_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ILocalIODevice_p.o: ../bamtools/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ILocalIODevice_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ILocalIODevice_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/ILocalIODevice_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/RollingBuffer_p.o: ../bamtools/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/RollingBuffer_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/RollingBuffer_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/RollingBuffer_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocket_p.o: ../bamtools/api/internal/io/TcpSocketEngine_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/internal/io/TcpSocketEngine_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/TcpSocketEngine_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/IBamIODevice.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/../bamtools/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/api/internal/io/NetUnix_p.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/toolkit/bamtools_sort.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/toolkit/bamtools_tool.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamConstants.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/api_global.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/shared/bamtools_global.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamMultiReader.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamReader.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamAlignment.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamAux.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamConstants.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamIndex.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamHeader.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamProgramChain.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamProgram.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamReadGroup.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamSequenceDictionary.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/SamSequence.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/BamWriter.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/api/algorithms/Sort.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/utils/bamtools_options.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/utils/bamtools_variant.h
$(BUILD)/bamtools/../bamtools/toolkit/bamtools_sort.o: ../bamtools/utils/utils_global.h
$(BUILD)/bamtools/../bamtools/utils/bamtools_options.o: ../bamtools/utils/bamtools_options.h
$(BUILD)/bamtools/../bamtools/utils/bamtools_options.o: ../bamtools/utils/bamtools_variant.h
$(BUILD)/bamtools/../bamtools/utils/bamtools_options.o: ../bamtools/utils/utils_global.h
$(BUILD)/bamtools/../bamtools/utils/bamtools_options.o: ../bamtools/shared/bamtools_global.h
