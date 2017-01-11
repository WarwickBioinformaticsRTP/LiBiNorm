#ifndef OPTIONS_H
#define OPTIONS_H


#define DEF_MAX_READS_FOR_PARAM_ESTIMATION 1000000
#define MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION 20000
#define MAX_READS_GENE 100
#define DEF_THREADS 3
#define MCMC_ITERATIONS 2000
#define DEFAULT_MODEL ModelBD
#define NUMBER_OF_MCMC_RUNS 3
#define DEFAULT_NORMALISATION_GENE_LENGTH 1000
#define MAX_GENE_LENGTH_FOR_NORM_PLOT 20000
#define DEFAULT_COUNT_MODE intersect_union
#define DEFAULT_FEATURE_TYPE_EXON "exon" 
#define DEFAULT_GTF_ID_ATTRIBUTE "gene_id"
#define DEFAULT_GFF_ID_ATTRIBUTE "Genbank"
#define END_LENGTH_SEARCHED_FOR_OPTIMAL_PARAMETERS 200

//	Use this to add the mode which creates fastq files based on bam files with artificial problems
// #define MAKE_FASTQ_MODE

//	Use this to add the mode where duplicates in bam files can be removed
// #define DEDUP_MODE

//	Use this mode to run a model with specific parameters.  This affects how the values are set
//	in ModelParameters.cpp and also ensures that the first N reads are use in data are loaded
//	in void rnaPosVec::selectAtMost(size_t s) in GeneCountData.cpp
//#define PRESET_VALUES {-0.0224725,	1.971188,	-4.373709,	-3.3618775,	0.8589705}

//	Some of the code in ModelData.cpp has also been writtent using vectors which is slower but the code
//	more closely matches the MATLAB code
//#define VECTOR_MATHS

//	The original MATLAB code had an error in setting the initial values for mcmc runs which this 
//	simulates (ModelData.cpp)
// #define SIMULATE_MATLAB_BUG

//
//	Use this option to use the parameters associated with the most likly parameter set
//	https://sciencehouse.wordpress.com/2010/06/23/mcmc-and-fitting-models-to-data/
//	rather than the median values
// #define USE_PARAMS_FROM_LOWEST_LL


#endif