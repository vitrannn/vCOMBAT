/**
 * @file   main.c
 * @version 5.0
 * @revision Vi Tran (original version written by T on 24 April 2015) 
 * @updated  Feb-2020
 */

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <yaml.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <math.h>
#include "carg_parser.h"
#include "full_model.h"
#include "addon.h"
#include "base_simulation.h"
#include "tuberculosis_simulation_config.h"

static const char* invocationName = NULL; ///< The full invocation text for the program including leading directories
static char* programName; ///< The basename of the program invocation
int verbose = 0; ///< Specifies whether the program should display verbose output


static void argParserShowError(const char* const msg, const int errorCode, const char help);
static void argParserInternalError(const char* const msg);
static const char* optname(const int code, const struct ap_Option options[]);
static void displayHelp(const char* programName);
static int outputResultsToFile(SimulationParameters sParam, ModelParameters mParam, SimulationResults results, FILE* oHandle);

#define DEFAULT_DUMMY -12345 ///< Definition of dummy value for marking those defaults which must be dynamically calculated
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define VERSION_STRING STR(TuberculosisSimulation_VERSION_MAJOR) "." STR(TuberculosisSimulation_VERSION_MINOR)
/**
 * Program entry procedure.
 * 
 * @param  argc Argument count
 * @param  argv Argument array
 *
 * @return     EXIT_SUCCESS on succesful completion
 */
int main(const int argc, char** argv) {
	// Set up argument parsing structures
	int i;
	clock_t t;
	int argIdx;
	struct Arg_parser parser;
	double* stateVector = NULL;
	const char* outputFile = NULL;
    const char* outputFileM = NULL;
    const char* inputFile = NULL;
	const char* tmpStr = NULL;
	int systemSize;
    
    //Vi changed the default to rk2
	//const gsl_odeiv2_step_type* steppingFunction = gsl_odeiv2_step_rkck;
    const gsl_odeiv2_step_type* steppingFunction = gsl_odeiv2_step_rk2;
	
	// Set-up model parameters with default arguments
	struct _ModelParameters mParam = {
		.transMembranePermeability = DEFAULT_TRANSMEMBRANE_PERMEABILITY,
		.intracellularVolume = DEFAULT_INTRACELLULAR_VOLUME,
		.targetMoleculeCount = DEFAULT_TARGET_MOLECULE_COUNT,
		.replicationThreshold = DEFAULT_DUMMY, // This will get defaulted based on targetMoleculeCount
		.killingThreshold = DEFAULT_THRESHOLD,   
		.baselineReplication = DEFAULT_BASELINE_REPLICATION,
		.maximumKillRate = DEFAULT_MAXIMUM_KILL_RATE,
		.targetAssociationRate = DEFAULT_TARGET_ASSOCIATION_RATE,
		.targetDissociationRate = DEFAULT_TARGET_DISSOCIATION_RATE,
		.nonSpecificAssociationRate = DEFAULT_NONSPECIFIC_ASSOCIATION_RATE,
		.nonSpecificDissociationRate = DEFAULT_NONSPECIFIC_DISSOCIATION_RATE,
		.timepoints = DEFAULT_TIMEPOINTS,
		.extendedModel = 0,  
                .molecularweight= DEFAULT_MOLECULARWEIGHT,
                .steptime = DEFAULT_STEPTIME,
                .carryingCapacity = DEFAULT_CARRYING_CAPACITY,
		.staticAntibioticConcentration = (double) DEFAULT_DUMMY,
		.hyperGeometricMatrix = NULL,
	};
    
    
    
	// Set-up simulation parameters with default arguments	
	struct _SimulationParameters sParam = {
		.startingAntibiotic = DEFAULT_STARTING_ANTIBIOTIC,
		.startingPopulation = DEFAULT_STARTING_POPULATION,
		.endTime = DEFAULT_SIMULATION_END_TIME,
		.stepSize = DEFAULT_SIMULATION_STEP_SIZE
	};
	
	struct _SimulationResults results;
	
	// Set-up command-line parameters for the carg_parse module
	const struct ap_Option options[] = {
		{ 'v', "verbose",                 ap_no  },
		{ 'h', "help",                    ap_no  },
		{ 'V', "intracellularVolume",     ap_yes },
		{ 'n', "targetMoleculeCount",     ap_yes },
		{ 'r', "replicationThreshold",    ap_yes },
		{ 'k', "killingThreshold",        ap_yes },
		{ 'R', "baselineReplicationRate", ap_yes },
		{ 'K', "maximumKillingRate",      ap_yes },
		{ 'A', "targetAssociationRate",   ap_yes },
		{ 'D', "targetDissociationRate",  ap_yes },
		{ 'C', "carryingCapacity",        ap_yes },
		{ 't', "time",                    ap_yes },
		{ 'd', "startingAntibiotic",      ap_yes },
                { 'M', "Antibiotic Mol. Weight",  ap_yes },
		{ 'p', "startingPopulation",      ap_yes },
		{ 'o', "outputFile",              ap_yes },
                { 'm', "outputFileM",             ap_yes },
                { 'i', "inputFile",               ap_yes },
		{ 'S', "steppingFunction",        ap_yes }
	};
	
	// Grab the invocation name from the command-line
	programName = basename(argv[0]);
	invocationName = argv[0];
	
	// Initialize the carg_parser module with command-line argument array
	if(!ap_init(&parser, argc, (const char *const *)argv, options, 0)) {
		argParserShowError("Not enough memory.", 0, 0);
		return EXIT_FAILURE;
	}
	
	// Check that there are no errors in the supplied arguments
	if(ap_error(&parser)) {
		argParserShowError(ap_error(&parser), 0, 1);
		return EXIT_FAILURE;
	}
	
	// Process options supplied on command-line, into machine variables
	for(argIdx = 0; argIdx < ap_arguments(&parser); ++argIdx) {
		const int code = ap_code(&parser, argIdx);
		if (!code)
			break;
		switch (code) {
		case 'v':
			verbose = 1;
			break;
		case 'h':
			displayHelp(programName);
			return EXIT_SUCCESS;
		case 'V':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.intracellularVolume);
			break;
		case 'n':
			sscanf(ap_argument(&parser, argIdx), "%d", &mParam.targetMoleculeCount);
			break;
		case 'r':
			sscanf(ap_argument(&parser, argIdx), "%d", &mParam.replicationThreshold);
			break;
		case 'k':
			sscanf(ap_argument(&parser, argIdx), "%d", &mParam.killingThreshold);
			break;
		case 'R':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.baselineReplication);
			break;
		case 'K':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.maximumKillRate);
			break;
		case 'A':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.targetAssociationRate);
			break;
		case 'D':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.targetDissociationRate);
			break;
		case 'C':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.carryingCapacity);
			break;
		case 'd':
			sscanf(ap_argument(&parser, argIdx), "%lg", &sParam.startingAntibiotic);
			break;
		case 'M':
			sscanf(ap_argument(&parser, argIdx), "%lg", &mParam.molecularweight);
			break;
		case 'p':
			sscanf(ap_argument(&parser, argIdx), "%lg", &sParam.startingPopulation);
			break;
		case 't':
			sscanf(ap_argument(&parser, argIdx), "%lg:%lg", &sParam.endTime, &sParam.stepSize);
			break;
		case 'o':
			outputFile = ap_argument(&parser, argIdx);
			break;
        case 'm':
            outputFileM = ap_argument(&parser, argIdx);
            break;
        case 'i':
            inputFile = ap_argument(&parser, argIdx);
            break;
        case 'S':
			tmpStr = ap_argument(&parser, argIdx);
			if (!strcmp(tmpStr, "rk4"))
				steppingFunction = gsl_odeiv2_step_rk4;
			else if (!strcmp(tmpStr, "rkf45"))
				steppingFunction = gsl_odeiv2_step_rkf45;
			else if (!strcmp(tmpStr, "rkck"))
				steppingFunction = gsl_odeiv2_step_rkck;
            else if (!strcmp(tmpStr, "msbdf"))//Vi added
				steppingFunction = gsl_odeiv2_step_msbdf;
            else if (!strcmp(tmpStr, "rk2"))
				steppingFunction = gsl_odeiv2_step_rk2;
            else if (!strcmp(tmpStr, "bsimp"))
				steppingFunction = gsl_odeiv2_step_bsimp;
            else if (!strcmp(tmpStr, "msadams"))
				steppingFunction = gsl_odeiv2_step_msadams;
			break;
		default:
			argParserInternalError("uncaught option.");
		}
	}
    
    //-------------------------------------------------------------------------
    // Define the total timepoint from the Simulation time and the interval
    mParam.timepoints = ((int)floorl(sParam.endTime /sParam.stepSize));
    //printf("%d\n",mParam.timepoints);
    double realantibioticconc[mParam.timepoints];
    
    mParam.steptime= sParam.stepSize; 

	if (inputFile != NULL) { 
		mParam.extendedModel = 1;
		printf("Running extended model with input file %s",inputFile);
		// Reading Antibiotic Concentartion from "input" file 
		FILE* myFile;
		myFile = fopen(inputFile, "r");
		if (myFile == NULL) {
			fprintf(stderr, "Could not open %s for reading\n", inputFile);
			return EXIT_FAILURE;
		}
		int x;
		
		for (x = 0; x<mParam.timepoints;x++){
			int parse_success;
			if (feof(myFile)) {
				fprintf(stderr,"File ended while attempting to scan entry %d\n",x+1);
				return EXIT_FAILURE;
			}
			parse_success = fscanf(myFile, "%lg", &mParam.realantibioticconc[x]);
			if (parse_success != 1){
				fprintf(stderr,"Parsing of entry %d failed\n",x+1);
				return EXIT_FAILURE;
			}
			// Old code: Changing nG/mL to number of molecules: A*1e3*1e-6*6.02e23/MW    //mParam.realantibioticconc[x]=mParam.realantibioticconc[x]*6.02e17*mParam.intracellularVolume/mParam.molecularweight;
				
			// Newcode by Vi 
			mParam.realantibioticconc[x]=mParam.realantibioticconc[x]*6.02e20*mParam.intracellularVolume/mParam.molecularweight;

		}
		fclose(myFile);
		myFile=NULL;
	} else {
		mParam.staticAntibioticConcentration = sParam.startingAntibiotic*6.02e20*mParam.intracellularVolume/mParam.molecularweight;
	}
    
    //-------------------------------------------------------------------------
    
	// More defaults; thresholds are set to half the total number of target molecules if no explicit
	// threshold was provided on the command-line
	/*if (mParam.replicationThreshold == DEFAULT_DUMMY)
        mParam.replicationThreshold = mParam.targetMoleculeCount / 2;
	if (mParam.killingThreshold == DEFAULT_DUMMY)
		mParam.killingThreshold = mParam.targetMoleculeCount / 2;*/
    //Dec. 12, Vi changed to, the first case is for the webpage where only killingThreshold is determined.
	if ((mParam.killingThreshold != DEFAULT_DUMMY)&&(mParam.replicationThreshold == DEFAULT_DUMMY))
        mParam.replicationThreshold = mParam.killingThreshold -1;
    else 
    {
        if (mParam.replicationThreshold == DEFAULT_DUMMY)
        mParam.replicationThreshold = mParam.targetMoleculeCount / 2;
	    if (mParam.killingThreshold == DEFAULT_DUMMY)
		mParam.killingThreshold = mParam.targetMoleculeCount / 2;
    }
	
    
	// Output the run parameters if verbose-mode is specified
	if (verbose) {
        printf("\n");
        printf("TB-Simulation: Version 6.0 \n");
        printf("------------------------------------\n\n");
		printf("\nStarting simulation with arguments\n");
		printf("------------------------------------\n\n");
		printf("Starting population     \t%lg\n", sParam.startingPopulation);
		if (!mParam.extendedModel) {
			printf("Antibiotic dose         \t%lg\n", sParam.startingAntibiotic);
		}
		printf("Time of simulation      \t%lg\n", sParam.endTime);
		printf("Step size               \t%lg\n\n", sParam.stepSize);		
		printf("Target molecules        \t%d\n", mParam.targetMoleculeCount);
		printf("Maximum kill rate       \t%lg\n", mParam.maximumKillRate);
		printf("Killing threshold       \t%d\n",  mParam.killingThreshold);
        printf("Replication threshold   \t%d\n",  mParam.replicationThreshold);
		printf("Baseline replication    \t%lg\n", mParam.baselineReplication);
		printf("Target association rate \t%lg\n", mParam.targetAssociationRate);
		printf("Target dissociation rate\t%lg\n", mParam.targetDissociationRate);
		printf("Drug Molecular Weight    \t%lg\n", mParam.molecularweight);
		printf("Carrying capacity       \t%lg\n", mParam.carryingCapacity);
		printf("Intracellular volume    \t%lg\n", mParam.intracellularVolume);
	}
	
	// Make sure that no weird parameters have been supplied
	if (sanityCheckModelParameters(&mParam) != 0) {
		fprintf(stderr, "Failure: Bad parameters supplied\n");
		return EXIT_FAILURE;
	}
    
    FILE* oHandleM;
    if (outputFileM != NULL) {        
        if (verbose)
            printf("Outputing compartmentBoundComplexState matrix to %s\n", outputFileM);
            printf("Output File order is as follows:\n");
        if ((oHandleM = fopen(outputFileM, "wb")) == NULL) {
            fprintf(stderr, "Could not open %s for writing\n", outputFileM);
            return EXIT_FAILURE;
        }
    } else{
        
    }
    
    // creating header for the output file
    
    int leng=mParam.targetMoleculeCount+4;
    const char *strs[leng];
            for (i = 1; i <= mParam.targetMoleculeCount; ++i){
            strs[i]="Li ";
        }
        strs[0]= "L0 ";
        strs[mParam.targetMoleculeCount]= "Ln ";
        strs[mParam.targetMoleculeCount+1]= "tm ";
        strs[mParam.targetMoleculeCount+2]= "BP ";
        strs[mParam.targetMoleculeCount+3]= "An ";
        strs[mParam.targetMoleculeCount+4]= "AT ";
    
    char headout[leng*3];
        strcpy(headout, strs[0]);
   for (i = 1; i <= leng; ++i){
        strcat(headout, strs[i]);
        }
   printf("%s\n",headout);
  
	// Initialize the initial state to all bacteria without bound targets and the initial dose in the extracellular medium
	stateVector = initializeStateVector(mParam.targetMoleculeCount, sParam.startingAntibiotic, sParam.startingPopulation);
	
	// Run the simulation itself, and measure its execution time
	t = clock();
	mParam.hyperGeometricMatrix = generateHypergeometricMatrix(mParam.targetMoleculeCount, mParam.replicationThreshold);
	if (runSimulation(steppingFunction, &mParam,  sParam.endTime, sParam.stepSize, stateVector, &results, outputFileM, oHandleM) != GSL_SUCCESS) {
		fprintf(stderr, "The simulation failed.\n");
		return EXIT_FAILURE;
	}
	t = clock() - t;
    
    // write the header to the output file
    fprintf (oHandleM, "%s\n", headout);
    fclose(oHandleM);
	
	// Output a summary of results, if verbose-mode is specified
	if (verbose) {
		double populationSum = 0.0;
		for (i = NUMBER_FREE_KINETIC_VARIABLES; i < mParam.targetMoleculeCount+NUMBER_FREE_KINETIC_VARIABLES + 1; ++i)
			populationSum += stateVector[i];
		printf("Results readout\n");
		printf("---------------\n\n");
		printf("Final population %g\n\n",populationSum);
		printf("It took me (%f milliseconds).\n\n",((float)t*1000.0)/CLOCKS_PER_SEC);
	}
	
	if (outputFile != NULL) {
		FILE* oHandle;
		
		if (verbose)
			printf("Outputing results to %s\n", outputFile);
		if ((oHandle = fopen(outputFile, "wb")) == NULL) {
			fprintf(stderr, "Could not open %s for writing\n", outputFile);
			return EXIT_FAILURE;
		}
		if (outputResultsToFile(&sParam, &mParam, &results, oHandle) == -1) {
			fprintf(stderr, "The was an error writing the output file\n");
			return EXIT_FAILURE;
		}
	}
    
	return EXIT_SUCCESS;
}


/**
 * Called if there is an error detected on the command-line such as a mismatched option
 * or the wrong number of options given.
 * 
 * @param msg       error message to display
 * @param errorCode numerical code for the error, if present
 * @param help      indicates whether help is available to provide hint to user
 */
static void argParserShowError(const char* const msg, const int errorCode, const char help) {
	if (msg && msg[0]) {
		fprintf(stderr, "%s: %s", programName, msg);
		if (errorCode > 0)
			fprintf(stderr, ": %s.", strerror(errorCode));
		fprintf(stderr, "\n");
	}
	if(help)
		fprintf(stderr, "Try '%s --help' for more information.\n", invocationName);
}

/**
 * Called during argument parsing if there was a problem in the arg_parser module.
 *
 * @param msg  message to display on exit
 */
static void argParserInternalError(const char* const msg) {
	fprintf(stderr, "%s: internal error: %s\n", programName, msg );
	exit(EXIT_FAILURE);
}

/**
 * Collect the option supplied to a command-line argument.
 *
 * @param code     Argument code
 * @param options  Array defined for the argument syntax of the program
 *
 * @return         The option as a string
 */
static const char* optname(const int code, const struct ap_Option options[]) {
	static char buf[2] = "?";
	int i;
	
	if( code != 0 )
		for(i = 0; options[i].code; ++i)
			if(code == options[i].code) {
				if(options[i].name)
					return options[i].name;
				else break;
			}
	
	if(code > 0 && code < 256)
		buf[0] = code;
	else
		buf[0] = '?';
	
	return buf;
}

//TODO ADD UNITS TODO
/**
 * Display the help annotation for the --help command-line option.
 *
 * @param programName  The invoked name of this program for display purposes
 */
static void displayHelp(const char* programName) {
	printf("\n"
	       "  %s - Tuberculosis simulation software Version 4.2b %s (2016)\n\n"
	       "  usage: %s [options]\n\n",
	       programName, VERSION_STRING, programName);
	
	printf("                                 GENERAL OPTIONS\n\n"
	       "   -h, --help    : Display this help message.\n"
	       "   -v, --verbose : Display extra information during program run.\n\n");
	
	printf("                              SIMULATION PARAMETERS\n\n"
	       "   -d, --startingAntibiotic [dose]   : Initial dose quantity.\n"
	       "                                         default: %lg\n"
	       "   -p, --startingPopulation [population]    : Initial bacterial population.\n"
	       "                                         default: %lg\n"
	       "   -S, --steppingFunction [function] : Stepping function to use for the numerical integration.\n"
	       "                                         where [function] is one of {rk4, rkf45, rkck, rk2}\n"
	       "                                         default: rk2\n"
	       "   -t, --time [etime (s)]:[intvl (s)]   : Specifies total simulation time [etime] and interval between time-points [intvl].\n"
	       "                                         default: %lg:%lg\n\n",
	       DEFAULT_STARTING_ANTIBIOTIC, DEFAULT_STARTING_POPULATION, DEFAULT_SIMULATION_END_TIME, DEFAULT_SIMULATION_STEP_SIZE);
	
	printf("                                 MODEL PARAMETERS\n\n"
	       "   -n, --targetMoleculeCount [Integer Number]     : Number of target molecules in a cell.\n"
	       "                                            default: %d\n"
		   "   -r, --replicationThreshold [FC*n=Integer Number]     : Number of bound target molecules in a cell to stop relications.\n"
	       "                                            default: n / 2\n"
		   "   -R, --baselineReplicationRate [rate (1/s)] : Rate of replication for those bacteria below replication threshold.\n"
	       "                                            default: %lg\n"
		   
	       "   -k, --killingThreshold [FC*n=Integer Number]     : Number of bound target molecules in a cell to cause death.\n"
	       "                                            default: n / 2\n"
	       "   -K, --maximumKillingRate [rate (1/s)]      : Rate of death for those bacteria above killing threshold.\n"
	       "                                            default: %lg\n"
           "   -M, --molecularweight [weight (gr/mol)]         : Drug Molecular Weight gram per mole.\n"
	       "                                            default: %lg\n" 
	       "   -A, --targetAssociationRate [rate (L/mol/s)]   : Rate constant for association between target and antibiotic.\n"
	       "                                            default: %lg\n"
	       "   -D, --targetDissociationRate [rate (1/s)]  : Rate constant for dissociation of target/antibiotic complex.\n"
	       "                                            default: %lg\n"
	       "   -V, --intracellularVolume [size (L)]     : Internal volume of a bacterium.\n"
	       "                                            default: %lg\n"
	       "   -C, --carryingCapacity [population]         : Carrying capacity (maximum population) of the system.\n"
	       "                                            default: %lg\n\n",
	       DEFAULT_TARGET_MOLECULE_COUNT, DEFAULT_BASELINE_REPLICATION, DEFAULT_MAXIMUM_KILL_RATE, DEFAULT_MOLECULARWEIGHT, DEFAULT_TARGET_ASSOCIATION_RATE, DEFAULT_TARGET_DISSOCIATION_RATE, DEFAULT_INTRACELLULAR_VOLUME, DEFAULT_CARRYING_CAPACITY);
	
	printf("                                DATA OUTPUT OPTIONS\n\n"
           "   -i, --inputFile [ofile]   : Read Drug Concentration from [ofile].\n\n"
           "   -m, --outputFileM [ofile] : Write intracellular compartment vectors to [ofile].\n\n");

}

#define EMIT_YAML_EVENT if (!yaml_emitter_emit(&emitter, &event)) goto error

static int writeNamedDouble(yaml_emitter_t* emitter, yaml_event_t* event, char* name, double value);

static int outputResultsToFile(SimulationParameters sParam, ModelParameters mParam, SimulationResults results, FILE* oHandle) {
	yaml_emitter_t emitter;
	yaml_event_t event;
	
	// Intialize the YAML emitter
	yaml_emitter_initialize(&emitter);
	yaml_emitter_set_output_file(&emitter, oHandle);
	yaml_emitter_set_unicode(&emitter, 1);
	
	// Intialize the YAML stream
	yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
	EMIT_YAML_EVENT;
	
	// Initialize the YAML document
	yaml_version_directive_t YAMLVersion = {
		.major = 1,
		.minor = 1
	};
	yaml_tag_directive_t versionTag = {
		.handle = (yaml_char_t*)("tuberculosis-simulator," VERSION_STRING),
		.prefix = (yaml_char_t*)"!"
	};

	yaml_document_start_event_initialize(&event, &YAMLVersion, &versionTag, &versionTag, 0);
	EMIT_YAML_EVENT;
	
	// Begin the main map, consisting of simulation parameters, model parameters, time-output header, time-output list and final-output
	yaml_mapping_start_event_initialize(&event, NULL, NULL, 1, YAML_BLOCK_MAPPING_STYLE);
	EMIT_YAML_EVENT;
	
	// TODO SIMULATION PARAMETERS
	
	// Write the mapping header for the simulation paramaters
	yaml_scalar_event_initialize(&event, NULL, NULL, (yaml_char_t*)"simulation-parameters", -1, 1, 1, YAML_PLAIN_SCALAR_STYLE);
	EMIT_YAML_EVENT;
	
	// Begin the map for the simulation parameters
	yaml_mapping_start_event_initialize(&event, NULL, NULL, 1, YAML_BLOCK_MAPPING_STYLE);
	EMIT_YAML_EVENT;
	
	if (writeNamedDouble(&emitter, &event, "starting-population", sParam->startingPopulation) == -1)
		goto error;
	if (writeNamedDouble(&emitter, &event, "starting-antibiotic", sParam->startingAntibiotic) == -1)
		goto error;
	
	// End the map for the simulation parameters
	yaml_mapping_end_event_initialize(&event);
	EMIT_YAML_EVENT;
	
	// End the main map
	yaml_mapping_end_event_initialize(&event);
	EMIT_YAML_EVENT;
	
	// End the document
	yaml_document_end_event_initialize(&event, 1);
	EMIT_YAML_EVENT;
	
	// End the stream
	yaml_stream_end_event_initialize(&event);
	EMIT_YAML_EVENT;

success:
	yaml_emitter_delete(&emitter);
	return 0;
error:
	yaml_emitter_delete(&emitter);
	return -1;
}

static int writeNamedDouble(yaml_emitter_t* emitter, yaml_event_t* event, char* name, double value) {
	char buffer[255];
	sprintf(buffer, "%lg", value);
	
	yaml_scalar_event_initialize(event, NULL, NULL, (yaml_char_t*)name, -1, 1, 1, YAML_PLAIN_SCALAR_STYLE);
	if (!yaml_emitter_emit(emitter, event))
		goto error;
	
	yaml_scalar_event_initialize(event, NULL, NULL, (yaml_char_t*)buffer, -1, 1, 1, YAML_PLAIN_SCALAR_STYLE);
	if (!yaml_emitter_emit(emitter, event))
		goto error;
	
success:
	return 0;
error:
	return -1;
}
