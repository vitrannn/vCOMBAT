/**
 * @file   full_model.h
 * @version 5
 * @revision Vi Tran (original version written by T on 24 April 2015, the Alireza in Oct 2016) 
 * @updated  Feb-2020
 * @brief  Default model arguments and definitions for full_model.c
 */

#define AVOGADRO_CONSTANT (6.02e23) ///< Definition of Avogadro's number for calculation of derivatives

/* Default options to be over-written by command-line
*/
#define DEFAULT_TRANSMEMBRANE_PERMEABILITY 0.0   ///< NOT IMPLEMENTED YET
//#define DEFAULT_INTRACELLULAR_VOLUME 3.0e-16        ///< For M.Tuberculosis
#define DEFAULT_INTRACELLULAR_VOLUME 1.0e-15        /// Vi changed to 1.0e-15
#define DEFAULT_TARGET_MOLECULE_COUNT 100         ///< From equationparserv2.R
#define DEFAULT_BASELINE_REPLICATION 8.34e-6         ///< From equationparserv2.R
#define DEFAULT_MAXIMUM_KILL_RATE 1.39e-5           ///< From equationparserv2.R
#define DEFAULT_MOLECULARWEIGHT 555.5  
#define DEFAULT_TARGET_ASSOCIATION_RATE 4450.0      ///< From equationparserv2.R
#define DEFAULT_TARGET_DISSOCIATION_RATE 0.0023   ///< From equationparserv2.R
#define DEFAULT_NONSPECIFIC_ASSOCIATION_RATE 0   ///< NOT IMPLEMENTED YET
#define DEFAULT_NONSPECIFIC_DISSOCIATION_RATE 0  ///< NOT IMPLEMENTED YET
#define DEFAULT_CARRYING_CAPACITY 1e9            ///< From equationparserv2.R
#define DEFAULT_TIMEPOINTS 360000.0
#define DEFAULT_STEPTIME 3600.0
#define DEFAULT_THRESHOLD 60
#define DEAULT_STARTING_ANTIBIOTIC 1.0
  


/**
 * Structure to hold all the per-simulation parameters of the model these will not be updated during the simulation but set
 * before each simulation begins. They will be passed to the derivative function in the "parameters" pointer.
 */
typedef struct  _ModelParameters {
	/* TODO : Add units for each to the comments
	** Cellular details.
	*/
	double transMembranePermeability;   ///< Permeability rate for the bacterial membrane.
	double intracellularVolume;         ///< Volume of each cell.
	
	/* Parameters related to target molecule effects on growth/death
	*/
	int targetMoleculeCount;            ///< Total number of target molecules per bacterial cell.
	int replicationThreshold;           ///< Threshold for simple step function, after this many targets are bound replication will be reduced to zero.
	int killingThreshold;               ///< Threshold for simple step function, after this many targets are bound killing will be raised to maximumKillRate.
    int timepoints;
    double steptime;
	double baselineReplication;         ///< The replication rate for cells that have not reached the replicationThreshold in  terms of bound targets.
	double maximumKillRate;             ///< The killing rate for cells that have passed the killingThrehsold in terms of bound  targets.
	
	/* Parameters related to anti-biotic kinetics
	*/
        double molecularweight;             ///< Molecular Weight of Antibiotic.
	double nonSpecificAssociationRate;  ///< The k parameter for the forward reaction of non-specific binding.
	double nonSpecificDissociationRate; ///< The k parameter for the backward reaction of non-specific binding.
	double targetAssociationRate;       ///< The k parameter for the forward reaction of specific binding.
	double targetDissociationRate;      ///< the k parameter for the backward reaction of specific binding.
	int extendedModel; // Specifies whether the program runs with an input file (extended model) or with a static total antibotic concentration (original model)
	
	double carryingCapacity;            ///< The total carrying capacity (maximum population) of the system.
	double* hyperGeometricMatrix;       ///< The matrix containing the hypergeometric sampling PDFs.
	double staticAntibioticConcentration; // In the case we are using the original model instead of the extended model
    double realantibioticconc[]; // In number of molecules, must be at the end of the struct because it is flexible. Note that this array is uninitialized and its contents are undefined in case we are using the original model
    //double *realantibioticconc; 
} ModelParameters;

/**
 * Used for type-casting the derivative function (which uses structures rather than vectors for lexical ease) into the format expected by GSL.
 */
typedef int (*GSLDerivCalcFunc)(double, const double*, double*, void*);

/**
 * Structure to hold all the within-simulation variables. This will be used both for the current-state and also in the calculation
 * of the derivatives in the GSL sub-function. This structure is for the deterministic model so uses concentrations of molecules.
 * This structure can both be interpreted as an array of doubles, but gives us the convenience of being able to refer to
 * each field by name
 */
typedef struct _ModelVariables {
	/* TODO : Add units for each to the comments
	*/
	double freeAntibiotic;               ///< Concentration of free (extracellular) unbound antibiotic.
	double freeTarget;                   ///< Concentration of free (extracellular) unbound target molecule.
	double freeBoundComplex;             ///< Concentration of free (extracellular) bound target/anti-biotic complex.
	double CompartmentBoundComplex[]; ///< Total number of cells in each compartment for intracellular bound complex. [GIVEN IN MOLECULES]
} ModelVariables;

///> Macro extracts the number of "free" compartments. I.e. those concentrations which are not in the array of cells-with-bound-targets
#define NUMBER_FREE_KINETIC_VARIABLES 3

int calculateModelDerivative_BindingOnly (double curTime,
                                          const ModelVariables* y,
                                          ModelVariables* dydt,
                                          ModelParameters* param);

int sanityCheckModelParameters(ModelParameters* param);
