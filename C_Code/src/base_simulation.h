/**
 * @file   base_simulation.h
 * @version 5
 * @revision Vi Tran (original version written by T on 24 April 2015) 
 * @updated  2020
 * @brief  Default simulation arguments and definitions for base_simulation.c
 */

#define DEFAULT_STARTING_ANTIBIOTIC 1e4  ///< From equationparserv2.R (for 700 microgram BDQ)
#define DEFAULT_STARTING_POPULATION 1e6    ///< From equationparserv2.R
#define DEFAULT_SIMULATION_END_TIME 360000.0  ///< From equationparserv2.R
#define DEFAULT_SIMULATION_STEP_SIZE 3600.0   ///< From equationparserv2.R


/**
 * Structure to hold all the parameters relating to how the model is simulated
 */
typedef struct _SimulationParameters {
	double startingAntibiotic;  ///< The initial dose of antibiotic.
	double startingPopulation; ///< The initial bacterial population of the system.
	double endTime;             ///< The time to run the simultion until for a SIMULATION_TYPE_FIXED_TIME
	double stepSize;            ///< The time-delta between time-points
} *SimulationParameters;

/**
 * Structure to hold the aggregated results of a simulation
 */
typedef struct _SimulationResults {
	double* timePoint;       ///< Vector containing the list of time-points measurements are taken at.
	double* totalPopulation; ///< Vector containing the list of population counts for each time-point.
    double* unboundantibiotic; ///< Vector containing the list of free antibiotic concentartion for each time-point.
	double finalTime;        ///< The final time-point of the system.
	double finalPopulation;  ///< The final population count of the system.
} *SimulationResults;

int runSimulation(const gsl_odeiv2_step_type* stepping, const ModelParameters mParam, const double endTime,
                  const double timeInterval, double* stateVector, SimulationResults results, const char* output, FILE* oHandleM);

double* generateHypergeometricMatrix(const int populationCount, const int replicationThreshold);

double* initializeStateVector(const int targetMoleculeCount, const double startingAntibiotic, const double startingPopulation);
