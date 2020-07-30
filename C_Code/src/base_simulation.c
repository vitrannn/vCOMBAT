/**
 * @file   base_simulation.c
 * @version 5
 * @revision Vi Tran (original version written by T on 24 April 2015) 
 * @updated  2020
 * @brief  Simulation initialization and main-loop procedures
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include "carg_parser.h"
#include "full_model.h"
#include "base_simulation.h"

extern int verbose;

/**
 * Simulation result aggregator function. Calculates the summary statistics which are to be calculated for each time-point and
 * stores them in the results structure, along with the temporal coordinate.
 *
 * @param mParam      Model parameters for the simulation.
 * @param state       Current state of the simulation (i.e. Structure of the dependent variables).
 * @param currentTime The absolute temporal coordinate for the current state.
 * @param counter     The number of time-points that have elapsed since the simulation started.
 * @param results     The results structure to be updated with the current simulation summary.
 */
static void updateSimulationResultsPerTick(const ModelParameters mParam, const ModelVariables state, const double currentTime,
                                           const int counter, SimulationResults results, const char* output, FILE* oHandleM) {
	int i;
	double populationSum = 0.0;
    int timetocon;
    timetocon=((int)floorl(currentTime/mParam->steptime));
    double plasmaconcenc = mParam->realantibioticconc[timetocon];
	double* compartmentBoundComplexState = &state->firstCompartmentBoundComplex;
    //double* compartmentfreeAntibiotic = &state->freeAntibiotic;
	double* compartmentfreeBoundComplex= &state->freeBoundComplex;
    double* compartmentfreetarget= &state->freeTarget;
    
	for (i = 0; i < mParam->targetMoleculeCount + 1; ++i)
		populationSum += compartmentBoundComplexState[i];
    if (populationSum < 0) {
		populationSum=0;
    }
//	printf("CompartmentBoundComplexState vector at time %lf is:\n", currentTime);
	if(output != NULL){
		for (i = 0; i < mParam->targetMoleculeCount + 1; ++i)
        //fprintf(oHandleM,"L \t%d\n", i);
        fprintf(oHandleM, "%lf ", compartmentBoundComplexState[i]);
        fprintf(oHandleM, "%lf ", currentTime);
        fprintf(oHandleM, "%lf ", populationSum);
        fprintf(oHandleM, "%lf ", mParam->realantibioticconc[timetocon]);
        fprintf(oHandleM, "%lf ", compartmentfreeBoundComplex[i]);
		fprintf(oHandleM, "\n");
        
	}
	
	results->timePoint[counter] = currentTime;
	results->totalPopulation[counter] = populationSum;
    results->unboundantibiotic[counter]= plasmaconcenc;
	
	if (verbose) {
		printf("%.4lg %.8lg %.8lg                   \r", currentTime, populationSum,plasmaconcenc);
		fflush(stdout);
	}
}

/**
 * The main simulation loop function. Will set up the ODE system with GSL and run the simulation within the specified
 * time bounds. Will dump the output to the specified file.
 *
 * @param stepping      GSL stepping function to use for the ODE solver.
 * @param mParam        Model parameters for the simulation.
 * @param endTime       The time to run the simulation until.
 * @param timeInterval  The amount of time between data-points.
 * @param stateVector   Initial starting conditions as input and final conditions as output.
 *
 * @return              GSL_SUCCESS if everything went well otherwise the GSL error code. 
 */
int runSimulation(const gsl_odeiv2_step_type* stepping, const ModelParameters mParam, const double endTime,
                  const double timeInterval, double* stateVector, SimulationResults results, const char* output, FILE* oHandleM) {
	int curTimePoint = 0;
	int systemSize = NUMBER_FREE_KINETIC_VARIABLES + mParam->targetMoleculeCount + 1;
	double nextTime = timeInterval;
	double curTime = 0.0;
	int measureRuntime = 0;
	int totalTimePoints = ((int)floorl(endTime / timeInterval)) + 1.0;
	
    
	results->timePoint = malloc(sizeof(double) * totalTimePoints);
	results->totalPopulation = malloc(sizeof(double) * totalTimePoints);
    results->unboundantibiotic = malloc(sizeof(double) * totalTimePoints);
	
	if (verbose)
		printf("\ncreating system with %d free variables\n", NUMBER_FREE_KINETIC_VARIABLES + mParam->targetMoleculeCount + 1);
	
	gsl_odeiv2_system sys = {(GSLDerivCalcFunc)calculateModelDerivative_BindingOnly, NULL, systemSize, mParam};
	gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new (&sys, stepping, timeInterval, 1e-5, 1e-5);
	
	updateSimulationResultsPerTick(mParam, (ModelVariables)stateVector, curTime, curTimePoint, results, output, oHandleM);
	
	while (nextTime < endTime) {
		int i;
		int status;
		
		++curTimePoint;
		
		if ((status = gsl_odeiv2_driver_apply (driver, &curTime, nextTime, stateVector)) != GSL_SUCCESS) {
			fprintf (stderr, "error in  gsl_odeiv2_driver_apply: %d (%s)\n", status, gsl_strerror (status));
			return status;
		}
		updateSimulationResultsPerTick(mParam, (ModelVariables)stateVector, curTime, curTimePoint, results, output, oHandleM);
		
		nextTime += timeInterval;
	}
	if (verbose)
		printf("\n\n");

	return GSL_SUCCESS;
}

/**
 * Function to pre-generate the hyper-geometric coefficients for the subsequent simulation. The entire matrix does not need
 * to be generated if a full step function is assumed, because reproduction will cease completely after the reproduction
 * threshold. Therefore the matrix is only calculated, one-sided, for those values below this threshold.
 *
 * @param populationCount       Total number of target molecules per cell.
 * @param replicationThreshold  Threshold at which replication ceases.
 *
 * @return                      One side of hypergeometric matrix (including diagonal) compressed into a linear vector.
 */
double* generateHypergeometricMatrix(const int populationCount, const int replicationThreshold) {
	int i,j;
	int matrixSize = (replicationThreshold + 1) * replicationThreshold / 2;
	double* matrix = (double*)malloc(matrixSize * sizeof(double));
	double* mPointer = matrix;
	
	double staticChoose = gsl_sf_lnchoose(2 * populationCount, populationCount);
	for (i = 0; i < replicationThreshold; ++i)
		for (j = i; j < replicationThreshold; ++j) {
			double tmpChoose1 = gsl_sf_lnchoose(j, i);
			double tmpChoose2 = gsl_sf_lnchoose(2 * populationCount - j, populationCount - i);
			
            /*double arg = tmpChoose1 + tmpChoose2 - staticChoose;
			double expOfArg = 0;
			if(arg > -700.) expOfArg = gsl_sf_exp(tmpChoose1 + tmpChoose2 - staticChoose);
			*mPointer++ = expOfArg;*/
            //Vi changes back to Trevor version
            *mPointer++ = gsl_sf_exp(tmpChoose1 + tmpChoose2 - staticChoose);
		}
		
	return matrix;
}

/**
 * Function to initialize the initial state of the simulation with a given population and quantity of antibiotic.
 *
 * @param targetMoleculeCount  Number of traget molecules per cell
 * @param startingAntibiotic   Starting dose of anti-biotic in the extracelullar medium
 * @param startingPopulation   Initial population of bacteria with no bound targets.
 *
 * @return                     The initialized state vector with antibiotic dose and zero-bound compartments filled with others zeroed.
 */
double* initializeStateVector(const int targetMoleculeCount, const double startingAntibiotic, const double startingPopulation) {
	int i;
	int systemSize = NUMBER_FREE_KINETIC_VARIABLES + targetMoleculeCount + 1;
	double* stateVector = (double*)malloc(sizeof(double) * systemSize);
	
	for (i=0; i< systemSize; i++)
		stateVector[i] = 0.0;
	stateVector[0] = startingAntibiotic;
	stateVector[NUMBER_FREE_KINETIC_VARIABLES] = startingPopulation;
	
	return stateVector;
}

