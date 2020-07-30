/**
 * @file   full_model.c
 * @version 5
 * @revision Vi Tran (original version written by T on 24 April 2015, then Alireza in Nov 2016) 
 * @updated  Feb-2020
 * @brief  Program entry, command-line processing and simulation invocation
 * @comments  Program entry, command-line processing and simulation invocation
 */

#include "full_model.h"
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <math.h>

/**
 * gsl_odeiv2_system inner function for calculating the derivative of the Bacteriostatic and Bactericidal action model
 * for a deterministic (concentration-based) simulation.
 *
 * @param curTime  The current time-point, provided by the GSL outer ODE function.
 * @param y        The input vector of state variables, interpreted as a structure for lexical ease.
 * @param dydt     The output vector of state derivatives, interpreted as a structure for lexical ease.
 * @param param    The model parameters. These do not change throughout the simulation.
 *
 * @return         GSL_SUCCESS on success. Failure not currently detected.
 */
int calculateModelDerivative_BindingOnly (double curTime,
                                          ModelVariables y,
                                          ModelVariables dydt,
                                          ModelParameters param) {
	int i,j,x;
	
	// Extract the intracellular compartment vectors from the state and the derivative structures
	double* compartmentBoundComplexState = &y->firstCompartmentBoundComplex;
	double* compartmentBoundComplexDeriv = &dydt->firstCompartmentBoundComplex;
	
	// TODO: This can be pre-calculated
	// Scratch space for calculation of $\frac{k_f}{n_AV_i}$
	double scratchVolumeModifiedK;
    
    // array of drug concentration in patients' body

     /* receiving inputs from user and storing it in array*/
 
    /* receiving inputs from user and storing it in array*/
    //for (x=0; x<=param->timepoints;x++)
   // {
       // printf("enter Antibiotic concentration %d\n", x);
      //  scanf("%lg", &realantibioticconc[x]);
   // }
	
	// Scratch space for calculation of $\frac{k_f}{n_AV_i}A \times B_x$
	double scratchForwardRateComponent[param->targetMoleculeCount];
	
	// Scratch space for calculation of $k_rxB_x$
	double scratchBackwardRateComponent[param->targetMoleculeCount];
	
	// Scratch space for calculation of $\frac{K - \sum_{j=0}^nB_j}{K}$
	double scratchReplicationSum = compartmentBoundComplexState[0];
	
	// Scratch space for calculation of sums for $\frac{dA}{dt}$
	double scratchSumForward = 0.0;
	double scratchSumReverse = 0.0;
	double scratchAntibioticTarget;
	double scratchSumDeath1 = 0.0;
	double scratchSumDeath2 = 0.0;
	double tmp;
	double tmpSum;
	double* incPointer;
	
     //Vi moved from line 150
    int timetocon;
    double yfreeAntibiotic;
    timetocon=((int)floorl(curTime/param->steptime));
    
    // Free Antibiotic Concentrations come from input file, change from yfreeAntibiotic to yfreeAntibiotic by declaring a new code variable
	yfreeAntibiotic = (curTime-timetocon*param->steptime)*(param->realantibioticconc[timetocon+1] - param->realantibioticconc[timetocon])/param->steptime +param->realantibioticconc[timetocon];
    //end change--------------
    
	// Calculation of $\frac{k_f}{n_AV_i}$
	scratchVolumeModifiedK = param->targetAssociationRate / (AVOGADRO_CONSTANT * param->intracellularVolume);
	
	// Calculation of $(\frac{k_f}{n_AV_i}A \times B_x for 0 < x < n$
	for (i = 0, j = 1; i < param->targetMoleculeCount; ++i, ++j) {
		scratchForwardRateComponent[i] = scratchVolumeModifiedK * yfreeAntibiotic * compartmentBoundComplexState[i];
		scratchBackwardRateComponent[i] =  param->targetDissociationRate * j * compartmentBoundComplexState[j];
		scratchReplicationSum += compartmentBoundComplexState[j];
		scratchSumForward += (param->targetMoleculeCount - i) * scratchForwardRateComponent[i];
		scratchSumReverse += scratchBackwardRateComponent[i];
		if (j >= param->killingThreshold) {
			scratchSumDeath1 +=  compartmentBoundComplexState[j] * (param->targetMoleculeCount - j);
			scratchSumDeath2 +=  compartmentBoundComplexState[j] * j;
		}
	}
	
	// Calculation of $\frac{K - \sum_{j=0}^nB_j}{K}
	scratchReplicationSum = (param->carryingCapacity - scratchReplicationSum) / param->carryingCapacity;
	
	// Calculation of $\frac{dB_0}{dt}$
	incPointer = param->hyperGeometricMatrix;
	tmpSum = 0.0;
	for (i=0; i < param->replicationThreshold; ++i)
		tmpSum += *incPointer++ * compartmentBoundComplexState[i];
	tmpSum = 2.0 * param->baselineReplication * (1.0-i/param->targetMoleculeCount) * tmpSum * scratchReplicationSum
	       - param->baselineReplication * (1.0-i/param->targetMoleculeCount) * scratchReplicationSum * compartmentBoundComplexState[0];
	
	// If the killing threshold is set to zero (unlikely but possible) then death must be factored into
	// bacteria with no bound target molecules
	if (param->killingThreshold == 0) {
		scratchSumDeath1 +=  compartmentBoundComplexState[0] * param->targetMoleculeCount;
		tmpSum -= param->maximumKillRate * compartmentBoundComplexState[0];
	}
	
	// scratchSumDeath can only be calculated after the above statement
	scratchSumDeath1 *= param->maximumKillRate;
	scratchSumDeath2 *= param->maximumKillRate;
	
	compartmentBoundComplexDeriv[0] = scratchBackwardRateComponent[0] 
	                                - param->targetMoleculeCount * scratchForwardRateComponent[0]
	                                + tmpSum;
	
	// Calculation of $\frac{dB_n}{dt}$
	compartmentBoundComplexDeriv[param->targetMoleculeCount] = scratchForwardRateComponent[param->targetMoleculeCount - 1] 
	                                                         - scratchBackwardRateComponent[param->targetMoleculeCount - 1]
	                                                         - param->maximumKillRate * compartmentBoundComplexState[param->targetMoleculeCount];
	
	// Calculation of $\frac{dB_i}{dt]$
	for (i=1; i < param->targetMoleculeCount; i++) {
		// This part calculates those components which depend on whether the compartment in question is below
		// the killing threshold.
		tmpSum = 0.0;
		if (i < param->replicationThreshold) {
			for (j=i; j < param->replicationThreshold; ++j)
				tmpSum += *incPointer++ * compartmentBoundComplexState[j];
			tmpSum = 2.0 * param->baselineReplication * (1.0-i/param->targetMoleculeCount) * tmpSum * scratchReplicationSum
			       - param->baselineReplication * (1.0-i/param->targetMoleculeCount) * scratchReplicationSum * compartmentBoundComplexState[i];
		}
		if (i >= param->killingThreshold)
			tmpSum -= param->maximumKillRate * compartmentBoundComplexState[i];
		
		compartmentBoundComplexDeriv[i] = (param->targetMoleculeCount - i + 1) * scratchForwardRateComponent[i - 1]
		                                - (param->targetMoleculeCount - i) * scratchForwardRateComponent[i]
		                                + scratchBackwardRateComponent[i]
		                                - scratchBackwardRateComponent[i - 1]
		                                + tmpSum;
	}
    
	// Calculation of $\frac{k_f}{n_AV_i}A.T - k_rAT$
	scratchAntibioticTarget = (yfreeAntibiotic * y->freeTarget * scratchVolumeModifiedK) - (param->targetDissociationRate * y->freeBoundComplex);
	
    //int timetocon;
    //timetocon=((int)floorl(curTime/param->steptime));
    // Free Antibiotic Concentrations come from input file
    //y->freeAntibiotic = param->realantibioticconc[timetocon];
    // Calculation of $\frac{dA}{dt}$
	//dydt->freeAntibiotic = (scratchSumReverse - scratchSumForward)- scratchAntibioticTarget;

	// Calculation of $\frac{dT}{dt}$
	dydt->freeTarget = scratchSumDeath1 - scratchAntibioticTarget;
	dydt->freeBoundComplex = scratchSumDeath2 + scratchAntibioticTarget;
	
	return GSL_SUCCESS;
}

/**
 * This function goes through all the parameters and checks whether any of them fall out of range.
 *
 * @param param  The ModelParamerers structure to check
 *
 * @return       0 on success, otherwise -1
 */
int sanityCheckModelParameters(ModelParameters param) {
	int goodFlag = 0;
	if (param->baselineReplication < 0) {
		fprintf(stderr, "Baseline replication was out of range. Must be R > 0.\n");
		--goodFlag;
	}
	if (param->targetMoleculeCount < 0) {
		fprintf(stderr, "Number of target molecules per-cell was out of range. Must be n > 0.\n");
		--goodFlag;
	}
	if (param->replicationThreshold < 0 || param->replicationThreshold > param->targetMoleculeCount) {
		fprintf(stderr, "Replication threshold was out of range. Must be 0 < r <= n.\n");
		--goodFlag;
	}
	if (param->maximumKillRate < 0 || param->maximumKillRate > 1.0) {
		fprintf(stderr, "Maximum killing rate was out of range. Must be 0 < K < 1.\n");
		--goodFlag;
	}
	if (param->killingThreshold < 0 || param->killingThreshold > param->targetMoleculeCount) {
		fprintf(stderr, "Killing threshold was out of range. Must be 0 < k <= n.\n");
		--goodFlag;
	}
	if (param->targetAssociationRate < 0.0) {
		fprintf(stderr, "Target association rate was out of range. Must be A > 0.\n");
		--goodFlag;
	}
	if (param->targetDissociationRate < 0.0) {
		fprintf(stderr, "Target dissociation rate was out of range. Must be D > 0.\n");
		--goodFlag;
	}
	if (param->carryingCapacity < 0.0) {
		fprintf(stderr, "Carrying capacity was out of range. Must be C > 0.\n");
		--goodFlag;
	}
	if (param->intracellularVolume < 0.0) {
		fprintf(stderr, "Intra-cellular volume was out of range. Must be V > 0.\n");
		--goodFlag;
	}
	return goodFlag;
}
