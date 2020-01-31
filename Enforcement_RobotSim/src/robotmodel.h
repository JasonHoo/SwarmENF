/* vim:set ts=4 sw=4 sts=4 et: */

/* Tools for handling the robot model
 */

#ifndef ROBOTMODEL_H
#define ROBOTMODEL_H

#include "utilities/param_utils.h"
#include "utilities/math_utils.h"
#include "utilities/datastructs.h"
#include "vizmode.h"
#include "sensors.h"
#include "algo.h"
#include "algo_gui.h"
#include "vizmode.h"
#include "utilities/debug_utils.h"
#include <math.h>
#include <string.h>

/* Step positions and velocities
 */
void Step(phase_t * OutputPhase, phase_t * GPSPhase, phase_t * GPSDelayedPhase,
        phase_t * PhaseData, unit_model_params_t * UnitParams,
        flocking_model_params_t * FlockingParams, sit_parameters_t * SitParams,
        vizmode_params_t * VizParams, int TimeStepLooped, int TimeStepReal,
        bool CountCollisions, bool * ConditionsReset, int *Collisions,
        bool * AgentsInDanger, double *WindVelocityVector,
        double *Accelerations);

/*  Initialization and killing
 */
// added by HC -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Swarm_activeNum(phase_t *Phase, flocking_model_params_t *FlockingParams );   //监控插装函数
void EnforceMonitor(int *Collisions,sit_parameters_t * SitParams,phase_t *Phase,flocking_model_params_t * FlockingParams,  int on_off);       //  TODO：监控接口和强制接口，需在Python上重新实现

int Micro_Observer1(int ob_Agents,phase_t *Phase, flocking_model_params_t * FlockingParams, const char *string);
int Macro_Observer1(int *Collisions, phase_t *Phase, flocking_model_params_t * FlockingParams, const char *string);                         //Runtime Monitors
int Macro_Observer2(flocking_model_params_t * FlockingParams, const char *string);

void Micro_Enforcer1(phase_t *Phase, flocking_model_params_t * FlockingParams);
void macro_Enforcer1();                                                                                                                                                                                                                                         // 1-step Enforcers
void macro_Enforcer2(sit_parameters_t * SitParams);

void Fault_injection1(phase_t *Phase, flocking_model_params_t *FlockingParams, unit_model_params_t * UnitParams );                   // 故障注入函数
void Fault_injection2(phase_t Phase, unit_model_params_t * UnitParams);
// added by HC -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void InitializePreferredVelocities(phase_t * Phase,
        flocking_model_params_t * FlockingParams, sit_parameters_t * SitParams,
        unit_model_params_t * UnitParams, double *WindVelocityVector);
void freePreferredVelocities(phase_t * Phase,
        flocking_model_params_t * FlockingParams, sit_parameters_t * SitParams);

#endif
