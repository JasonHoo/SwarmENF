/* vim:set ts=4 sw=4 sts=4 et: */

/* Tools for handling the robot model
 */

#include "robotmodel.h"
#include "algo_spp_evol.h"

double **PreferredVelocities;
phase_t LocalActualPhase;
phase_t LocalActualDelayedPhase;
phase_t SteppedPhase;
phase_t TempPhase;

double *ChangedInnerStateOfActualAgent;

int injecFlag = 1;
static int NumberOfNeighbours = 0;      /* Number of units observed by the actual agent */
static int ob_AgentsEnforceArr[100];   /*ID array of enforcement needed agents*/
static int previous_ob1[100] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};  /*the previous truth value of enforce_monitor observed on  FORMULA 1 */

/* For passing debug information to CalculatePreferredVelocity function */
agent_debug_info_t DebugInfo;

/* Caulculating the phase space observed by the "WhichAgent"th unit.  */
void CreatePhase(phase_t * LocalActualPhaseToCreate,
        phase_t * GPSPhase,
        phase_t * GPSDelayedPhase,
        phase_t * Phase,
        phase_t * DelayedPhase,
        const int WhichAgent,
        const double R_C,
        const double packet_loss_ratio,
        const double packet_loss_distance, const bool OrderByDistance) {

    int i, j;
    LocalActualPhaseToCreate->NumberOfAgents = Phase->NumberOfAgents;   // ???
    LocalActualPhaseToCreate->NumberOfInnerStates = Phase->NumberOfInnerStates; // ???

    /* Setting up order by distance from actual unit */

    for (i = 0; i < Phase->NumberOfAgents; i++) {

        LocalActualPhaseToCreate->RealIDs[i] = Phase->RealIDs[i];

        for (j = 0; j < 3; j++) {
            LocalActualPhaseToCreate->Coordinates[i][j] =
                    Phase->Coordinates[i][j];
            LocalActualPhaseToCreate->Velocities[i][j] =
                    Phase->Velocities[i][j];
        }

        for (j = 0; j < Phase->NumberOfInnerStates; j++) {
            LocalActualPhaseToCreate->InnerStates[i][j] =
                    Phase->InnerStates[i][j];
        }

    }

    static double ActualAgentsPosition[3];
    GetAgentsCoordinates(ActualAgentsPosition, Phase, WhichAgent);

    if (OrderByDistance) {
        /* OrderByDistance selects the nearby agents, and places their phase to the beginning of the full phase array. */
        NumberOfNeighbours =
                SelectNearbyVisibleAgents(LocalActualPhaseToCreate,
                ActualAgentsPosition, R_C,
                packet_loss_distance >
                0 ? packet_loss_ratio / packet_loss_distance /
                packet_loss_distance : 0);
    } else {
        SwapAgents(LocalActualPhaseToCreate, WhichAgent, 0);
        NumberOfNeighbours = 1;
    }

    /* Setting up delay and GPS inaccuracy for positions and velocities */

    /* Velocities */
    static double RealVelocity[3];
    NullVect(RealVelocity, 3);
    static double GPSVelocityToAdd[3];
    NullVect(GPSVelocityToAdd, 3);
    /* Positions */
    static double RealPosition[3];
    NullVect(RealPosition, 3);
    static double GPSPositionToAdd[3];
    NullVect(GPSPositionToAdd, 3);

    for (i = 1; i < NumberOfNeighbours; i++) {

        GetAgentsCoordinates(RealPosition, DelayedPhase,
                LocalActualPhaseToCreate->RealIDs[i]);
        GetAgentsVelocity(RealVelocity, DelayedPhase,
                LocalActualPhaseToCreate->RealIDs[i]);
        InsertAgentsCoordinates(LocalActualPhaseToCreate, RealPosition, i);
        InsertAgentsVelocity(LocalActualPhaseToCreate, RealVelocity, i);

    }

    GetAgentsCoordinates(RealPosition, Phase, WhichAgent);

    /* Adding term of GPS errors (XY) */
    GetAgentsCoordinates(GPSPositionToAdd, GPSPhase, WhichAgent);       // WhichAgent = RealIDs[0] = 0
    GetAgentsVelocity(GPSVelocityToAdd, GPSPhase, WhichAgent);
    GetAgentsVelocity(RealVelocity, LocalActualPhaseToCreate, 0);
    GetAgentsCoordinates(RealPosition, LocalActualPhaseToCreate, 0);
    VectSum(RealVelocity, RealVelocity, GPSVelocityToAdd);
    VectSum(RealPosition, RealPosition, GPSPositionToAdd);
    InsertAgentsCoordinates(LocalActualPhaseToCreate, RealPosition, 0);
    InsertAgentsVelocity(LocalActualPhaseToCreate, RealVelocity, 0);

    for (i = 1; i < NumberOfNeighbours; i++) {

        GetAgentsVelocity(RealVelocity, LocalActualPhaseToCreate, i);
        GetAgentsCoordinates(RealPosition, LocalActualPhaseToCreate, i);

        GetAgentsCoordinates(GPSPositionToAdd, GPSDelayedPhase,
                LocalActualPhaseToCreate->RealIDs[i]);
        GetAgentsVelocity(GPSVelocityToAdd, GPSDelayedPhase,
                LocalActualPhaseToCreate->RealIDs[i]);

        VectSum(RealVelocity, RealVelocity, GPSVelocityToAdd);
        VectSum(RealPosition, RealPosition, GPSPositionToAdd);

        InsertAgentsCoordinates(LocalActualPhaseToCreate, RealPosition, i);
        InsertAgentsVelocity(LocalActualPhaseToCreate, RealVelocity, i);

    }

    LocalActualPhaseToCreate->NumberOfAgents = NumberOfNeighbours;

}

/* Adding outer noise to final velocity vector */
/* Diffusive noise is a more-or-less effective model of the unknown properties of the control algorithm on the robots. */
void AddNoiseToVector(double *NoisedVector, double *NoiselessVector,
        double *RealVelocity, unit_model_params_t * UnitParams,
        const double DeltaT, double *WindVelocityVector) {

    int i;
    static double NoiseToAdd[3];
    NullVect(NoiseToAdd, 3);

    /* Random noise with Gaussian distribution */
    NoiseToAdd[0] = randomizeGaussDouble(0, 1);
    NoiseToAdd[1] = randomizeGaussDouble(0, 1);
    NoiseToAdd[2] = randomizeGaussDouble(0, 1);

    MultiplicateWithScalar(NoiseToAdd, NoiseToAdd,
            sqrt(2 * UnitParams->Sigma_Outer_XY.Value) * sqrt(DeltaT), 2);
    NoiseToAdd[2] *= sqrt(2 * UnitParams->Sigma_Outer_Z.Value) * sqrt(DeltaT);

    /* Noise is an additive term to the final output velocity */
    VectSum(NoiseToAdd, NoiselessVector, NoiseToAdd);
    for (i = 0; i < 2; i++) {
        NoisedVector[i] = NoiseToAdd[i];
    }
    NoisedVector[2] = NoiseToAdd[2];

}

/* Calculating wind effect to add to final output acceleration */
/* It is based on a simple Stokes-like force-law. */
void StepWind(unit_model_params_t * UnitParams, const double DeltaT,
        double *WindVelocityVector) {

    static double NoiseToAdd[2];
    NoiseToAdd[0] = randomizeGaussDouble(0.0, 1.0);
    NoiseToAdd[1] = randomizeGaussDouble(0.0, 1.0);

    WindVelocityVector[0] =
            WindVelocityVector[0] +
            NoiseToAdd[0] * sqrt(2 * UnitParams->Wind_StDev.Value * DeltaT);
    WindVelocityVector[1] =
            WindVelocityVector[1] +
            NoiseToAdd[1] * sqrt(2 * UnitParams->Wind_StDev.Value * DeltaT);

    //TODO: Potential for Wind speed...

}

/* Force law contains specific features of a real robot
*/
void RealCoptForceLaw(double *OutputVelocity, double *OutputInnerState,
        phase_t * Phase, double *RealVelocity, unit_model_params_t * UnitParams,
        flocking_model_params_t * FlockingParams, vizmode_params_t * VizParams,
        const double DeltaT, const int TimeStepReal, const int TimeStepLooped,
        const int WhichAgent, double *WindVelocityVector) {

    int i;

    /*
     * The units in the array "Phase" ara ordered by the distance from "WhichAgent"th unit.
     * Therefore, Position and velocity of "WhichAgent"th unit is stored in Phase[0] - Phase[5].
     */

    static double PreviousVelocity[3];
    GetAgentsVelocity(PreviousVelocity, Phase, 0);

    static double TempTarget[3];

    for (i = 0; i < Phase->NumberOfInnerStates; i++) {
        OutputInnerState[i] = Phase->InnerStates[0][i];
    }

    if (TimeStepLooped % ((int) (UnitParams->t_GPS.Value / DeltaT)) == 0) {

        /*for (i=0; i<FlockingParams->NumberOfInnerStates; i++) {
           OutputInnerState[i] = Phase->InnerStates[0][i];
           } */

        /* Calculating target velocity */
        NullVect(TempTarget, 3);

        FILE *log;

       //     log=fopen("flockingLogs.txt","a+");

     //   fprintf(log,"AgentID:%d ",11);
   //    fclose(log);

        CalculatePreferredVelocity(TempTarget, OutputInnerState, Phase, 0,
                FlockingParams, VizParams, UnitParams->t_del.Value,
                TimeStepReal * DeltaT, &DebugInfo);

        for (i = 0; i < 3; i++) {

            PreferredVelocities[WhichAgent][i] = TempTarget[i];

        }

    }

    for (i = 0; i < 2; i++) {

        OutputVelocity[i] =
                RealVelocity[i] +
                (DeltaT / UnitParams->Tau_PID_XY.Value) *
                (PreferredVelocities[WhichAgent][i] - PreviousVelocity[i]);

    }

    OutputVelocity[2] =
            RealVelocity[2] +
            (DeltaT / UnitParams->Tau_PID_Z.Value) *
            (PreferredVelocities[WhichAgent][2] - PreviousVelocity[2]);

}

/* Step positions and velocities and copy real IDs */
void Step(phase_t * OutputPhase, phase_t * GPSPhase, phase_t * GPSDelayedPhase,
        phase_t * PhaseData, unit_model_params_t * UnitParams,
        flocking_model_params_t * FlockingParams, sit_parameters_t * SitParams,
        vizmode_params_t * VizParams, int TimeStepLooped, int TimeStepReal,
        bool CountCollisions, bool * ConditionsReset, int *Collisions,
        bool * AgentsInDanger, double *WindVelocityVector,
        double *Accelerations) {

    int i, j, k;

    static double CheckVelocityCache[3];
    NullVect(CheckVelocityCache, 3);
    static double CheckAccelerationCache[3];
    NullVect(CheckAccelerationCache, 3);
    static double CheckDifferenceCache[3];
    NullVect(CheckDifferenceCache, 3);
    static double UnitVectDifference[3];
    NullVect(UnitVectDifference, 3);

    static double DelayStep;
    DelayStep = (UnitParams->t_del.Value / SitParams->DeltaT);

    /* Getting phase of actual TimeStepfrom PhaseData */
    LocalActualPhase = PhaseData[TimeStepLooped];
    LocalActualDelayedPhase = PhaseData[TimeStepLooped - (int) DelayStep];

    /* Counting Collisions */
    static int PreviousColl = 0;
    if (CountCollisions == true) {
        *Collisions +=
                HowManyCollisions(&LocalActualPhase, AgentsInDanger,
                CountCollisions, SitParams->Radius);
        PreviousColl = *Collisions;
    }

    /* Step coordinates (with velocity of previous TimeStepLooped) */
    static double Velocity[3];
    static double CoordinatesToStep[3];
    for (j = 0; j < SitParams->NumberOfAgents; j++) {

        GetAgentsVelocity(Velocity, &LocalActualPhase, j);
        GetAgentsCoordinates(CoordinatesToStep, &LocalActualPhase, j);

        for (i = 0; i < 3; i++) {

            CoordinatesToStep[i] += Velocity[i] * SitParams->DeltaT;

        }

        InsertAgentsCoordinates(&SteppedPhase, CoordinatesToStep, j);

    }

    /* Step GPS coordinates and velocities (in every "t_gps"th second) */
    if ((TimeStepLooped) % ((int) (UnitParams->t_GPS.Value /
                            SitParams->DeltaT)) == 0) {
        StepGPSNoises(GPSPhase, UnitParams);
        StepGPSNoises(GPSDelayedPhase, UnitParams);
    }

    /* Step Wind vector */
    //StepWind (UnitParams, SitParams->DeltaT, WindVelocityVector);

    /* "Realcopt" force law */

    static double RealCoptForceVector[3];
    NullVect(RealCoptForceVector, 3);
    static double ActualRealVelocity[3];
    NullVect(ActualRealVelocity, 3);

    for (j = 0; j < SitParams->NumberOfAgents; j++) {

        /* Constructing debug information about the actual agent */
        DebugInfo.AgentsSeqNumber = j;
        GetAgentsCoordinates(DebugInfo.AgentsRealPosition, &LocalActualPhase,
                j);
        GetAgentsVelocity(DebugInfo.AgentsRealVelocity, &LocalActualPhase, j);
        DebugInfo.RealPhase = &LocalActualPhase;

        /* Creating phase from the viewpoint of the actual agent */
        CreatePhase(&TempPhase, GPSPhase, GPSDelayedPhase, &LocalActualPhase,
                &LocalActualDelayedPhase, j, UnitParams->R_C.Value,
                UnitParams->packet_loss_ratio.Value,
                UnitParams->packet_loss_distance.Value,
                (TimeStepLooped % ((int) (UnitParams->t_GPS.Value /
                                        SitParams->DeltaT)) == 0));

        GetAgentsVelocity(ActualRealVelocity, &LocalActualPhase, j);

        /* Solving Newtonian with Euler-Naruyama method */
        NullVect(RealCoptForceVector, 3);

        RealCoptForceLaw(RealCoptForceVector, ChangedInnerStateOfActualAgent,
                &TempPhase,ActualRealVelocity, UnitParams, FlockingParams,
                VizParams, SitParams->DeltaT, TimeStepReal, TimeStepLooped, j,
                WindVelocityVector);

        NullVect(CheckVelocityCache, 3);
        VectSum(CheckVelocityCache, CheckVelocityCache, RealCoptForceVector);

        /* Updating inner states */
        InsertAgentsVelocity(&SteppedPhase, CheckVelocityCache, j);
        for (k = 0; k < PhaseData[0].NumberOfInnerStates; k++) {
            SteppedPhase.InnerStates[j][k] = ChangedInnerStateOfActualAgent[k];
        }

        //     Instrument of the Fault Inject and Enforcers---------------------------------------------------------------------------------------------------------------------------------------------------------------
        Fault_injection1(&SteppedPhase, FlockingParams, UnitParams );

        Micro_Enforcer1(&SteppedPhase, FlockingParams);
        //    Instrument of the Fault Inject and Enforcers----------------------------------------------------------------------------------------------------------------------------------------------------------------

    }
    double OnePerDeltaT = 1. / SitParams->DeltaT;
    /* The acceleration saturates at a_max. We save out the acceleration magnitude values before we add OuterNoise to the velocities and lose the possibility to derivate numerically */
    for (i = 0; i < SitParams->NumberOfAgents; i++) {

        GetAgentsVelocity(CheckAccelerationCache, &LocalActualPhase, i);
        GetAgentsVelocity(CheckVelocityCache, &SteppedPhase, i);
        VectDifference(CheckDifferenceCache, CheckVelocityCache,
                CheckAccelerationCache);
        UnitVect(UnitVectDifference, CheckDifferenceCache);
        Accelerations[i] = VectAbs(CheckDifferenceCache) * OnePerDeltaT;

        if (Accelerations[i] > UnitParams->a_max.Value) {

            for (k = 0; k < 3; k++) {

                CheckVelocityCache[k] =
                        CheckAccelerationCache[k] +
                        UnitParams->a_max.Value * SitParams->DeltaT *
                        UnitVectDifference[k];

            }

            InsertAgentsVelocity(&SteppedPhase, CheckVelocityCache, i);
            Accelerations[i] = UnitParams->a_max.Value;
        }

    }
    //    Instrument of the Fault Inject and Enforcers----------------------------------------------------------------------------------------------------------------------------------------------------------------
    Fault_injection2(SteppedPhase, UnitParams);
    //    Instrument of the Fault Inject and Enforcers----------------------------------------------------------------------------------------------------------------------------------------------------------------

    /* Outer Noise Term - Gaussian white noise */
    for (j = 0; j < SitParams->NumberOfAgents; j++) {
        if (true == Noises[j]) {
            GetAgentsVelocity(CheckAccelerationCache, &LocalActualPhase, j);
            GetAgentsVelocity(CheckVelocityCache, &SteppedPhase, j);
            AddNoiseToVector(CheckVelocityCache, CheckVelocityCache,
                    CheckAccelerationCache, UnitParams, SitParams->DeltaT,
                    WindVelocityVector);
            InsertAgentsVelocity(&SteppedPhase, CheckVelocityCache, j);
        }
    }

    /* Redistribution of agents when pressing F12 */
    if (ConditionsReset[0] == true) {
        DestroyPhase(&SteppedPhase, FlockingParams, SitParams);
        RandomizePhase(&SteppedPhase, SitParams->InitialX, SitParams->InitialY,
                SitParams->InitialZ, VizParams->CenterX, VizParams->CenterY,
                VizParams->CenterZ, 0, SitParams->NumberOfAgents,
                SitParams->Radius);
        InitializePhase(&SteppedPhase, FlockingParams, SitParams);
        ModelSpecificReset(&SteppedPhase, SitParams->InitialX,
                SitParams->InitialY, SitParams->InitialZ, VizParams,
                FlockingParams, SitParams->Radius);
        ConditionsReset[0] = false;
    } else if (ConditionsReset[1] == true) {
        if (VizParams->MapSizeXY < SitParams->InitialX
                || VizParams->MapSizeXY < SitParams->InitialY) {
            DestroyPhase(&SteppedPhase, FlockingParams, SitParams);
            RandomizePhase(&SteppedPhase, SitParams->InitialX,
                    SitParams->InitialY, SitParams->InitialZ,
                    VizParams->CenterX, VizParams->CenterY, VizParams->CenterZ,
                    0, SitParams->NumberOfAgents, SitParams->Radius);
            InitializePhase(&SteppedPhase, FlockingParams, SitParams);
            ModelSpecificReset(&SteppedPhase,
                    SitParams->InitialX, SitParams->InitialY,
                    SitParams->InitialZ, VizParams, FlockingParams,
                    SitParams->Radius);
        } else {
            DestroyPhase(&SteppedPhase, FlockingParams, SitParams);
            RandomizePhase(&SteppedPhase, VizParams->MapSizeXY,
                    VizParams->MapSizeXY, VizParams->MapSizeZ,
                    VizParams->CenterX, VizParams->CenterY, VizParams->CenterZ,
                    0, SitParams->NumberOfAgents, SitParams->Radius);
            InitializePhase(&SteppedPhase, FlockingParams, SitParams);
            ModelSpecificReset(&SteppedPhase, VizParams->MapSizeXY,
                    VizParams->MapSizeXY, VizParams->MapSizeZ, VizParams,
                    FlockingParams, SitParams->Radius);
        }
        ConditionsReset[1] = false;
    }

    /* Insert Phase into PhaseData... */
    for (j = 0; j < SitParams->NumberOfAgents; j++) {
        for (i = 0; i < 3; i++) {
            OutputPhase->Coordinates[j][i] = SteppedPhase.Coordinates[j][i];
            OutputPhase->Velocities[j][i] = SteppedPhase.Velocities[j][i];
        }
        for (i = 0; i < SteppedPhase.NumberOfInnerStates; i++) {
            OutputPhase->InnerStates[j][i] = SteppedPhase.InnerStates[j][i];
        }
    }

//     Instrument of the Swarm_innerDistance calculation ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    Swarm_activeNum( OutputPhase, FlockingParams);
       //printf("Now the  active robot number of this swarm is  %d  \n", FlockingParams->NumberofActiveAgents);

//     Instrument of the Enforcement Observer---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    EnforceMonitor(Collisions, SitParams,&SteppedPhase, FlockingParams, 1);
}

//added by HC
void Swarm_activeNum(phase_t *Phase, flocking_model_params_t *FlockingParams ){
	    int i;
	    int flag =0;
//        double CheckVelocityCache[3];
//        NullVect(CheckVelocityCache, 3);
//	    InsertAgentsVelocity(&SteppedPhase, CheckVelocityCache, 2);
		for( i = 0; i<100;  i++){
			    double *tmpVelocities = Phase->Velocities[i];
			    double tmpVeloAbs = VectAbs(tmpVelocities);
			    double veloFilter= 20.0;
		   	    if( (tmpVeloAbs < veloFilter) && FlockingParams->WatchTimer > 10000){
		   	                 	FlockingParams->NumberofActiveAgents =   FlockingParams->NumberofActiveAgents  -1;
		   	                    FlockingParams->WatchTimer = 0;
		    	                printf("Now the Agent Number of this swarm is  %d .  The Agent ID %d  is  broken! \n", FlockingParams->NumberofActiveAgents, i);
		    	                flag = 1;
		        }
		   	    else if((tmpVeloAbs < veloFilter) && FlockingParams->WatchTimer <= 10000){
		   	    				FlockingParams->WatchTimer = FlockingParams->WatchTimer +1;
		   			    	   // printf("Now the WatchTimer is  %d \n", FlockingParams->WatchTimer);
		   			    	    flag = 1;
		   	   }
		    }
		if (flag == 0){
			FlockingParams->WatchTimer = 0;
		}
}

void EnforceMonitor(int *Collisions,sit_parameters_t * SitParams,phase_t *Phase,flocking_model_params_t * FlockingParams,  int on_off){
		if(on_off == 1){
			int i;
			int ob1[100];
			for(i=0; i<100; i++){
				ob1[i] = Micro_Observer1(i, Phase, FlockingParams,"(V_Flock<500 and activeAgentNum >=3) --> ( all Phase.neighborNum >= 2)");

				if(ob1[i] ==0 && previous_ob1[i] == 1){
						ob_AgentsEnforceArr[i] = 1;
						printf("ob of %d  =  %d   ,  ", i,  ob1[i]);
					    printf("ob_previous of %d  =  %d \n", i,  previous_ob1[i]);
				}
//				else if (ob1[i] ==1 && previous_ob1[i] == 1){
//						ob_AgentsEnforceArr[i] = 0;
//				}
				const int tmp = ob1[i];
				previous_ob1[i] = tmp;  //1-step enforcement and only deal with the first  deviation
			}



			int ob2 = Macro_Observer1(Collisions, Phase, FlockingParams, "(Collisions >=3  or V_Flock>=500 ) --> ( FlockingParams->innerDistance >= 200)");
			if(ob2 == 0 && V_Flock>=500 ){
						macro_Enforcer1();
						printf("here is the Macro_enforcer1 !  Enforcing the V_Flock!   \n");
			}

			int ob3 = Macro_Observer2(FlockingParams,"( FlockingParams->SwarmLocation !=  Dest ) --> (FlockingParams->NumberofActiveAgents >= 75)");
			if(ob3 == 0 ){
						printf("here is the Macro_enforcer2 !  Enforcing the  regrouping!   \n");
						macro_Enforcer2(SitParams);
			}
	   	}
		else {

		}

}

void Micro_Enforcer1(phase_t *Phase, flocking_model_params_t * FlockingParams ){

//			double *AgentsVelocity;
			double *TmpCenterCoord;
			double *TmpAgentCoord;
			int i;
			for(i=0; i<100; i++) {
				if (ob_AgentsEnforceArr[i] == 1) {
					double *TmpVelocity = Phase->Velocities[i];
//					AgentsVelocity= Phase.Velocities[i];
//							TmpVelocity[0] = -10*AgentsVelocity[0];
//							TmpVelocity[1] = -10*AgentsVelocity[1];
//							TmpVelocity[2] = 0;

					TmpAgentCoord = Phase->Coordinates[i];
					TmpCenterCoord = FlockingParams->SwarmLocation;
					TmpVelocity[0] = TmpCenterCoord[0] - TmpAgentCoord[0];
					TmpVelocity[1] = TmpCenterCoord[1] - TmpAgentCoord[1];
					TmpVelocity[2] = 0;
					UnitVect(TmpVelocity, TmpVelocity);
					MultiplicateWithScalar(TmpVelocity, TmpVelocity, V_Flock*20, 2);

					Phase->Velocities[i] = TmpVelocity;

					ob_AgentsEnforceArr[i] = 0;

					//printf("here is the enforcer1 !  Enforcing the Robot %d ------------------------------------------------- \n", i);
				}
			}
}

void macro_Enforcer1(){

		V_Flock = V_Flock/2;
}

void macro_Enforcer2(sit_parameters_t * SitParams ){

		//SitParams->NumberOfAgents = SitParams->NumberOfAgents -88;

}

void Fault_injection1(phase_t *Phase, flocking_model_params_t *FlockingParams, unit_model_params_t * UnitParams ){
	    double *tmpVelocity = Phase->Velocities[1];

	    int i;
	    for(i=0; i<100; i++) {
	    				if (ob_AgentsEnforceArr[i] == 1) {
	    							injecFlag = 0;
	    							break;
	    				}
	    }
		if(injecFlag ==1 &&
				UnitParams->Robot1_Coord_X.Value == 20.0){
				tmpVelocity[0] = -900;
				tmpVelocity[1] = -900;
				tmpVelocity[2] = 0;
		}
		else{

		}
		Phase->Velocities[1] = tmpVelocity;
		//printf("here is the Robot1 coord !  x = %f   ,   y = %f ------------------------------------------------- \n", Phase->Coordinates[1][0] ,Phase->Coordinates[1][1] );
}

void Fault_injection2(phase_t Phase, unit_model_params_t * UnitParams ){
		double tmpVelocity[3];
		tmpVelocity[0] = 0;
		tmpVelocity[1] = 0;
		tmpVelocity[2] = 0;
		int ie;
	    if(UnitParams->Robot1_Coord_X.Value == 30.0){
	    	for(ie=0; ie<30; ie++) {

						Phase.Velocities[ie] = tmpVelocity;
	    	}
	    }
}

int Micro_Observer1(int ob_Agents,phase_t *Phase, flocking_model_params_t * FlockingParams, const char *string){

		double *TmpAgentsCoordinates;
	    double *NeighboursCoordinates;
	    double DifferenceVector[3];
	    NullVect(DifferenceVector, 3);
		double DistanceFromNeighbour = 0.0;
		//int deviation_flag = 0;
		char stringFormula[] =  "(V_Flock<500 and activeAgentNum >=3) --> ( all Phase.neighborNum >= 2)";
		if(strcmp(string,stringFormula)==0 ){
				if(V_Flock >= 500 ){
                        return 1;
				}
				else if (FlockingParams->NumberofActiveAgents < 3){
					    return 1;
				}
				else{
                        // The calculation for every robot agent -- here, n denotes the neighborNum
					    TmpAgentsCoordinates = Phase->Coordinates[ob_Agents];
                        int it;
                        int n = 0;
                        for (it = 0; it < 100; it++) {
                        	if (it == ob_Agents)
                        				continue;
                        	NeighboursCoordinates = Phase->Coordinates[it];
                        	VectDifference(DifferenceVector, TmpAgentsCoordinates, NeighboursCoordinates);
                        	if (2 == Dim) {
                        	        	DifferenceVector[2] = 0.0;
                        	}
                        	DistanceFromNeighbour = VectAbs(DifferenceVector);

                        	        /* Compare with  the distance of inner-communication  */
                        	if (DistanceFromNeighbour >=10 * 2*0.5*R_0) {
                        	        	continue;
                        	}
                            n += 1;
                        }

                        if (n<2){
                        			//printf("The deviation agent is  %d , neighborNum = %d \n", ob_Agents, n);
                        			return 0;
                         }
                        else{
                        			return 1;
                        }
			    }
		}
		else{
				return -1;
		}
};

int Macro_Observer1(int *Collisions, phase_t *Phase, flocking_model_params_t * FlockingParams, const char *string){

	//int deviation_flag = 0;
	char stringFormula[] =  "(Collisions >=3  or V_Flock>=500 ) --> ( FlockingParams->innerDistance >= 200)";
	if(strcmp(string,stringFormula)==0 ){
			if(V_Flock < 500 && *Collisions < 3 ){
                    return 1;
			}
			else{
                    // The calculation for innerDistance --
                    if (FlockingParams->innerDistance < 200){
                    			//printf("The deviation agent is  %d , neighborNum = %d \n", ob_Agents, n);
                    			return 0;
                     }
                    else{
                    			return 1;
                    }
		    }
	}
	else{
			return -1;
	}
};

int Macro_Observer2(flocking_model_params_t * FlockingParams, const char *string){

	double *TmpSwarmLoc = FlockingParams->SwarmLocation;
	//int deviation_flag = 0;
	char stringFormula[] =  "( FlockingParams->SwarmLocation !=  Dest ) --> (FlockingParams->NumberofActiveAgents >= 75)";
	if(strcmp(string,stringFormula)==0 ){

			if((40000.0 < TmpSwarmLoc[0] &&  TmpSwarmLoc[0] < 60000.0) && (40000.0 < TmpSwarmLoc[1] && TmpSwarmLoc[1]  < 60000.0) ){
                    return 1;
			}
			else{
                    // The calculation for innerDistance --
                    if (FlockingParams->NumberofActiveAgents < 75){
                    			//printf("The deviation agent is  %d , neighborNum = %d \n", ob_Agents, n);
                    			return 0;
                     }
                    else{
                    			return 1;
                    }
		    }
	}
	else{
			return -1;
	}
};

/* Some global variables are allocated and initialized in this function */
void InitializePreferredVelocities(phase_t * Phase,
        flocking_model_params_t * FlockingParams, sit_parameters_t * SitParams,
        unit_model_params_t * UnitParams, double *WindVelocityVector) {

    /* Preferred Velocities */
    PreferredVelocities = doubleMatrix(SitParams->NumberOfAgents, 3);
    Noises = BooleanData(SitParams->NumberOfAgents);

    int i, j;
    for (i = 0; i < SitParams->NumberOfAgents; i++) {
        for (j = 0; j < 3; j++) {
            PreferredVelocities[i][j] = 0.0;
        }
        Noises[i] = true;
    }

    /* Some helper dynamic arrays should be allocated here */
    AllocatePhase(&SteppedPhase, SitParams->NumberOfAgents,
            Phase->NumberOfInnerStates);
    AllocatePhase(&TempPhase, SitParams->NumberOfAgents,
            Phase->NumberOfInnerStates);
    ChangedInnerStateOfActualAgent =
            (double *) calloc(Phase->NumberOfInnerStates, sizeof(double));

    /* Setting up wind velocity vector */
    WindVelocityVector[0] = cos(UnitParams->Wind_Angle.Value);
    WindVelocityVector[1] = sin(UnitParams->Wind_Angle.Value);

}

/* Free every variable allocated in InitializePreferredVelocities function */
void freePreferredVelocities(phase_t * Phase,
        flocking_model_params_t * FlockingParams,
        sit_parameters_t * SitParams) {

    freeMatrix(PreferredVelocities, SitParams->NumberOfAgents, 3);
    free(Noises);

    /* Freeing memory owned by helper arrays */
    freePhase(&SteppedPhase);
    freePhase(&TempPhase);

    if (Phase->NumberOfInnerStates != 0) {
        free(ChangedInnerStateOfActualAgent);
    }

}
