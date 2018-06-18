#ifndef FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_H
#define FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_H
//Imcompressible solver for constant density
//and without energy equation.
//The viscosity is assumed to be a constant.
//All the input and output are dimensional,
//but the inner solver is dimensionless.
#include <sstream>

#include <dynamicMesh.h>
#include <distributedObjBase.h>
#include <matrix.h>
namespace fps
{
  template <class T>
  class imcompNSDynamic:
    public distributedObjBase
  {
    private:
      dynamicMesh<T> mesh_;
      //primitive variables: u,v,w,p
      std::vector<math::matrix> u;
      std::vector<math::matrix> up[3];
      std::vector<real> dt;
      std::vector<math::matrix> invMassMatrix;
      std::vector<math::matrix> invVolJacob;
      std::vector<math::matrix> poissionMatrixCell;
      std::vector<math::matrix> poissionMatrixFace;
      std::vector<math::matrix> rhsConvective;
      std::vector<math::matrix> rhsConvectiveP[3];
      std::vector<math::matrix> rhsSource;
      int start;

      //The Four basic dimensional variables
      //used to nondimensionalize the N-S Equations.

      real density;
      real refL;
      real refV;
      real refTemp;
      real refPre;
      real refTime;
      real Reynolds;

      //Some control variables
      bool viscous;
      //Whether to use the fixed time step or not
      bool fixedStep;
      int timeAccuracy;
      const real cofGamma[3] = {
        1.0, 1.5, 11.0/6.0
      };
      const real cofAlpha[3][3] = {
        { 1.0,    0.0,     0.0},
        { 2.0,   -0.5,     0.0},
        { 3.0,   -1.5, 1.0/3.0}
      };
      const real cofBeta[3][3] = {
        {1.0, 0.0, 0.0},
        {2.0, -1.0, 0.0},
        {3.0,-3.0,1.0},
      };

      real physicalDt;
      //If not fixedStep the CFL number will be used
      //to determine the forward timestep
      real physicalCFL;
      //Current time
      real ttime;
      //Current step
      int nStep;
      //End time
      real maxTime;
      //Max physical time steps
      int maxStep;
      //Funcitons used to determine whether to output or not
      bool readyToOutPut();
      //Funcitons used to determine whether to backup or not
      bool readyToBackup();

      int minLevel;
      int maxLevel;
    private:
      struct boundaryCondition{
        real u[3];
        real du[3][3];
        real dudn[3];
        real p;
        real dp[3];
        real dpdn[3];
      };
      std::vector<boundaryCondition> boundaryConditionList_;
    public:
      real                time_refineCoarsen;
      real                time_balance;
      real                time_run;
      real                time_mapData;

      void updateReconMatrix(void);
      void evaluateFlux(real time = 0);
      //Evaluate the delta T for inner iteration
      real evaluateTimeStep(real CFL, bool globalTimeStep);
      //For implicit method, evaluate the delta U for each equation
      //LUSGS or GMRES can be used here
      void steadyIterationLUSGS();
      void refineCoarsen();
      real time(){return ttime*refTime;}
      bool end() {return ttime >= maxTime || nStep > maxStep;}
      int nCells(){return mesh_.nCells();}
      enum bcType
      {
        bc_null = 0,
        bc_user_defined1 = 1,
        bc_extrapolate = 2,
        bc_farfield = 3,
        bc_inflow = 4,
        bc_inflow_subsonic = 5,
        bc_inflow_supersonic = 6,
        bc_outflow = 7,
        bc_outflow_pressure = 22,
        bc_outflow_subsonic = 8,
        bc_outflow_supersonic = 9,
        bc_symmetry_plane = 10,
        bc_tunnel_inflow = 11,
        bc_tunnel_outflow =12,
        bc_wall = 13,
        bc_wall_inviscid = 14,
        bc_wall_viscous = 15,
        bc_wall_viscous_heatflux = 16,
        bc_wall_viscous_isothermal = 17,
        bc_user_defined2 = 19,
        bc_inner_face = 20,
      };
    private:
      std::vector<bcType>           boundaryTypeList_;
    public:
      imcompNSDynamic(MPI_Comm mpiComm = MPI_COMM_WORLD)
        :distributedObjBase(mpiComm),mesh_(mpiComm)
      {
      }

      ~imcompNSDynamic()
      {
      }
      //Things to do:
      //  1. Read reference variables (velocity and pressure) and Re number
      //  2. Read solver control variables. dt/CFL, timeforward method, etc.
      //  3. Read mesh file name/predefined mesh file name
      //     and the reference length
      //  4: Read mesh and boundary conditions
      //  5: Allocate variables' memory
      //  6: Initialize matrices
      void initSolver(const std::string&);
      void updateMatrices();

      //Temporarily set as public function
      void explicitConvectiveStep();
      void poissonStep();
      void projectionStep();
      void viscousStep();

      //------------To be implemented later------------
      //void initFlow(const std::fstream&);
      //void doBalance();
      //void doAMR();

      //void steadyTimeForwardRK(int nMax, bool globalTimeStep);
      //void steadyIteration(int nMax, bool globalTimeStep = false, int nStepsOutPut = 1, double epson = 1e-6);
      //void outPut(const std::string& filename, bool tag = 0);
      //void backUp(const std::string& filename, bool tag = 0);
      //void postProcess(const std::string& filename);
  };
}

#include "./imcompNSDynamic_initSolver.h"
#include "./imcompNSDynamic_updateMatrices.h"
#include "./imcompNSDynamic_explicitConvectiveStep.h"
#endif
