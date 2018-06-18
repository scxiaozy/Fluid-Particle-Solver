#ifndef FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_INITSOLVER_H
#define FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_INITSOLVER_H
#include <libconfig.h++>
#include <cstdlib>

#include <imcompNSDynamic.h>
#include <dataType.h>
#include <matrix.h>
namespace fps
{
  template<class T>
  void imcompNSDynamic<T>::initSolver(const std::string& cfgname){
    using namespace libconfig;
    using namespace std;
    Config cfg;
    try
    {
      cfg.readFile(cfgname.c_str());
    }
    catch(const FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file." << std::endl;
      exit(-1);
    }
    catch(const ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << std::endl;
      exit(-1);
    }
    std::string predefinedMesh;
    std::string meshFileName;
    bool e1 = false, e2 = false;
    e1 = cfg.lookupValue("predefinedMesh", predefinedMesh);
    e2 = cfg.lookupValue("meshFileName",   meshFileName);
    if(e1&&e2){
      std::cerr << "Duplicate specification of predefinedMesh and meshFileName"
        << "\nOnly one of these two can be specified\n";
      exit(-1);
    }else if(!(e1||e2)){
      std::cerr << "Mesh name must be specified by defining predefinedMesh or meshFileName" << std::endl;
    }
    int initLevel;
    if(!cfg.lookupValue("initLevel", initLevel)) initLevel = 1;
    if(!cfg.lookupValue("minLevel", minLevel)) minLevel = initLevel;
    if(!cfg.lookupValue("maxLevel", maxLevel)) maxLevel = minLevel;
    if(!cfg.lookupValue("referenceLength", refL)) refL = 1.0;
    if(!cfg.lookupValue("referenceVelocity", refV)) refV = 1.0;
    refTime = refL/refV;

    if(!cfg.lookupValue("density", density)) density = 1.0;
    if(!cfg.lookupValue("Reynolds", Reynolds)){
      std::cerr << Reynolds << " Reynolds must be specified!" << std::endl;
      exit(-1);
    }
    if(!cfg.lookupValue("viscous", viscous)) viscous = true;

    bool t1 = cfg.lookupValue("physicalCFL", physicalCFL);
    bool t2 = cfg.lookupValue("physicalDt", physicalDt);

    if(!(t1||t2)){
      std::cerr << "CFL number (keyword: physicalCFL) or physical dt (keyword: physialDt) for time forward is required!"
        << std::endl;
      exit(-1);
    }
    if(t2&&t1){
      std::cerr << "CFL number (keyword: physicalCFL) and physical dt (keyword: physialDt) cannot be defined together!"
        << std::endl;
      exit(-1);
    }
    if(t1){
      fixedStep = false;
    }else{
      fixedStep = true;
    }

    if(!cfg.lookupValue("currentTime", ttime)) ttime = 0.0;
    else ttime/=refTime;

    if(!cfg.lookupValue("currentStep", nStep)) nStep = 0;

    bool e3 = cfg.lookupValue("maximumTime", maxTime);
    bool e4 = cfg.lookupValue("maximumStep", maxStep);
    if(!cfg.lookupValue("timeAccuracy", timeAccuracy)) timeAccuracy = 1;
    if(!(e3||e4)){
      std::cerr << "Either maximum time ( keyworld: maximumTime ) "
        << "or maximum step (keyword: maximumStep) should be defined for the simulation."  << std::endl;
      exit(-1);
    }
    if(!e3){
      maxTime = REAL_MAX;
    }else{
      maxTime/=refTime;
    }
    if(!e4) maxStep = INT_MAX;

    if(e1){
      std::cout << "predefinedMesh " << predefinedMesh << std::endl;
      if(predefinedMesh == "periodicCube"){
        mesh_.initMesh(dynamicMesh<T>::predefinedMesh::periodicCube, initLevel, refL);
      }else{
        std::cout << "Undefined predefinedMesh type" << std::endl;
        exit(-1);
      }
    }
    if(e2){
      mesh_.initMesh(meshFileName, initLevel, refL);
    }

    //Allocate the memory for variables
    const unsigned int & ndofs = T::refFrame::nSolutionPoints;
    llabel nCells = mesh_.nCells();
    llabel nGhostCells = mesh_.nGhostCells();
    llabel nTCells = nCells + nGhostCells;

    u.resize(nTCells);
    rhsConvective.resize(nTCells);
    rhsSource.resize(nTCells);

    for (int ss = 0; ss < timeAccuracy; ++ss){
      up[ss].resize(nTCells);
      rhsConvectiveP[ss].resize(nTCells);
    }

    for (llabel nc = 0; nc < nTCells; ++nc){
      u[nc].resize(T::refFrame::nSolutionPoints, 4);
      rhsConvective[nc].resize(T::refFrame::nSolutionPoints, 3);
      rhsSource[nc].resize(T::refFrame::nSolutionPoints, 3);

      for (int ss = 0; ss < timeAccuracy; ++ss){
        up[ss][nc].resize(T::refFrame::nSolutionPoints, 4);
        rhsConvectiveP[ss][nc].resize(T::refFrame::nSolutionPoints, 3);
      }
    }
    //Set boundary conditions
    //To be determined later
    //Init fluid field
    start = 1;
    for (int ii = 0; ii < mesh_.nCells(); ++ii){
        ///Init Section
        ///To be determined later
    }
    //Init the ghost data
    ghostDataSynchronizeDeepCopy<math::matrix>(mesh_, 1, &u);

    //All the DG matrices will be handled here, including allocation and initialization
    updateMatrices();
  }
}
#endif
