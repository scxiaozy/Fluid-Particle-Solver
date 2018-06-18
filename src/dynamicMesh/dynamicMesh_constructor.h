#ifndef FPS_DYNAMICMESH_DYNAMICMESH_CONSTRUCTOR_H
#define FPS_DYNAMICMESH_DYNAMICMESH_CONSTRUCTOR_H
#include <string>
#include <mpi.h>
#include <p8est_connectivity.h>
#include <dynamicMesh.h>
namespace fps
{
  template <class T>
    dynamicMesh<T>::dynamicMesh(MPI_Comm mpiComm)
      :mesh(mpiComm),
      initialized(false),
      connPtr_(nullptr),
      withPeriodicBoundary_(false),
      p8estPtr_(nullptr),
      p8estGhostPtr_(nullptr),
      ghostData_(nullptr),
      referenceLength(1.0)
  {
    int k = T::refFrame::order;
    int o1 = k+1;
    int o2 = 3*k/2+1;
    if(o1==o2) oversetFaceIntPt = false;
    else oversetFaceIntPt = true;
  }

  template<class T>
    dynamicMesh<T>::dynamicMesh
    (
     typename dynamicMesh<T>::predefinedMesh t,
     int initLevel,
     real refL,
     MPI_Comm mpicomm
     )
    : mesh(mpicomm)
    {
      int k = T::refFrame::order;
      int o1 = k+1;
      int o2 = 3*k/2+1;
      if(o1==o2) oversetFaceIntPt = false;
      else oversetFaceIntPt = true;
      initialized = false;
      referenceLength = refL;
      initMesh(t, initLevel, referenceLength);
    }

  template <class T>
    dynamicMesh<T>::dynamicMesh
    (
     const std::string& fname,
     int initLevel,
     real refL,
     MPI_Comm mpiComm
    )
    : mesh(mpiComm)
    {
      int k = T::refFrame::order;
      int o1 = k+1;
      int o2 = 3*k/2+1;
      if(o1==o2) oversetFaceIntPt = false;
      else oversetFaceIntPt = true;
      initialized = false;
      referenceLength = refL;
      initMesh(fname, initLevel);
    }

  template<class T>
    void dynamicMesh<T>::resetMesh(){
      this->mesh::resetMesh();
      //Free the memory of p8est and set the pointer to nullptr
      P4EST_FREE(ghostData_);
      ghostData_ = nullptr;
      p8est_ghost_destroy(p8estGhostPtr_);
      p8estGhostPtr_ = nullptr;
      p8est_connectivity_destroy(connPtr_);
      connPtr_=nullptr;
      p8est_destroy(p8estPtr_);
      p8estPtr_ = nullptr;
      this->initialized = false;
      referenceLength = 1.0;
    }

  template<class T>
    dynamicMesh<T>::~dynamicMesh()
    {
      //Free the pointer by p8est
      P4EST_FREE(ghostData_);
      if(p8estGhostPtr_)
        p8est_ghost_destroy(p8estGhostPtr_);
      if(connPtr_)
        p8est_connectivity_destroy(connPtr_);
      if(p8estPtr_)
        p8est_destroy(p8estPtr_);
    };

  template<class T>
    void p8estInitGeom(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q)
    {
      typename dynamicMesh<T>::p8estData *data
        = (typename dynamicMesh<T>::p8estData*)q->p.user_data;
      for (int ii = 0; ii < P8EST_FACES; ++ii){
        data->patchIdx[ii] = -1;
      }
      data->p.fineCoarsenTag = -2;//New cell
    }

  template<class T>
  int dynamicMesh<T>::initMeshPeriodicCube
  (
   int initLevel,
   real refL
   )
  {
    if(initialized){
      std::cout << "Mesh has already been initialized! " << std::endl;
      exit(-1);
    }
    connPtr_       = p8est_connectivity_new_periodic();
    referenceLength = refL;
    const real3 len(1.0,1.0,1.0);
    const double        vertices[8 * 3] = {
      0, 0, 0,
      len[0]/refL, 0, 0,
      0, len[1]/refL, 0,
      len[0]/refL, len[1]/refL, 0,
      0, 0, len[2]/refL,
      len[0]/refL, 0, len[2]/refL,
      0, len[1]/refL, len[2]/refL,
     len[0]/refL, len[1]/refL, len[2]/refL,
    };

    for(int ii = 0; ii < 24; ++ii) connPtr_->vertices[ii] = vertices[ii];
    this->highOrderNodes.reserve(8);
    for(int kk = 0; kk < 8; ++kk){
      this->highOrderNodes.push_back(
          real3(
            vertices[3*kk+0],
            vertices[3*kk+1],
            vertices[3*kk+2]
            )
          );
    }
    this->cellToHiOrderNodes.resize(1);
    this->cellToHiOrderNodes[0].resize(8);
    for (int ii = 0; ii < 8; ++ii){
      this->cellToHiOrderNodes[0][ii] = ii;
    }
    //Boundary Property for dummy -1 face
    boundaryPropertyList_[-1] = boundaryProperty();

    std::vector< std::vector<std::pair<int,int> > > boundary_list;
    std::map<int,int> physicalId_;
    withPeriodicBoundary_ = true;
    std::vector< std::vector<std::pair<int,int> > > periodic_list;
    const int numPatches = 0;
    const int numPeriodic = 6;
    periodic_list.resize(numPeriodic);
    for(int ii = 0; ii < numPeriodic; ++ii){
      int pid = ii;
      physicalId_[ii] = pid;

      real3 trans;
      int iface = ii/2;
      int direc = ii%2==0?+1:-1;
      trans[iface] += direc*len[iface];

      boundaryPropertyList_[pid] = boundaryProperty(1, _periodicValue(1, real3(0,0,0), real3(0,0,0), trans), nullptr);
    }
    p8estPtr_      = p8est_new_ext(this->mpiComm(), connPtr_, 0, initLevel, 1, sizeof(p8estData), p8estInitGeom<T>, this);

    sc_array_t         *trees = p8estPtr_->trees;
    p8est_tree_t       *tree;
    p8est_quadrant_t   *quad;
    p8estData          *ud;
    size_t              si, n_quads;
    sc_array_t         *quadrants;

    for (int which_face = 0; which_face < 6; ++which_face){
      int t = 0;
      tree = p8est_tree_array_index (trees, t);
      quadrants = &(tree->quadrants);
      n_quads = quadrants->elem_count;
      for (si = 0; si < n_quads; si++) {
        quad = p8est_quadrant_array_index (quadrants, si);
        //Set the physical Id and offset of this quadrant
        ud = (typename dynamicMesh<T>::p8estData*)quad->p.user_data;
        if(quad->x==0&&which_face==0)
          ud->patchIdx[0] = physicalId_[which_face];
        if(quad->y==0&&which_face==2)
          ud->patchIdx[2] = physicalId_[which_face];
        if(quad->z==0&&which_face==4)
          ud->patchIdx[4] = physicalId_[which_face];
        if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
          ud->patchIdx[1] = physicalId_[which_face];
        if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
          ud->patchIdx[3] = physicalId_[which_face];
        if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
          ud->patchIdx[5] = physicalId_[which_face];
      }
    }
    MPI_Barrier(this->mpiComm());
    p8estGhostPtr_ = p8est_ghost_new(p8estPtr_, P8EST_CONNECT_FACE);
    if(p8estGhostPtr_->ghosts.elem_count){
      ghostData_     = P4EST_ALLOC(p8estData, p8estGhostPtr_->ghosts.elem_count);
      p8est_ghost_exchange_data(p8estPtr_, p8estGhostPtr_, ghostData_);
    }
    initialized = true;
    oldNCells_ = 0;
    oldNGhostCells_ = 0;
    updateMesh();
    updateMapVector();
    return 0;
  }
  template<class T>
  int dynamicMesh<T>::initMesh
  (
   predefinedMesh t, int initLevel, real refL
   )
  {
    switch(t){
      case periodicCube:
        initMeshPeriodicCube(initLevel, refL);
        break;
      default:
        std::cout << "Unimplemented mesh type! Abort!" << std::endl;
        break;
    }
    return 0;
  }

  template<class T>
    int dynamicMesh<T>::initMesh
    (
     const std::string& fname,
     int initLevel,
     real refL
     )
    {
      if(initialized){
        std::cout << "Mesh has already been initialized! " << std::endl;
        MPI_Barrier(this->mpiComm());
        MPI_Abort(this->mpiComm(), -1);
      }
      referenceLength = refL;
      const std::string p8estFileName = fname+".p8estConn";
      size_t connSize;
      if(!this->mpiRank())std::cout << "Loading p8est connectivity file" << std::endl;
      connPtr_         = p8est_connectivity_load(p8estFileName.c_str(),&connSize);
      if(!this->mpiRank())std::cout << "Loading over!" << std::endl;

      //scale
      for(int ii = 0; ii < connPtr_->num_vertices*3; ++ii) connPtr_->vertices[ii] /= refL;

      std::string hiONodesFileName = fname+".hiONodes";
      std::ifstream hiONodesStream;
      hiONodesStream.open(hiONodesFileName.c_str());
      //Element type
      int type;
      hiONodesStream >> type;
      if(T::refFrame::type != type){
        std::cout << "The element type of the mesh file does not coordinate with the declared dynamicMesh class"
          << std::endl;
      }

      int numHiONodesPerCell = T::refFrame::nVertices;

      long int num_vertices;
      int num_trees;
      hiONodesStream >> num_vertices;
      hiONodesStream >> num_trees;
      highOrderNodes.resize(num_vertices);
      for (int ii = 0; ii < num_vertices; ++ii){
        real3 temp;
        hiONodesStream >> temp[0] >> temp[1] >> temp[2];
        //scale
        highOrderNodes[ii] = temp/refL;
      }
      cellToHiOrderNodes.resize(num_trees);
      for (int ii = 0; ii < num_trees; ++ii){
        cellToHiOrderNodes[ii].resize(numHiONodesPerCell);
        long int idxTemp;
        for (int jj = 0; jj < numHiONodesPerCell; ++jj){
          hiONodesStream >> idxTemp;
          cellToHiOrderNodes[ii][jj] = idxTemp;
        }
      }
      hiONodesStream.close();

      std::string   bcFileName = fname+".bc";
      std::ifstream bcFileStream;
      bcFileStream.open(bcFileName.c_str());
      int numPatches;
      int numPeriodic;
      bcFileStream >> numPatches;
      for(int ii = 0; ii < numPatches; ++ii){
        int physicalId;
        std::string name("");
        bcFileStream >> physicalId;
        char ch;
        bcFileStream >> std::noskipws;
        do{
          bcFileStream >> ch;
        }while(ch!='"');
        while(1){
          bcFileStream >> ch;
          if(ch!='"') name += ch;
          else break;
        }
        bcFileStream >> std::skipws;
      }

      //Boundary Property for dummy -1 face
      boundaryPropertyList_[-1] = boundaryProperty();
      boundaryPropertyList_[-1].isPeriodic = false;
      boundaryPropertyList_[-1].periodicValue.type = 1;
      boundaryPropertyList_[-1].periodicValue.rotateCenter = real3(0,0,0);
      boundaryPropertyList_[-1].periodicValue.rotateAngle  = real3(0,0,0);
      boundaryPropertyList_[-1].periodicValue.translation  = real3(0,0,0);

      std::vector< std::vector<std::pair<int,int> > > boundary_list;
      std::map<int,int> physicalId_;
      int maxPid = -1;
      if(numPatches){
        boundary_list.resize(numPatches);
        for (int ii = 0; ii < numPatches; ++ii){
          int numbcs;
          int pid;
          int temp;
          bcFileStream >> temp >>  pid  >> numbcs;
          maxPid = maxPid < pid ? pid : maxPid;
          physicalId_[ii] = pid;
          boundaryPropertyList_[pid] = boundaryProperty();
          boundaryPropertyList_[pid].isPeriodic = false;
          boundary_list[ii].resize(numbcs);
          for (int jj = 0; jj < numbcs; ++jj){
            bcFileStream >> boundary_list[ii][jj].first >> boundary_list[ii][jj].second;
          }
        }
      }
      bcFileStream >> numPeriodic;
      std::vector< std::vector<std::pair<int,int> > > periodic_list;
      if(numPeriodic)
      {
        withPeriodicBoundary_ = true;
        periodic_list.resize(numPeriodic);
        for(int ii = 0; ii < numPeriodic; ++ii){
          int temp;
          bcFileStream >> temp;
          int pid = maxPid + ii + 1;
          physicalId_[ii+numPatches] = pid;
          boundaryPropertyList_[pid] = boundaryProperty();
          boundaryPropertyList_[pid].isPeriodic = true;
          real3 rc;
          real3 ra;
          real3 trans;
          int numPeriodicface;
          bcFileStream >> boundaryPropertyList_[pid].periodicValue.type
            >> rc[0] >> rc[1] >> rc[2] >> ra[0] >> ra[1] >> ra[2] >> trans[0] >> trans[1] >> trans[2] >> numPeriodicface;
          this->boundaryPropertyList_[pid].periodicValue.rotateCenter = rc;
          this->boundaryPropertyList_[pid].periodicValue.rotateAngle = ra;
          this->boundaryPropertyList_[pid].periodicValue.translation = trans;
          periodic_list[ii].resize(numPeriodicface);
          for (int jj = 0; jj < numPeriodicface; ++jj){
            bcFileStream >> periodic_list[ii][jj].first >> periodic_list[ii][jj].second;
          }
        }
      }
      bcFileStream.close();
      p8estPtr_      = p8est_new_ext(this->mpiComm(), connPtr_, 0, initLevel, 1, sizeof(p8estData), p8estInitGeom<T>, this);

      p4est_topidx_t      t;
      p4est_topidx_t      first_local_tree = p8estPtr_->first_local_tree;
      p4est_topidx_t      last_local_tree = p8estPtr_->last_local_tree;
      sc_array_t         *trees = p8estPtr_->trees;
      p8est_tree_t       *tree;
      p8est_quadrant_t   *quad;
      p8estData          *ud;
      size_t              si, n_quads;
      sc_array_t         *quadrants;

      for (int ii = 0; ii < numPatches; ++ii){
        for (int jj = 0; jj < boundary_list[ii].size(); ++jj){
          int t = boundary_list[ii][jj].first;
          int which_face = boundary_list[ii][jj].second;
          tree = p8est_tree_array_index (trees, t);
          quadrants = &(tree->quadrants);
          n_quads = quadrants->elem_count;
          for (si = 0; si < n_quads; si++) {
            quad = p8est_quadrant_array_index (quadrants, si);
            //Set the physical Id and offset of this quadrant
            ud = (typename dynamicMesh<T>::p8estData*)quad->p.user_data;
            if(quad->x==0&&which_face==0)
              ud->patchIdx[0] = physicalId_[ii];
            if(quad->y==0&&which_face==2)
              ud->patchIdx[2] = physicalId_[ii];
            if(quad->z==0&&which_face==4)
              ud->patchIdx[4] = physicalId_[ii];
            if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
              ud->patchIdx[1] = physicalId_[ii];
            if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
              ud->patchIdx[3] = physicalId_[ii];
            if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
              ud->patchIdx[5] = physicalId_[ii];
          }
        }
      }

      for (int ii = 0; ii < numPeriodic; ++ii){
        for (int jj = 0; jj < periodic_list[ii].size(); ++jj){
          int t = periodic_list[ii][jj].first;
          int which_face = periodic_list[ii][jj].second;
          tree = p8est_tree_array_index (trees, t);
          quadrants = &(tree->quadrants);
          n_quads = quadrants->elem_count;
          for (si = 0; si < n_quads; si++) {
            quad = p8est_quadrant_array_index (quadrants, si);
            //Set the physical Id and offset of this quadrant
            ud = (typename dynamicMesh<T>::p8estData*)quad->p.user_data;
            if(quad->x==0&&which_face==0)
              ud->patchIdx[0] = physicalId_[ii+numPatches];
            if(quad->y==0&&which_face==2)
              ud->patchIdx[2] = physicalId_[ii+numPatches];
            if(quad->z==0&&which_face==4)
              ud->patchIdx[4] = physicalId_[ii+numPatches];
            if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
              ud->patchIdx[1] = physicalId_[ii+numPatches];
            if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
              ud->patchIdx[3] = physicalId_[ii+numPatches];
            if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
              ud->patchIdx[5] = physicalId_[ii+numPatches];
          }
        }
      }
      MPI_Barrier(this->mpiComm());
      p8estGhostPtr_ = p8est_ghost_new(p8estPtr_, P8EST_CONNECT_FACE);
      if(p8estGhostPtr_->ghosts.elem_count){
        ghostData_     = P4EST_ALLOC(p8estData, p8estGhostPtr_->ghosts.elem_count);
        p8est_ghost_exchange_data(p8estPtr_, p8estGhostPtr_, ghostData_);
      }
      initialized = true;
      oldNCells_ = 0;
      oldNGhostCells_ = 0;
      updateMesh();
      updateMapVector();
      return 0;
    }
}
#endif
