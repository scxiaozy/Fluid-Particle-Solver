#ifndef FPS_DYNAMICMESH_DYNAMICMESH_H
#define FPS_DYNAMICMESH_DYNAMICMESH_H
#include <sstream>
#include <fstream>

#include <sc.h>
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#include <p8est_geometry.h>

#include <dataType.h>
#include <real3.h>
#include <mesh.h>
#include <faceElement.h>

namespace fps
{
  //Forwad declaration
  //template<unsigned int O, unsigned int MO>
  //class cfvDynamicMesh;
  template<class T>
  void updateMeshIterVol (p8est_iter_volume_info_t *info, void *user_data);
  template<class T>
  void updateMeshIterFace(p8est_iter_face_info_t *info, void *user_data);
  //template<unsigned int O,unsigned int MO>
  //void meshIterFacePartition(p8est_iter_face_info_t *info, void *user_data);
  //template<unsigned int O,unsigned int MO>
  //void p8estInitGeom(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
  //template<unsigned int O,unsigned int MO>
  //void p8estInitGeomPeriodicCube(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
  //template <unsigned int O,unsigned int MO>
  //void geomInit(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
  //template<unsigned int O,unsigned int MO>
  //void cfvHighOrderNodes (p8est_geometry_t * geom,
    //p4est_topidx_t which_tree,
    //const double abc[3],
    //double xyz[3]);

  //Derived class of <class distributedObjBase>
  template<class T>
  class dynamicMesh
    :public mesh
  {
    private:
      bool initialized;
      real3Vector highOrderNodes;
      std::vector<std::vector <llabel> > cellToHiOrderNodes;
    //May should be put into private section??????
    public:
      class p8estData
      {
        public:
        //At any stage, localIndex should be kept updated.
        llabel localIndex;
        glabel globalIndex;
        //Flag to tag the boundary type of a cell;
        //Used to mantain the boundary information more efficiently
        //int boundaryType[6];
        //int physicalId[6];
        int patchIdx[6];
        //Mapping information for old mesh
        //fineCoarsenTag == -2: new cell //Used at init
        //fineCoarsenTag ==  0: equal
        //fineCoarsenTag == -1: coarser
        //fineCoarsenTag == +1: finer
        struct map{
          int fineCoarsenTag;
          union
          {
            struct
            {
              llabel oldLocalIdx;
              glabel oldGlobalIdx;
            } equal;
            struct
            {
              llabel oldLocalIdx;
              glabel oldGlobalIdx;
              int childID;
            } finer;
            struct
            {
              llabel oldLocalIdx[8];
              glabel oldGlobalIdx[8];
            } coarser;
          } fcTag;
        } p;
      };
    private:
      struct invMap{
        int fineCoarsenTag;
        union
        {
          struct
          {
            llabel newLocalIdx;
          } equal;
          struct
          {
            llabel newLocalIdx[8];
            int childID[8];
          } finer;
          struct
          {
            llabel newLocalIdx;
          } coarser;
        } fcTag;
      };

      //Data section for p8est
      p8est_connectivity_t *connPtr_;
      bool withPeriodicBoundary_;
      p8est_t              *p8estPtr_;
      p8est_ghost_t        *p8estGhostPtr_;
      p8estData            *ghostData_;
      std::vector<typename p8estData::map>     mapVector_;
      std::vector<invMap>     invMapVector_;

      //If cell "A" or one of the neighboring cells of cell "A"
      //is a new cell (after refining and coarsening), the newCellFlag_[a] = 1;
      //This is used for PCG procedure;
      std::vector<int8_t>         newCellFlag_;
      std::vector<p8estData*>     p8estDataPtrVector_;

      //Used for partition fuctions
      std::vector<glabel> oldGlobalFirstQuadrant;
      std::vector<glabel> newGlobalFirstQuadrant;
      //These two varibles store the information for the old mesh
      std::vector<std::vector<llabel> > send_pattern_partition;
      std::vector<int> rank_for_cell_partition;

      llabel oldNCells_;
      llabel oldNGhostCells_;

      //Face element and volume element;
      //Using two sets of faceElems;
      //the point number of the faceElems1_ is k+1,
      //the point number of the faceElems2_ is floor(3k/2)+1,
      //k is the polynomial order
      //
      //If k+1 == floor(3k/2) +1, only one faceElems are used
      bool oversetFaceIntPt;
      std::vector<faceElement>            faceElems1_;
      std::vector<faceElement>            faceElems2_;
      std::vector<T>                      volElems_;

      //Old information used in updateMesh function
      faceVector oldfaces_;
      std::vector<faceElement> oldfaceElems1_;
      std::vector<faceElement> oldfaceElems2_;
      std::vector<T>           oldvolElems_;

      llabelVectorVector oldcefa_;

      //class used to store the structured hierachy of the mesh.
      //Since the matrices usually depends on the
      //relative position of neighboring cells, this class is used
      //to keep the relative positions unchanged.
      class neighbor{
        public:
          int finer;
          llabel index[P8EST_HALF];
          llabel pidx[P8EST_HALF];
      };
      struct ceceStruct{
        neighbor ne[P8EST_FACES];
      };
      std::vector<ceceStruct> p8estCece_;


    public:
      const std::vector<faceElement>& faceElems1(){return faceElems1_;}
      const std::vector<faceElement>& faceElems2(){
        if(oversetFaceIntPt) return faceElems2_;
        else return faceElems1_;
        }
      const std::vector<T>& volElems(){return volElems_;}
    public:
      //void restoreFromeFile(const std::string &fname, const std::string& meshfile, real refL);
      //Called after refineCoarsening or initMesh
      //Generate the map and invMap
      //Used for mapping data
      void updateMapVector();

      const std::vector<typename p8estData::map>& mapVector(){return    mapVector_;}
      const std::vector<p8estData*>& p8estDataVector(){return    p8estDataPtrVector_;}
      const std::vector<int8_t>& newCellFlag(){return    newCellFlag_;}

      //int refineCoarsen(std::vector<int>& flag, int minLevel, int maxLevel);
      //void checkGeom();

      //Partition funciton
      //Rebalance the mesh and redistribute the user defined cell based variables:
      //udReal, udHigh, and udPointer
      //The udSize defines the size of real varibles for each cell, and its length is n;
      //*T could not contain other pointers, if it has, the memory it points to would not be redistributed
      //void partitionMesh();

      //template<class TT>
      //void partitionUserData(int n, const std::vector<int>& udRealSize, std::vector<real> **udReal,
          //int p, std::vector<TT*> **udPointer
          //);
    public:
      enum predefinedMesh
      {
        periodicCube,
      };
    private:
      int initMeshPeriodicCube(int initLevel = 0, real refL = 1.0);
      real referenceLength;
    public:
      //Function needed by the constructor
      int initMesh(const std::string& fname, int initLevel = 0, real refL = 1.0);
      int initMesh(predefinedMesh t, int initLevel = 0, real refL = 1.0);
      //Constructor
      dynamicMesh(MPI_Comm mpiComm = MPI_COMM_WORLD);
      //Read mesh from file
      dynamicMesh(const std::string& fname, int initLevel = 0, real refL = 1.0, MPI_Comm mpiComm = MPI_COMM_WORLD);
      //Init mesh from predefined mesh type
      dynamicMesh(predefinedMesh configMesh, int initLevel = 0, real refL = 1.0, MPI_Comm mpiComm = MPI_COMM_WORLD);

      //Function used to reset mesh();
      void resetMesh();


      //Update the mesh information (except the node information)
      //including cece cefa faces.firstCell and faces.secondCell
      //faceElement information, etc.

      //The p8est, p8est_ghost_t and p8estData information must be
      //up-to-date
      void updateMesh();

      ~dynamicMesh();

      //void saveMesh(const std::string& filename);
      //void mapData (int n, std::vector<real>& data);
      friend void updateMeshIterVol <T> (p8est_iter_volume_info_t *info, void *user_data);
      friend void updateMeshIterFace<T> (p8est_iter_face_info_t *info, void *user_data);
      //friend void meshIterFacePartition<O,MO> (p8est_iter_face_info_t *info, void *user_data);
      //friend void p8estInitGeom<O,MO>(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
      //friend void p8estInitGeomPeriodicCube<O,MO>(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
      //friend void geomInit<O,MO>(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
      //friend void cfvHighOrderNodes<O,MO>(p8est_geometry_t * geom, p4est_topidx_t which_tree, const double abc[3], double xyz[3]);
  };

  //template<unsigned int O, unsigned int MO>
  //void meshIterVolAllocUserData(p8est_iter_volume_info_t *info, void *user_data)
  //{
    //p8est_quadrant_t *q = info->quad;
    //q->p.user_data = sc_mempool_alloc(info->p4est->user_data_pool);
    //p8estInitGeom<O,MO>(info->p4est, info->treeid, info->quad);
  //};

  //template<unsigned int O, unsigned int MO>
  //void cfvDynamicMesh<O, MO>::saveMesh( const std::string& filename)
  //{
    //std::string p8estname = filename+".p8est";

    //p8estPtr_->connectivity = connNPPtr_;
    //p8est_save(p8estname.c_str(), p8estPtr_, 0);
    //p8estPtr_->connectivity = connPtr_;

    //std::string coordname = filename+".coord";
    //std::string momentname = filename+".moment";
    //std::string patchname = filename+".patch";
    //int width = floor(log10(this->mpiSize()))+1;
    //if(!this->mpiRank()){
      //std::ofstream out;
      //out.open(coordname.c_str());
      //out << MO <<"\n";
      //out << this->mpiSize() << "\n";
      //for (int ii = 0; ii < this->mpiSize()+1; ++ii){
        //out << p8estPtr_->global_first_quadrant[ii] << "\n";
      //}
      //for (int ii = 0; ii < this->mpiSize(); ++ii){
        //std::stringstream cnameii;
        //cnameii << filename << ".coord" << std::setw(width) << std::setfill('0') << ii;
        //out << cnameii.str() << "\n";
      //}
      //out.close();

      //out.open(momentname.c_str());
      //out << O <<"\n";
      //out << this->mpiSize() << "\n";
      //for (int ii = 0; ii < this->mpiSize()+1; ++ii){
        //out << p8estPtr_->global_first_quadrant[ii] << "\n";
      //}
      //for (int ii = 0; ii < this->mpiSize(); ++ii){
        //std::stringstream mnameii;
        //mnameii << filename << ".moment" << std::setw(width) << std::setfill('0') << ii;
        //out << mnameii.str() << "\n";
      //}
      //out.close();

      //out.open(patchname.c_str());
      //out << this->mpiSize() << "\n";
      //for (int ii = 0; ii < this->mpiSize()+1; ++ii){
        //out << p8estPtr_->global_first_quadrant[ii] << "\n";
      //}
      //for (int ii = 0; ii < this->mpiSize(); ++ii){
        //std::stringstream pnameii;
        //pnameii << filename << ".patch" << std::setw(width) << std::setfill('0') << ii;
        //out << pnameii.str() << "\n";
      //}
      //out.close();
    //}
    //std::stringstream cname_this_rank;
    //cname_this_rank << filename << ".coord" << std::setw(width) << std::setfill('0') << this->mpiRank();
    //std::stringstream mname_this_rank;
    //mname_this_rank << filename << ".moment" << std::setw(width) << std::setfill('0') << this->mpiRank();
    //std::stringstream pname_this_rank;
    //pname_this_rank << filename << ".patch" << std::setw(width) << std::setfill('0') << this->mpiRank();
    //std::ofstream out;
    //out.open(cname_this_rank.str(), std::ios::out | std::ios::binary);
    //if(!out.good()){
      //printf("FILE open error in %s, line %d\n", __FILE__,__LINE__);
      //exit(-1);
    //}
    //for (int ii = 0; ii < this->nCells(); ++ii){
      //out.write( (char*)((p8estDataPtrVector_[ii])->coord), sizeof(p8estData::coord));
    //}
    //out.close();

    //out.open(mname_this_rank.str(), std::ios::out | std::ios::binary);
    //if(!out.good()){
      //printf("FILE open error in %s, line %d\n", __FILE__,__LINE__);
      //exit(-1);
    //}
    //for (int ii = 0; ii < this->nCells(); ++ii){
      //out.write((char*)(&(this->cellMoment_[ii])), sizeof(moment<O>));
    //}
    //out.close();

    //out.open(pname_this_rank.str(), std::ios::out);
    //if(!out.good()){
      //printf("FILE open error in %s, line %d\n", __FILE__,__LINE__);
      //exit(-1);
    //}
    //for (int ii = 0; ii < this->nFaces(); ++ii){
      //const face& f = this->faces()[ii];
      //if(f.secondCell==-1){
        //out << f.firstCell << " " << (int)f.firstSide << " " << f.patchIdx << "\n";
      //}
    //}
    //out.close();
  //}

  //template<unsigned int O, unsigned int MO>
  //void p8estInitGeomPeriodicCube(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q)
  //{
    //const double pi = atan(1.0L)*4;
    //cfvDynamicMesh<O,MO> *mesh = (cfvDynamicMesh<O,MO>*)p8est->user_pointer;
    //typename cfvDynamicMesh<O,MO>::p8estData *data
      //= (typename cfvDynamicMesh<O, MO>::p8estData*)q->p.user_data;

    //for (int ii = 0; ii < P8EST_FACES; ++ii){
      //data->patchIdx[ii] = -1;
      //data->periodicIdx[ii] = -1;
    //}
    //data->p.fineCoarsenTag = -2;//New cell
    //real3 xyz_root[8];
    //real3 off(0,0,0);
    //for (int ii = 0; ii < 6; ++ii){
      //data->faceOffsets[ii] = off;
      //if(q->x==0){
        //data->faceOffsets[0].x = +1;
        //data->periodicIdx[0] = 0;
      //}
      //if(q->y==0){
        //data->faceOffsets[2].y = +1;
        //data->periodicIdx[2] = 2;
      //}
      //if(q->z==0){
        //data->faceOffsets[4].z = +1;
        //data->periodicIdx[4] = 4;
      //}
      //if(q->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(q->level)){
        //data->faceOffsets[1].x = -1;
        //data->periodicIdx[1] = 1;
      //}
      //if(q->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(q->level)){
        //data->faceOffsets[3].y = -1;
        //data->periodicIdx[3] = 3;
      //}
      //if(q->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(q->level)){
        //data->faceOffsets[5].z = -1;
        //data->periodicIdx[5] = 5;
      //}
    //}
    //for (int kk = 0; kk < 2; ++kk){
      //for (int jj = 0; jj < 2; ++jj){
        //for (int ii = 0; ii < 2; ++ii){
          //p4est_topidx_t v = mesh->connPtr_->tree_to_vertex[8*which_tree+4*kk+2*jj+ii];
          //xyz_root[kk*4+2*jj+ii].x = mesh->connPtr_->vertices[3*v  ];
          //xyz_root[kk*4+2*jj+ii].y = mesh->connPtr_->vertices[3*v+1];
          //xyz_root[kk*4+2*jj+ii].z = mesh->connPtr_->vertices[3*v+2];
        //}
      //}
    //}

    //real3 xyz[8];
    //double intsize = 1.0/P8EST_ROOT_LEN;
    //double eta_x[2];
    //double eta_y[2];
    //double eta_z[2];
    //eta_x[0] = 2*(q->x*intsize)-1;
    //eta_y[0] = 2*(q->y*intsize)-1;
    //eta_z[0] = 2*(q->z*intsize)-1;

    //eta_x[1] = 2*(q->x+P8EST_QUADRANT_LEN(q->level))*intsize-1;
    //eta_y[1] = 2*(q->y+P8EST_QUADRANT_LEN(q->level))*intsize-1;
    //eta_z[1] = 2*(q->z+P8EST_QUADRANT_LEN(q->level))*intsize-1;

    //for (int kk = 0; kk < MO+1; ++kk){
      //for (int jj = 0; jj < MO+1; ++jj){
        //for (int ii = 0; ii < MO+1; ++ii){
          //double xf = (ii+0.0)/MO*(eta_x[1]-eta_x[0])+eta_x[0];
          //double yf = (jj+0.0)/MO*(eta_y[1]-eta_y[0])+eta_y[0];
          //double zf = (kk+0.0)/MO*(eta_z[1]-eta_z[0])+eta_z[0];
          //data->coord[kk*(MO+1)*(MO+1)+jj*(MO+1)+ii] =
            //geom::getCoordAtPointHex<1>(xyz_root, real3(xf,yf,zf));
        //}
      //}
    //}
  //}


  //template<unsigned int O, unsigned int MO>
  //void cfvDynamicMesh<O, MO>::initP8estPeriodicCube(int initLevel, real refL)
  //{
    //const double pi = atan(1.0L)*4;
    //connNPPtr_       = p8est_connectivity_new_unitcube();
    //connPtr_       = p8est_connectivity_new_periodic();
    //const double xlen = 1;
    //const double ylen = 1;
    //const double zlen = 1;
    //const double        vertices[8 * 3] = {
      //0, 0, 0,
      //xlen, 0, 0,
      //0, ylen, 0,
      //xlen, zlen, 0,
      //0, 0, zlen,
      //xlen, 0, zlen,
      //0, ylen, zlen,
      //xlen, ylen, zlen,
    //};

    //for(int ii = 0; ii < 24; ++ii) connPtr_->vertices[ii] = vertices[ii];
    //for (int ii = 0; ii < 3*connPtr_->num_vertices; ++ii){
      ////connPtr_->vertices[ii] = (connPtr_->vertices[ii]-0.5)*2*pi;
      //connPtr_->vertices[ii]/=refL;
      ////printf("%6.3e \n",connPtr_->vertices[ii]);
    //}
    ////connPtr_       = p8est_connectivity_new_unitcube();
    //highOrderMesh<O>::meshOrder_ = MO;
    //highOrderMesh<O>::highOrderNodes = new real[(MO+1)*(MO+1)*(MO+1)*3];
    //for(int kk = 0; kk < MO+1; ++kk){
      //for (int jj = 0; jj < MO+1; ++jj){
        //for (int ii = 0; ii < MO+1; ++ii){
          //int v = (MO+1)*(MO+1)*kk + (MO+1)*jj + ii;
          //this->highOrderNodes[3*v+0] = ((ii+0.0)/MO)*xlen;
          //this->highOrderNodes[3*v+1] = ((jj+0.0)/MO)*ylen;
          //this->highOrderNodes[3*v+2] = ((kk+0.0)/MO)*zlen;
        //}
      //}
    //}
    //highOrderMesh<O>::cellToHiOrderNodes.resize(1);
    //highOrderMesh<O>::cellToHiOrderNodes[0].resize((MO+1)*(MO+1)*(MO+1));
    //for (int ii = 0; ii < (MO+1)*(MO+1)*(MO+1); ++ii){
      //highOrderMesh<O>::cellToHiOrderNodes[0][ii] = ii;
    //}
    //cfvMesh<O>::physicalId_.resize(0);
    //cfvMesh<O>::periodicType_.assign(6,1);
    //cfvMesh<O>::periodicValue_.resize(6);
    //for(int ii = 0; ii < 6; ++ii){
      //this->periodicValue_[ii].rotateCenter = real3(0,0,0);
      //this->periodicValue_[ii].rotateAngle  = real3(0,0,0);
    //}
    //this->periodicValue_[0].translation  = real3(xlen,0,0);
    //this->periodicValue_[1].translation  = real3(-xlen,0,0);
    //this->periodicValue_[2].translation  = real3(0,ylen,0);
    //this->periodicValue_[3].translation  = real3(0,-ylen,0);
    //this->periodicValue_[4].translation  = real3(0,0,zlen);
    //this->periodicValue_[5].translation  = real3(0,0,-zlen);
    //cfvMesh<O>::boundaryList_.resize(0);
    ////cfvMesh<O>::physicalId_[0] = 1;

    //p8estPtr_      = p8est_new_ext(primitiveMesh::mpiComm_, connPtr_, 0,
        //initLevel, 1, sizeof(p8estData), p8estInitGeomPeriodicCube<O,MO>, this);
    //MPI_Barrier(primitiveMesh::mpiComm_);
    //p8estGhostPtr_ = p8est_ghost_new(p8estPtr_, P8EST_CONNECT_FACE);
    //MPI_Barrier(primitiveMesh::mpiComm_);
    //ghostData_     = P4EST_ALLOC(p8estData, 2*p8estGhostPtr_->ghosts.elem_count);
    //p8est_ghost_exchange_data(p8estPtr_, p8estGhostPtr_, ghostData_);
  //}

  //template<unsigned int O,unsigned int MO>
  //void p8estInitGeom(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q)
  //{
    //const double pi = atan(1.0L)*4;
    //cfvDynamicMesh<O,MO> *mesh = (cfvDynamicMesh<O,MO>*)p8est->user_pointer;
    //typename cfvDynamicMesh<O,MO>::p8estData *data
      //= (typename cfvDynamicMesh<O,MO>::p8estData*)q->p.user_data;

    //for (int ii = 0; ii < P8EST_FACES; ++ii){
      //data->patchIdx[ii] = -1;
      //data->periodicIdx[ii] = -1;
    //}
    //data->p.fineCoarsenTag = -2;//New cell
    //real3* xyz_root;
    //real3 off(0,0,0);
    //for (int ii = 0; ii < 6; ++ii){
      //data->faceOffsets[ii] = off;
    //}
    //int numHiONodesPerCell = (MO+1)*(MO+1)*(MO+1);
    //xyz_root = new real3[numHiONodesPerCell];
    //for (int kk = 0; kk < MO+1; ++kk){
      //for (int jj = 0; jj < MO+1; ++jj){
        //for (int ii = 0; ii < MO+1; ++ii){
          //if(MO==1){
            //p4est_topidx_t v = mesh->connPtr_->tree_to_vertex[8*which_tree+4*kk+2*jj+ii];
            //xyz_root[kk*(MO+1)*(MO+1)+(MO+1)*jj+ii].x = mesh->connPtr_->vertices[3*v  ];
            //xyz_root[kk*(MO+1)*(MO+1)+(MO+1)*jj+ii].y = mesh->connPtr_->vertices[3*v+1];
            //xyz_root[kk*(MO+1)*(MO+1)+(MO+1)*jj+ii].z = mesh->connPtr_->vertices[3*v+2];
          //}else{
            //long int v = mesh->cellToHiOrderNodes[which_tree][(MO+1)*(MO+1)*kk+(MO+1)*jj+ii];
            //xyz_root[kk*(MO+1)*(MO+1)+(MO+1)*jj+ii].x = mesh->highOrderNodes[3*v  ];
            //xyz_root[kk*(MO+1)*(MO+1)+(MO+1)*jj+ii].y = mesh->highOrderNodes[3*v+1];
            //xyz_root[kk*(MO+1)*(MO+1)+(MO+1)*jj+ii].z = mesh->highOrderNodes[3*v+2];
          //}
        //}
      //}
    //}

    //real3 xyz[8];
    //double intsize = 1.0/P8EST_ROOT_LEN;
    //double eta_x[2];
    //double eta_y[2];
    //double eta_z[2];
    //eta_x[0] = 2*(q->x*intsize)-1;
    //eta_y[0] = 2*(q->y*intsize)-1;
    //eta_z[0] = 2*(q->z*intsize)-1;

    //eta_x[1] = 2*(q->x+P8EST_QUADRANT_LEN(q->level))*intsize-1;
    //eta_y[1] = 2*(q->y+P8EST_QUADRANT_LEN(q->level))*intsize-1;
    //eta_z[1] = 2*(q->z+P8EST_QUADRANT_LEN(q->level))*intsize-1;

    //for (int kk = 0; kk < MO+1; ++kk){
      //for (int jj = 0; jj < MO+1; ++jj){
        //for (int ii = 0; ii < MO+1; ++ii){
          //double xf = (ii+0.0)/MO*(eta_x[1]-eta_x[0])+eta_x[0];
          //double yf = (jj+0.0)/MO*(eta_y[1]-eta_y[0])+eta_y[0];
          //double zf = (kk+0.0)/MO*(eta_z[1]-eta_z[0])+eta_z[0];
          //data->coord[kk*(MO+1)*(MO+1)+jj*(MO+1)+ii] =
            //geom::getCoordAtPointHex(MO, xyz_root, real3(xf,yf,zf));
        //}
      //}
    //}
    //delete[] xyz_root;
  //}

  //template<unsigned int O, unsigned int MO>
  //void cfvDynamicMesh<O, MO>::initP8est(const std::string& fname, int initLevel, real refL)
  //{
    //std::string p8estFileName = fname+".p8estConn";
    //std::string p8estNPFileName = fname+".p8estConnNP";
    //size_t connSize;
    //if(!this->mpiRank())std::cout << "Loading p8est connectivity file" << std::endl;
    //connPtr_         = p8est_connectivity_load(p8estFileName.c_str(),&connSize);
    //connNPPtr_       = p8est_connectivity_load(p8estNPFileName.c_str(),&connSize);
    //if(!this->mpiRank())std::cout << "End loading p8est connectivity file" << std::endl;
    //for (int ii = 0; ii < 3*connPtr_->num_vertices; ++ii){
      //connPtr_->vertices[ii]/=refL;
    //}

    ////for (int jj = 0; 196+4*jj+1<296;++jj){
      ////int idx1 = 196+4*jj+0;
      ////int idx2 = 196+4*jj+1;
      ////connPtr_->vertices[idx1*3+0] = -0.4177786794919244;
      ////connPtr_->vertices[idx2*3+0] = -0.4177786794919244;
    ////}
    ////const double pi = atan(1.0)*4;
    ////real theta = 1*pi/4;
////#ifdef CIRCLE_VISCOUS
    ////for (int ii = 0; ii < connPtr_->num_vertices; ++ii){
      ////connPtr_->vertices[3*ii+2] *= 0.01;
      ////connNPPtr_->vertices[3*ii+2] *= 0.01;
    ////}
////#endif

    //std::string hiONodesFileName = fname+".hiONodes";
    //std::ifstream hiONodesStream;
    //hiONodesStream.open(hiONodesFileName.c_str());
    //int meshOrder;
    //hiONodesStream >> meshOrder;
    //highOrderMesh<O>::meshOrder_ = MO;
    //if(MO > 1 && MO != meshOrder){
      //std::cerr << "The mesh file order does not coordinate with the declared cfvDynamicMesh class" << std::endl;
      //exit(-1);
    //}

    //int numHiONodesPerCell = pow(highOrderMesh<O>::meshOrder_+1,3);
    //if(MO==1)
      //highOrderMesh<O>::highOrderNodes = nullptr;
    //else{
      //long int num_vertices;
      //int num_trees;
      //hiONodesStream >> num_vertices;
      //hiONodesStream >> num_trees;
      //highOrderMesh<O>::highOrderNodes = new real[num_vertices*3];
      //for (int ii = 0; ii < num_vertices*3; ++ii){
        //double temp;
        //hiONodesStream >> temp;
        //highOrderMesh<O>::highOrderNodes[ii]=temp/refL;
      //}
////#ifdef CIRCLE_VISCOUS
      ////for (int ii = 0; ii < num_vertices; ++ii){
        ////highOrderMesh<O>::highOrderNodes[3*ii+2]*=0.01;
      ////}
////#endif
      //highOrderMesh<O>::cellToHiOrderNodes.resize(num_trees);
      //for (int ii = 0; ii < num_trees; ++ii){
        //highOrderMesh<O>::cellToHiOrderNodes[ii].resize(numHiONodesPerCell);
        //long int idxTemp;
        //for (int jj = 0; jj < numHiONodesPerCell; ++jj){
          //hiONodesStream >> idxTemp;
          //highOrderMesh<O>::cellToHiOrderNodes[ii][jj] = idxTemp;
        //}
      //}
    //}

    //hiONodesStream.close();
    //std::string   bcFileName = fname+".bc";
    //std::ifstream bcFileStream;
    //bcFileStream.open(bcFileName.c_str());
    //int numPatches;
    //int numPeriodic;
    //bcFileStream >> numPatches;
    //for(int ii = 0; ii < numPatches; ++ii){
      //int physicalId;
      //std::string name("");
      //bcFileStream >> physicalId;
      //char ch;
      //bcFileStream >> std::noskipws;
      //do{
        //bcFileStream >> ch;
      //}while(ch!='"');
      //while(1){
        //bcFileStream >> ch;
        //if(ch!='"') name += ch;
        //else break;
      //}
      //bcFileStream >> std::skipws;
    //}
    //cfvMesh<O>::physicalId_.resize(numPatches);
    ////Just init, the detailed boundaryType is processed by the specified PDE solver
    //cfvMesh<O>::boundaryList_.resize(numPatches);
    //std::vector< std::vector<std::pair<int,int> > > boundary_list;
    //boundary_list.resize(numPatches);
    //for (int ii = 0; ii < numPatches; ++ii){
      //int patchIdx;
      //int numbcs;
      //bcFileStream >> patchIdx >> cfvMesh<O>::physicalId_[ii] >> numbcs;
      //boundary_list[ii].resize(numbcs);
      //for (int jj = 0; jj < numbcs; ++jj){
        //bcFileStream >> boundary_list[ii][jj].first >> boundary_list[ii][jj].second;
      //}
    //}
    //bcFileStream >> numPeriodic;
    //cfvMesh<O>::periodicType_.resize(numPeriodic);
    //cfvMesh<O>::periodicValue_.resize(numPeriodic);
    //std::vector< std::vector<std::pair<int,int> > > periodic_list;
    //periodic_list.resize(numPeriodic);
    //for(int ii = 0; ii < numPeriodic; ++ii){
      //int id;
      //real3 rc;
      //real3 ra;
      //real3 trans;
      //int numPeriodicface;
      //bcFileStream >> id >> cfvMesh<O>::periodicType_[ii]
        //>> rc.x >> rc.y >> rc.z >> ra.x >> ra.y >> ra.z >> trans.x >> trans.y >> trans.z >> numPeriodicface;
////#ifdef CIRCLE_VISCOUS
      ////trans.z*=0.01;
////#endif
      //this->periodicValue_[ii].rotateCenter = rc;
      //this->periodicValue_[ii].rotateAngle = ra;
      //this->periodicValue_[ii].translation = trans;
      //periodic_list[ii].resize(numPeriodicface);
      //for (int jj = 0; jj < numPeriodicface; ++jj){
        //bcFileStream >> periodic_list[ii][jj].first >> periodic_list[ii][jj].second;
      //}
    //}
    //bcFileStream.close();

    //p8estPtr_      = p8est_new_ext(primitiveMesh::mpiComm_, connPtr_, 0, initLevel, 1, sizeof(p8estData), p8estInitGeom<O,MO>, this);
    //p4est_topidx_t      t;
    //p4est_topidx_t      first_local_tree = p8estPtr_->first_local_tree;
    //p4est_topidx_t      last_local_tree = p8estPtr_->last_local_tree;
    //sc_array_t         *trees = p8estPtr_->trees;
    //p8est_tree_t       *tree;
    //p8est_quadrant_t   *quad;
    //p8estData          *ud;
    //size_t              si, n_quads;
    //sc_array_t         *quadrants;

    //for (int ii = 0; ii < numPatches; ++ii){
      //for (int jj = 0; jj < boundary_list[ii].size(); ++jj){
        //int t = boundary_list[ii][jj].first;
        //int which_face = boundary_list[ii][jj].second;
        //tree = p8est_tree_array_index (trees, t);
        //quadrants = &(tree->quadrants);
        //n_quads = quadrants->elem_count;
        //for (si = 0; si < n_quads; si++) {
          //quad = p8est_quadrant_array_index (quadrants, si);
          ////Set the physical Id and offset of this quadrant
          //ud = (typename cfvDynamicMesh<O,MO>::p8estData*)quad->p.user_data;
          //if(quad->x==0&&which_face==0)
            //ud->patchIdx[0] = ii;
          //if(quad->y==0&&which_face==2)
            //ud->patchIdx[2] = ii;
          //if(quad->z==0&&which_face==4)
            //ud->patchIdx[4] = ii;
          //if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
            //ud->patchIdx[1] = ii;
          //if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
            //ud->patchIdx[3] = ii;
          //if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
            //ud->patchIdx[5] = ii;
        //}
      //}
    //}

    //for (int ii = 0; ii < numPeriodic; ++ii){
      //for (int jj = 0; jj < periodic_list[ii].size(); ++jj){
        //int t = periodic_list[ii][jj].first;
        //int which_face = periodic_list[ii][jj].second;
        //tree = p8est_tree_array_index (trees, t);
        //quadrants = &(tree->quadrants);
        //n_quads = quadrants->elem_count;
        //for (si = 0; si < n_quads; si++) {
          //quad = p8est_quadrant_array_index (quadrants, si);
          ////Set the physical Id and offset of this quadrant
          //ud = (typename cfvDynamicMesh<O,MO>::p8estData*)quad->p.user_data;
          //if(quad->x==0&&which_face==0)
            //ud->periodicIdx[0] = ii;
          //if(quad->y==0&&which_face==2)
            //ud->periodicIdx[2] = ii;
          //if(quad->z==0&&which_face==4)
            //ud->periodicIdx[4] = ii;
          //if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
            //ud->periodicIdx[1] = ii;
          //if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
            //ud->periodicIdx[3] = ii;
          //if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
            //ud->periodicIdx[5] = ii;
        //}
      //}
    //}
    //MPI_Barrier(primitiveMesh::mpiComm_);
    //p8estGhostPtr_ = p8est_ghost_new(p8estPtr_, P8EST_CONNECT_FACE);
    //MPI_Barrier(primitiveMesh::mpiComm_);
    //ghostData_     = P4EST_ALLOC(p8estData, 2*p8estGhostPtr_->ghosts.elem_count);
    //p8est_ghost_exchange_data(p8estPtr_, p8estGhostPtr_, ghostData_);
  //}

  //template<unsigned int O, unsigned int MO>
  //void cfvDynamicMesh<O, MO>:: restoreFromeFile(const std::string& fname, const std::string& meshFile, real refL)
  //{
    //std::string p8estFileName = meshFile+".p8estConn";
    //std::string p8estNPFileName = meshFile+".p8estConnNP";
    //size_t connSize;
    //if(!this->mpiRank())std::cout << "Loading p8est connectivity file" << std::endl;
    //connPtr_         = p8est_connectivity_load(p8estFileName.c_str(),&connSize);
    ////connNPPtr_       = p8est_connectivity_load(p8estNPFileName.c_str(),&connSize);
    //if(!this->mpiRank())std::cout << "End loading p8est connectivity file" << std::endl;
    //for (int ii = 0; ii < 3*connPtr_->num_vertices; ++ii){
      //connPtr_->vertices[ii]/=refL;
    //}
    //std::string hiONodesFileName = meshFile+".hiONodes";
    //std::ifstream hiONodesStream;
    //hiONodesStream.open(hiONodesFileName.c_str());
    //int meshOrder;
    //hiONodesStream >> meshOrder;
    //highOrderMesh<O>::meshOrder_ = MO;
    //if(MO > 1 && MO != meshOrder){
      //std::cerr << "The mesh file order does not coordinate with the declared cfvDynamicMesh class" << std::endl;
      //exit(-1);
    //}

    //int numHiONodesPerCell = pow(highOrderMesh<O>::meshOrder_+1,3);
    //if(MO==1)
      //highOrderMesh<O>::highOrderNodes = nullptr;
    //else{
      //long int num_vertices;
      //int num_trees;
      //hiONodesStream >> num_vertices;
      //hiONodesStream >> num_trees;
      //highOrderMesh<O>::highOrderNodes = new real[num_vertices*3];
      //for (int ii = 0; ii < num_vertices*3; ++ii){
        //double temp;
        //hiONodesStream >> temp;
        //highOrderMesh<O>::highOrderNodes[ii]=temp/refL;
      //}
      //highOrderMesh<O>::cellToHiOrderNodes.resize(num_trees);
      //for (int ii = 0; ii < num_trees; ++ii){
        //highOrderMesh<O>::cellToHiOrderNodes[ii].resize(numHiONodesPerCell);
        //long int idxTemp;
        //for (int jj = 0; jj < numHiONodesPerCell; ++jj){
          //hiONodesStream >> idxTemp;
          //highOrderMesh<O>::cellToHiOrderNodes[ii][jj] = idxTemp;
        //}
      //}
    //}

    //hiONodesStream.close();
    //std::string   bcFileName = meshFile+".bc";
    //std::ifstream bcFileStream;
    //bcFileStream.open(bcFileName.c_str());
    //int numPatches;
    //int numPeriodic;
    //bcFileStream >> numPatches;
    //for(int ii = 0; ii < numPatches; ++ii){
      //int physicalId;
      //std::string name("");
      //bcFileStream >> physicalId;
      //char ch;
      //bcFileStream >> std::noskipws;
      //do{
        //bcFileStream >> ch;
      //}while(ch!='"');
      //while(1){
        //bcFileStream >> ch;
        //if(ch!='"') name += ch;
        //else break;
      //}
      //bcFileStream >> std::skipws;
    //}
    //cfvMesh<O>::physicalId_.resize(numPatches);
    ////Just init, the detailed boundaryType is processed by the specified PDE solver
    //cfvMesh<O>::boundaryList_.resize(numPatches);
    //std::vector< std::vector<std::pair<int,int> > > boundary_list;
    //boundary_list.resize(numPatches);
    //for (int ii = 0; ii < numPatches; ++ii){
      //int patchIdx;
      //int numbcs;
      //bcFileStream >> patchIdx >> cfvMesh<O>::physicalId_[ii] >> numbcs;
      //boundary_list[ii].resize(numbcs);
      //for (int jj = 0; jj < numbcs; ++jj){
        //bcFileStream >> boundary_list[ii][jj].first >> boundary_list[ii][jj].second;
      //}
    //}
    //bcFileStream >> numPeriodic;
    //cfvMesh<O>::periodicType_.resize(numPeriodic);
    //cfvMesh<O>::periodicValue_.resize(numPeriodic);
    //std::vector< std::vector<std::pair<int,int> > > periodic_list;
    //periodic_list.resize(numPeriodic);
    //for(int ii = 0; ii < numPeriodic; ++ii){
      //int id;
      //real3 rc;
      //real3 ra;
      //real3 trans;
      //int numPeriodicface;
      //bcFileStream >> id >> cfvMesh<O>::periodicType_[ii]
        //>> rc.x >> rc.y >> rc.z >> ra.x >> ra.y >> ra.z >> trans.x >> trans.y >> trans.z >> numPeriodicface;
      //this->periodicValue_[ii].rotateCenter = rc;
      //this->periodicValue_[ii].rotateAngle = ra;
      //this->periodicValue_[ii].translation = trans;
      //periodic_list[ii].resize(numPeriodicface);
      //for (int jj = 0; jj < numPeriodicface; ++jj){
        //bcFileStream >> periodic_list[ii][jj].first >> periodic_list[ii][jj].second;
      //}
    //}
    //bcFileStream.close();

    //std::string p8estName = fname + ".p8est";
    //p8estPtr_ = p8est_load_ext(p8estName.c_str(), primitiveMesh::mpiComm_, sizeof(p8estData), 0, 1, 0, nullptr, &connNPPtr_);
    //p8estPtr_->connectivity = connPtr_;
    //p8estPtr_->user_pointer = this;
    ////Init the p8estData;
    //p8estPtr_->data_size = sizeof(p8estData);
    //p8estPtr_->user_data_pool = sc_mempool_new(sizeof(p8estData));
    ////Allocate the user data for the p8est
    //p8est_iterate(p8estPtr_, nullptr, nullptr, meshIterVolAllocUserData<O,MO>, nullptr, nullptr, nullptr);

    //p4est_topidx_t      t;
    //p4est_topidx_t      first_local_tree = p8estPtr_->first_local_tree;
    //p4est_topidx_t      last_local_tree = p8estPtr_->last_local_tree;
    //sc_array_t         *trees = p8estPtr_->trees;
    //p8est_tree_t       *tree;
    //p8est_quadrant_t   *quad;
    //p8estData          *ud;
    //size_t              si, n_quads;
    //sc_array_t         *quadrants;

    //for (int ii = 0; ii < numPatches; ++ii){
      //for (int jj = 0; jj < boundary_list[ii].size(); ++jj){
        //int t = boundary_list[ii][jj].first;
        //int which_face = boundary_list[ii][jj].second;
        //tree = p8est_tree_array_index (trees, t);
        //quadrants = &(tree->quadrants);
        //n_quads = quadrants->elem_count;
        //for (si = 0; si < n_quads; si++) {
          //quad = p8est_quadrant_array_index (quadrants, si);
          ////Set the physical Id and offset of this quadrant
          //ud = (typename cfvDynamicMesh<O,MO>::p8estData*)quad->p.user_data;
          //if(quad->x==0&&which_face==0)
            //ud->patchIdx[0] = ii;
          //if(quad->y==0&&which_face==2)
            //ud->patchIdx[2] = ii;
          //if(quad->z==0&&which_face==4)
            //ud->patchIdx[4] = ii;
          //if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
            //ud->patchIdx[1] = ii;
          //if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
            //ud->patchIdx[3] = ii;
          //if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
            //ud->patchIdx[5] = ii;
        //}
      //}
    //}

    //for (int ii = 0; ii < numPeriodic; ++ii){
      //for (int jj = 0; jj < periodic_list[ii].size(); ++jj){
        //int t = periodic_list[ii][jj].first;
        //int which_face = periodic_list[ii][jj].second;
        //tree = p8est_tree_array_index (trees, t);
        //quadrants = &(tree->quadrants);
        //n_quads = quadrants->elem_count;
        //for (si = 0; si < n_quads; si++) {
          //quad = p8est_quadrant_array_index (quadrants, si);
          ////Set the physical Id and offset of this quadrant
          //ud = (typename cfvDynamicMesh<O,MO>::p8estData*)quad->p.user_data;
          //if(quad->x==0&&which_face==0)
            //ud->periodicIdx[0] = ii;
          //if(quad->y==0&&which_face==2)
            //ud->periodicIdx[2] = ii;
          //if(quad->z==0&&which_face==4)
            //ud->periodicIdx[4] = ii;
          //if(quad->x==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==1)
            //ud->periodicIdx[1] = ii;
          //if(quad->y==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==3)
            //ud->periodicIdx[3] = ii;
          //if(quad->z==P8EST_ROOT_LEN-P8EST_QUADRANT_LEN(quad->level)&&which_face==5)
            //ud->periodicIdx[5] = ii;
        //}
      //}
    //}
    //MPI_Barrier(primitiveMesh::mpiComm_);
    //p8estGhostPtr_ = p8est_ghost_new(p8estPtr_, P8EST_CONNECT_FACE);
    //MPI_Barrier(primitiveMesh::mpiComm_);
    //ghostData_     = P4EST_ALLOC(p8estData, 2*p8estGhostPtr_->ghosts.elem_count);
    //p8est_ghost_exchange_data(p8estPtr_, p8estGhostPtr_, ghostData_);
  //}
  //template<unsigned int O, unsigned int MO>
  //void cfvDynamicMesh<O, MO>::updateP8estDataPtrVector()
  //{
    //label nCells = p8estPtr_->local_num_quadrants;
    //label nGhostCells = p8estGhostPtr_->ghosts.elem_count;
    //p8estDataPtrVector_.resize(nCells+nGhostCells);
    //int count = 0;

    //for (int ii = p8estPtr_->first_local_tree; ii <= p8estPtr_->last_local_tree; ++ii){
      //p8est_tree_t *tree = p8est_tree_array_index(p8estPtr_->trees, ii);
      //for (unsigned int jj = 0; jj < tree->quadrants.elem_count; ++jj){
        //p8est_quadrant_t *q = p8est_quadrant_array_index(&(tree->quadrants), jj);
        //p8estData *data = (p8estData*)q->p.user_data;
        //p8estDataPtrVector_[count] = data;
        //++count;
      //}
    //}
    //for (int ii = 0; ii < p8estGhostPtr_->ghosts.elem_count; ++ii){
      //p8estData *data = ghostData_ + ii;
      //p8estDataPtrVector_[ii+nCells] = data;
    //}
    //MPI_Barrier(primitiveMesh::mpiComm_);
  //}

  //template<unsigned int O, unsigned int MO>
  //void cfvDynamicMesh<O, MO>::updateMapVector()
  //{
    //label nCells = p8estPtr_->local_num_quadrants;
    //label nGhostCells = p8estGhostPtr_->ghosts.elem_count;
    //mapVector_.resize(nCells+nGhostCells);
    //int count = 0;
    //invMapVector_.resize(oldNCells_);
    //std::vector<std::vector<glabel> > temp1(oldNCells_);
    //std::vector<std::vector<int> > temp2(oldNCells_);
    //std::vector<int> tag(oldNCells_);

    //for (int ii = p8estPtr_->first_local_tree; ii <= p8estPtr_->last_local_tree; ++ii){
      //p8est_tree_t *tree = p8est_tree_array_index(p8estPtr_->trees, ii);
      //for (unsigned int jj = 0; jj < tree->quadrants.elem_count; ++jj){
        //p8est_quadrant_t *q = p8est_quadrant_array_index(&(tree->quadrants), jj);
        //p8estData *data = (p8estData*)q->p.user_data;
        //mapVector_[count] = data->p;
        //if(mapVector_[count].fineCoarsenTag==0){
          //glabel oldLocalIndex = mapVector_[count].fcTag.equal.oldLocalIndex;
          //temp1[oldLocalIndex] .push_back(count);
          //tag[oldLocalIndex] = 0;
        //}else if(mapVector_[count].fineCoarsenTag==-1){
          //for (int h = 0; h < P8EST_CHILDREN; ++h){
            //glabel oldLocalIndex = mapVector_[count].fcTag.coarser.oldLocalIndices[h];
            //temp1[oldLocalIndex].push_back(count);
            //tag[oldLocalIndex] = -1;
          //}
        //}else if(mapVector_[count].fineCoarsenTag==+1){
           //glabel oldLocalIndex = mapVector_[count].fcTag.finer.oldLocalIndex;
           //int    cid           = mapVector_[count].fcTag.finer.childID;
           //temp1[oldLocalIndex].push_back(count);
           //temp2[oldLocalIndex].push_back(cid);
           //tag[oldLocalIndex] = +1;
        //}
        //++count;
      //}
    //}
    //for (int ii = 0; ii < this->oldNCells_; ++ii){
      //if(tag[ii] == 0){
        //invMapVector_[ii].fineCoarsenTag = 0;
        //invMapVector_[ii].fcTag.equal.newLocalIndex = temp1[ii][0];
      //}else if(tag[ii] == 1){
        //invMapVector_[ii].fineCoarsenTag = 1;
        //assert(temp1[ii].size()==P8EST_CHILDREN);
        //for (int h = 0; h < P8EST_CHILDREN; ++h){
          //invMapVector_[ii].fcTag.finer.newLocalIndices[h] = temp1[ii][h];
          //invMapVector_[ii].fcTag.finer.childID[h] = temp2[ii][h];
        //}
      //}else if(tag[ii] == -1){
        //invMapVector_[ii].fineCoarsenTag = -1;
        //invMapVector_[ii].fcTag.coarser.newLocalIndex = temp1[ii][0];
      //}
    //}
    //for (int ii = 0; ii < p8estGhostPtr_->ghosts.elem_count; ++ii){
      //p8estData *data = ghostData_ + ii;
      //mapVector_[ii+nCells] = data->p;
    //}
    //MPI_Barrier(primitiveMesh::mpiComm_);
  //}
}//End for namespace fps

#include "./dynamicMesh_constructor.h"
#include "./dynamicMesh_updateMesh.h"
#include "./dynamicMesh_updateMapVector.h"
#endif
