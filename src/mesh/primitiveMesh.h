#ifndef FPS_MESH_PRIMITIVEMESH_H
#define FPS_MESH_PRIMITIVEMESH_H
#include <vector>

#include <config.inc>
#include <distributedObjBase.h>
#include <dataType.h>
#include <real3.h>
#include <face.h>
#include <cell.h>

#include <mpi.h>

namespace fps
{
  class primitiveMesh:
      public distributedObjBase
  {
    protected:
    //Primitive size data
      //nPoints only refer to the linear points, i.e. the vertices
      llabel nPoints_;
      llabel nFaces_;
      llabel nCells_;
      //Local points
      real3Vector  points_;


      //MPI_INFORMATION
      //Ghost cells
      //Information needed for FV computations
      llabel nGhostCells_;
      //Local + ghost cell's global indices
      glabelVector  cellGlobalIdx_;

      //Ghost cells' ranks
      std::vector<int>  ghostCellMPIRank_;

      //Send Pattern
      //sendPatternPtr->size() == mpiSize
      //The j-th cell needed to be sent to rank i is:
      //  sendPatternPtr[i][j](Local index)
      llabelVectorVector sendPattern_;

      //Receive Pattern
      //receivePatternPtr->size() == mpiSize
      //The j-th cell ready to be received from rank i is:
      //  reveivePatternPtr[i][j](global index)
      glabelVectorVector receivePattern_;

    //Connectivity
    //The physical id of a face is stored in the <class face>.
    //And the detailed boundary conditions are defined by the solver
    //including the periodic boundary conditions

      //Cell-cells
      //If index = (cece_)[i][j] >= nCells_, index refers to a ghost cell;
      //cece_.size() == nCells.
      //if cece_[i][j] < 0; then -cece_[i][j] is the patch index;
      llabelVectorVector cece_;
      //Cell-faces
      //The expression (cefa_)[i][j] < nFaces is always satisfied
      llabelVectorVector cefa_;
      //Face-cells
      //Every face has two faces;
      //For boundary face, the second cell index is -1;
      //If index = (face_)[i][0:1] > nCells_, the second cell is ghost.
      cellVector        cells_;
      faceVector        faces_;

    private:
      primitiveMesh(const primitiveMesh&){};
      //virtual void virtualfun() const = 0;
    protected:
      primitiveMesh(MPI_Comm mpicomm = MPI_COMM_WORLD)
          :distributedObjBase(mpicomm),
          nPoints_(0), nFaces_(0), nCells_(0), nGhostCells_(0)
      {};

    public:
      ~primitiveMesh(){};
    public:
      void resetMesh()
      {
        nCells_ = 0;
        nPoints_= 0;
        nFaces_ = 0;

        nGhostCells_ = 0;
        cellGlobalIdx_.clear();
        ghostCellMPIRank_.clear();
        sendPattern_.clear();
        receivePattern_.clear();

        points_.clear();
        cece_.clear();
        cefa_.clear();
        faces_.clear();
        cells_.clear();
      }
      const faceVector& faces() const{
        return faces_;
      }
      const cellVector& cells() const{
        return cells_;
      }

      const llabelVectorVector& cece()const{
        return cece_;
      }
      const llabelVectorVector& cefa()const{
        return cefa_;
      }
      llabel nCells() const{
        return nCells_;
      }
      llabel nFaces() const{
        return nFaces_;
      }

      llabel nGhostCells() const{
        return nGhostCells_;
      }
      const glabelVector& cellGlobalIdx() const{
        return cellGlobalIdx_;
      }
      const llabelVectorVector& sendPattern() const{
        return sendPattern_;
      }
      const glabelVectorVector& receivePattern() const{
        return receivePattern_;
      }
  };
}
#endif
