#ifndef FPS_MESH_MESH_H
#define FPS_MESH_MESH_H
#include <map>
#include <dataType.h>
#include <real3.h>
#include <primitiveMesh.h>
namespace fps
{
  class mesh
    : public primitiveMesh
  {
    protected:
      //Periodic peroperty of the mesh
      struct _periodicValue{
        //Periodic type:
        //1, Only with translation;
        //2, Rotate + (Translate)
        int type;
        real3 rotateCenter;
        real3 rotateAngle;
        real3 translation;
        _periodicValue(
            int t = 0,
            const real3& rc = real3(0,0,0),
            const real3& ra = real3(0,0,0),
            const real3& tr = real3(0,0,0)
            )
          :type(t),
          rotateCenter(rc),
          rotateAngle(ra),
          translation(tr)
        {
        }
      };
      //fieldValuePtr_ used for detailed solver;
      struct boundaryProperty{
          bool isPeriodic;
          _periodicValue periodicValue;
          void* fieldValuePtr_;
          boundaryProperty(
              bool t1 = false,
              _periodicValue t2 = _periodicValue(),
              void* t3 = nullptr
              )
            :isPeriodic(t1),
            periodicValue(t2),
            fieldValuePtr_(t3)
        {
        }
      };
      std::map<int, boundaryProperty> boundaryPropertyList_;
      std::map<int, std::vector<llabel> > boundaryFaceIdxList_;
    public:
      mesh(MPI_Comm mpicomm = MPI_COMM_WORLD):primitiveMesh(mpicomm){}
      //This list only stores the boundary indices (does not include periodic boundaries)
      const std::map<int, std::vector<llabel> >&
        boundaryFaceIdxList()const{
        return boundaryFaceIdxList_;
      }

      inline void resetMesh(){
        boundaryPropertyList_.clear();
        boundaryFaceIdxList_.clear();
        primitiveMesh::resetMesh();
      }
      ~mesh(){
      }

      std::map<int, boundaryProperty>& boundaryPropertyList(){
        return boundaryPropertyList_;
      }
      const std::map<int, std::vector<llabel> >& boundaryFaceIdxList(){
        return boundaryFaceIdxList_;
      };
  };
}
#endif
