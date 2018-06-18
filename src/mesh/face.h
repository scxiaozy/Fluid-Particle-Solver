#ifndef FPS_MESH_FACE_H
#define FPS_MESH_FACE_H
#include <cstdint>

#include <dataType.h>
#include <real3.h>
namespace fps
{
  class face
  {
    public:
      llabelPair lcellIdx;
      glabelPair gcellIdx;
      //The side of the left face
      int8_t firstSide;
      //The side of the right face
      int8_t secondSide;
      //The child id of the left face
      //Only the left face could be a hanging face
      //Used for dynamic mesh.
      //If ChildIdx == -1, this is not a hanging face
      //If ChildIdx == 0,1,2,3 or others. (Defined by detailed cell type)
      int8_t firstFaceChildIdx;
      int8_t secondFaceChildIdx;
      int8_t orientation;
      int8_t firstCellType;
      int8_t secondCellType;
      int8_t isHanging;
      //For the outside of a boundary, patchIdx.second = -1;
      intPair patchIdx;
      face (const llabelPair& lcell = llabelPair(-1,-1),
          const glabelPair& gcell = glabelPair(-1,-1),
          const intPair& pIdx= intPair(-1,-1))
        : lcellIdx(lcell), gcellIdx(gcell),
        firstSide(-1), secondSide(-1), firstFaceChildIdx(-1), secondFaceChildIdx(-1),
        orientation(-1), firstCellType(-1), secondCellType(-1),isHanging(0), patchIdx(pIdx)
      {
      };
  };
  using faceVector = std::vector<face>;
}
#endif
