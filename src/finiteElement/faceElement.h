#ifndef FPS_FINITEELEMENT_FACEELEMENT_H
#define FPS_FINITEELEMENT_FACEELEMENT_H
#include <real3.h>
#include <finiteElementBase.h>
namespace fps{
  //Supply class for <class finiteElementBase>
  //Used for face integration
  class faceElement
  {
    public:
      int order;
      int nFaceIntPoints;
      realVector WxJ;
      //intPhyCoord Exact values for the left cell
      //For the right cell, the actual positon values = rotataionMatrix * intPhyCoord + transitionVector
      real3Vector intPhyCoord;
      real3Vector intNormVector;
      //The index of the integration points
      //(first for the left cell and second for the right)
      std::vector<int> leftFaceIntIdx;
      std::vector<int> rightFaceIntIdx;
    public:
      faceElement(int oo, const face&f,
          const finiteElementBase& lc,
          const finiteElementBase& rc)
        :
          order(oo)
      {
        lc.determineFaceIntParameters(f, order, *this);
      }
      //Dummy Constructor
      faceElement(){}
      //Initializer
      void initiate(int oo, const face&f,
          const finiteElementBase& lc
          )
      {
        order = oo;
        lc.determineFaceIntParameters(f, order, *this);
      }
  };
}
#endif
