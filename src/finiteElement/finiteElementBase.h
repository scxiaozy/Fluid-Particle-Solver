#ifndef FPS_FINITEELEMENT_FINITE_ELEMENT_BASE_H
#define FPS_FINITEELEMENT_FINITE_ELEMENT_BASE_H
#include <face.h>
#include <real3.h>
namespace fps
{
  //Forward declaration
  class faceElement;

  //Only the properties in the physical frame are defined
  class finiteElementBase
  {
    public:
      int type;
      int order;
      std::vector<real3> vertices;

      //The coordinates of the degrees of freedom in physical frame
      real3Vector solutionPhyCoord;

      //The weightXJacobian for each volume integration points
      //Vertices for the finite element
      realVector volIntWxJ;
      real3Vector volIntPhyCoord;
    public:
      //Tyep = -1 means the type of this FE is undefined.
      finiteElementBase(int t_ = -1, int o_ = -1)
        :type(t_), order (o_){}


    public:
      virtual int determineFaceIntParameters(
          const face& f,
          int order,
          faceElement&) const{return 0;}

      friend size_t serializedSize(const finiteElementBase&);
      friend void serialize(const finiteElementBase&, char*);
      friend void unserialize(finiteElementBase&, char*);
  };
}
#endif
