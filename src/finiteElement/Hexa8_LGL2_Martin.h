#ifndef FPS_FINITEELEMENT_HEXA8_LGL2_MARTIN_H
#define FPS_FINITEELEMENT_HEXA8_LGL2_MARTIN_H
#include <finiteElementBase.h>
#include <faceElement.h>
namespace fps
{
  class Hexa8_LGL2_Martin:
    public finiteElementBase
  {
    public:
      struct refFrame
      {
        static const unsigned int nSolutionPoints;
        static const unsigned int nVertices;
        static const unsigned int type;
        static const unsigned int order;
        static const real3Vector solutionCoord;

        static const unsigned int nVolIntPoints;
        static const real3Vector volIntRefCoord;
        static const real3Vector faceIntRefCoord;
        static const realVector volIntWeight;
        static const realVector faceIntWeight;

        //Phi is the base for the solution point
        static const realVectorVector phiAtVolIntPoint;
        static const realVectorVector phiAtFaceIntPoint;
        static const real3VectorVector dphiAtVolIntPoint;
        static const real3VectorVector dphiAtFaceIntPoint;

        //psi is the base for the vertices
        static const realVectorVector psiAtVolIntPoint;
        static const realVectorVector psiAtFaceIntPoint;
        static const real3VectorVector dpsiAtVolIntPoint;
        static const real3VectorVector dpsiAtFaceIntPoint;
        //determinFaceIntIdx: supposted that both sides are Hexa8_LGL2 cells
        //and the Gaussian points are both supposed to be 2
        static void determineFaceIntIdx(
            //Input parameter
            const face& f,
            //Output parameter
            std::vector<int>& lFaceIntIdx,std::vector<int>& rFaceIntIdx);
      };

    private:
      int determineFaceIntIdx(
          const face& f, int nGaussianPoints,
          std::vector<int>& lFaceIntIdx,
          std::vector<int>& rFaceIntIdx
          )const;

    public:
      int determineFaceIntParameters(
          const face& f,
          int order,
          faceElement& faceElem
          )const;

      Hexa8_LGL2_Martin(const real3Vector& vertices);
      Hexa8_LGL2_Martin():finiteElementBase(1,1){}

      void initiate(const real3Vector& vertices);
      static void obtainCoord(const real3Vector& nodes, const llabelVector& idx,
        const real3Vector& refCoord, real3Vector& phyCoord);

      friend size_t serializedSize(const Hexa8_LGL2_Martin&);
      friend void serialize(const Hexa8_LGL2_Martin&, char*);
      friend void unserialize(Hexa8_LGL2_Martin&, char*);
  };
}
#endif
