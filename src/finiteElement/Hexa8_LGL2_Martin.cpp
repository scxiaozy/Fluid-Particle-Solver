#include <cassert>
#include <Hexa8_LGL2_Martin.h>
#include <face.h>
#include <real3.h>
namespace fps{
  void Hexa8_LGL2_Martin::refFrame::determineFaceIntIdx
    (
     const face& f,
     std::vector<int>& lFaceIntIdx,
     std::vector<int>& rFaceIntIdx
     )
    {
      //R Q T P SS: assistant matrices of transformation
      //The left face's x coord corresponds with
      //right face's
      // P[Q[R[f.firstSide][f.secondSide]][f.orientation]][0]
      //and direction is
      // S[Q[R[f.firstSide][f.secondSide]][f.orientation]][0]
      //Same for the y coord.
      //The left face's corner $m$ corresponds with right face's
      //corner
      // T[Q[R[f.firstSide][f.secondSide]][f.orientation]][m]
      // and vice versa.
      const int R[6][6] = {
        {0,1,1,0,0,1},
        {2,0,0,1,1,0},
        {2,0,0,1,1,0},
        {0,2,2,0,0,1},
        {0,2,2,0,0,1},
        {2,0,0,2,2,0},
      };
      const int Q[3][4] = {
        {1,2,5,6},
        {0,3,4,7},
        {0,4,3,7},
      };
      //const int T[8][4] = {
        //{0,1,2,3},
        //{0,2,1,3},
        //{1,0,3,2},
        //{1,3,0,2},
        //{2,0,3,1},
        //{2,3,0,1},
        //{3,1,2,0},
        //{3,2,1,0},
      //};
      //The corresponding coordinates of the face-local x & y
      //1: x coordinate
      //2: y coordinate
      const int P[8][2] = {
        {0,1},
        {1,0},
        {0,1},
        {1,0},
        {1,0},
        {0,1},
        {1,0},
        {0,1}
      };
      //The corresponding coordinate direction
      //of the face-local x & y
      //1: x coordinate
      //2: y coordinate
      const int S[8][2] = {
        {+1,+1},
        {+1,+1},
        {-1,+1},
        {+1,-1},
        {-1,+1},
        {+1,-1},
        {-1,-1},
        {-1,-1}
      };
      //Note: the hanging face always locates at right
      int l_strip = 2; //The length of Gaussian point strip in 1-D of order 2
      int l_shift = 0; //The shift index for the hanging face;
      int l_shift_i = 0;
      int l_shift_j = 0;
      if (f.isHanging){
        l_strip = 4;
        //The start index of the sub_face Gaussian points is 6\times 2\times 2
        l_shift = 6*2*2;
        //Since the hanging face always locates at right, the left face may be a sub_face.
        //f.firstFaceChildIdx denotes the sub_face number.
        //If the face is not a sub_face, this number is 0
        l_shift_i = f.firstFaceChildIdx%2;
        l_shift_j = f.firstFaceChildIdx/2;
      }
      lFaceIntIdx.clear();
      rFaceIntIdx.clear();
      for (int jj = 0; jj < 2; ++jj){
        for (int ii = 0; ii < 2; ++ii){
          lFaceIntIdx.push_back(l_shift + f.firstSide *l_strip*l_strip
              + (jj + l_shift_j)*l_strip+(ii+l_shift_i));
        }
      }
      //If the right side is a boundary assign lFaceIntIdx to rFaceIntIdx
      if(f.secondSide == -1){
        rFaceIntIdx = lFaceIntIdx;
        return;
      }

      //Notes, hanging face always locates at right
      //So the r_strip is always 2 for HEXA8_LGL2
      const int r_strip = 2;
      int tdx[2] = {0,0};
      for (int jj = 0; jj < 2; ++jj){
        for (int ii = 0; ii < 2; ++ii){
          int t = Q[R[f.firstSide][f.secondSide]][f.orientation];
          tdx[P[t][0]] = (1-S[t][0])/2*(r_strip-1) + S[t][0]*ii;
          tdx[P[t][1]] = (1-S[t][1])/2*(r_strip-1) + S[t][1]*jj;
          //Then transform the neighboring face's dx dy into indices
          rFaceIntIdx.push_back(f.secondSide*r_strip*r_strip + tdx[1]*r_strip+tdx[0]);
        }
      }
    }



  int Hexa8_LGL2_Martin::determineFaceIntIdx(
      const face& f, int order,
      std::vector<int>& lFaceIntIdx,
      std::vector<int>& rFaceIntIdx
      )const
  {
    //Currently only Hexa8 with LGL2_Martin(type == 1, order == 2) cell is supported!
    assert(f.secondCellType == 1);
    assert(order == 3);
    assert(0 <= f.firstSide && f.firstSide < 6);
    Hexa8_LGL2_Martin::refFrame::determineFaceIntIdx(f, lFaceIntIdx, rFaceIntIdx);
    return 0;
  }

  int Hexa8_LGL2_Martin::determineFaceIntParameters(
      const face& f,
      //Order of Gaussian Legendre points
      //For LGL_2 this number should be 2;
      int order,
      faceElement& faceElem
      )const
  {
    assert(f.firstSide < 6 && f.firstSide >=0);
    assert(f.secondSide < 6 && f.secondSide >=0);
    this->determineFaceIntIdx(
        f,order,
        faceElem.leftFaceIntIdx, faceElem.rightFaceIntIdx
        );
    //The integration precision
    faceElem.order = 3;
    faceElem.nFaceIntPoints = 4;
    //The normal vectors should always point from left to right
    int fi,fj;
    switch(f.firstSide){
      case 0://I left face
        fi = 2;
        fj = 1;
        break;
      case 1://I right face
        fi = 1;
        fj = 2;
        break;
      case 2://J left face
        fi = 0;
        fj = 2;
        break;
      case 3:
        fi = 2;
        fj = 0;
        break;
      case 4:
        fi = 1;
        fj = 0;
        break;
      case 5:
        fi = 0;
        fj = 1;
        break;
    }
    faceElem.intNormVector.clear();
    for (unsigned int ii  = 0; ii < faceElem.leftFaceIntIdx.size(); ++ii){
      real3 location(0,0,0);
      real3 dxi(0,0,0);
      real3 deta(0,0,0);
      int& fidx = faceElem.leftFaceIntIdx[ii];
      for (unsigned int jj = 0; jj < refFrame::nVertices; ++jj){
        location += refFrame:: psiAtFaceIntPoint[jj][fidx]    *this->vertices[jj];
        dxi      += refFrame::dpsiAtFaceIntPoint[jj][fidx][fi]*this->vertices[jj];
        deta     += refFrame::dpsiAtFaceIntPoint[jj][fidx][fj]*this->vertices[jj];
      }
      faceElem.intPhyCoord.push_back(location);
      real3 norm = crossMult(dxi,deta);
      real  dA = norm.mod();
      norm/=dA;
      faceElem.WxJ.push_back(dA*refFrame::faceIntWeight[fidx]);
      faceElem.intNormVector.push_back(norm);
    }
    return 0;
  }

  //Hexa8_LGL2_Martin::
    //Hexa8_LGL2_Martin(void)
    //:finiteElementBase(Hexa8_LGL2_Martin::refFrame::type)
  //{
    //this->vertices = refFrame::solutionCoord;
  //}

  void Hexa8_LGL2_Martin:: initiate(const real3Vector& vertices )
  {
    this->vertices = vertices;
    this->solutionPhyCoord.resize(refFrame::nSolutionPoints);
    this->volIntPhyCoord.resize(refFrame::nVolIntPoints);
    this->volIntWxJ.assign(refFrame::nVolIntPoints, 0);
    ////For Hexa8_LGL2 the solution points overlap with the vertices
    for (unsigned int jj = 0; jj < refFrame::nSolutionPoints; ++jj){
      this->solutionPhyCoord[jj]= this->vertices[jj];
    }

    for (unsigned int jj = 0; jj < refFrame::nVolIntPoints; ++jj){
      real3 dxi (0,0,0), deta(0,0,0), dzeta(0,0,0);
      for (unsigned int ii = 0; ii < refFrame::nVertices; ++ii){
        this->volIntPhyCoord[jj] += refFrame::psiAtVolIntPoint[ii][jj]*this->vertices[ii];
        dxi   += refFrame::dpsiAtVolIntPoint[ii][jj][0]*this->vertices[ii];
        deta  += refFrame::dpsiAtVolIntPoint[ii][jj][1]*this->vertices[ii];
        dzeta += refFrame::dpsiAtVolIntPoint[ii][jj][2]*this->vertices[ii];
      }
      this->volIntWxJ[jj] = dxi*crossMult(deta,dzeta)*refFrame::volIntWeight[jj];
      assert(this->volIntWxJ[jj] >= 0);
    }
  }

  Hexa8_LGL2_Martin::
    Hexa8_LGL2_Martin(const real3Vector& vertices )
    :finiteElementBase(1,1)
  {
    initiate(vertices);
  }

  void Hexa8_LGL2_Martin::obtainCoord(const real3Vector& nodes, const llabelVector& idx,
    const real3Vector& refCoord, real3Vector& phyCoord)
  {
    size_t ss = refCoord.size();
    phyCoord.resize(ss);
    //Init the linear finite elment
    for (unsigned int kk = 0; kk < ss; ++kk){
      real x0 = (1.0-refCoord[kk][0])/2.0;
      real x1 = (1.0+refCoord[kk][0])/2.0;
      real y0 = (1.0-refCoord[kk][1])/2.0;
      real y1 = (1.0+refCoord[kk][1])/2.0;
      real z0 = (1.0-refCoord[kk][2])/2.0;
      real z1 = (1.0+refCoord[kk][2])/2.0;
      phyCoord[kk] =
        z0*(
            y0*
            (
              x0*nodes[idx[0]]
             +x1*nodes[idx[1]]
            )
            +
            y1*
            (
              x0*nodes[idx[2]]
             +x1*nodes[idx[3]]
            )
            )
        +
        z1*(
            y0*
            (
              x0*nodes[idx[4]]
             +x1*nodes[idx[5]]
            )
            +
            y1*
            (
              x0*nodes[idx[6]]
             +x1*nodes[idx[7]]
            )
            )
        ;
    }//End for kk
  }
  size_t serializedSize(const Hexa8_LGL2_Martin& elem)
  {
    return serializedSize(static_cast<finiteElementBase>(elem));
  }
  void serialize(const Hexa8_LGL2_Martin& elem, char* ptr)
  {
    serialize(static_cast<const finiteElementBase&>(elem), ptr);
  }
  void unserialize(Hexa8_LGL2_Martin& elem, char* ptr)
  {
     unserialize(static_cast<finiteElementBase&>(elem), ptr);
  }
}
