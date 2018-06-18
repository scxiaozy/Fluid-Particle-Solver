#ifndef FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_UPDATEMATRICES_H
#define FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_UPDATEMATRICES_H
#include <lapacke.h>
#include <imcompNSDynamic.h>
namespace fps
{
  template<class T>
  void imcompNSDynamic<T>::updateMatrices()
  {
    llabel nCells = mesh_.nCells();
    llabel nFaces = mesh_.nFaces();
    //volume integration points
    const unsigned int vips = T::refFrame::nVolIntPoints;
    //solution points
    const unsigned int sps = T::refFrame::nSolutionPoints;
    int* ipiv = new int[sps+1];
    std::vector<math::matrix> matrixTemp;
    std::vector<math::matrix> invVolJacobTemp;
    std::vector<math::matrix> poissionMatrixCellTemp;
    std::vector<math::matrix> poissionMatrixFaceTemp;

    matrixTemp.swap(invMassMatrix);
    invVolJacobTemp.swap(invVolJacob);
    poissionMatrixCell.swap(poissionMatrixCellTemp);
    poissionMatrixFace.swap(poissionMatrixFaceTemp);
    
    invMassMatrix.resize(nCells);
    invVolJacob.resize(nCells);
    poissionMatrixCell.resize(nCells);
    poissionMatrixFace.resize(nCells);

    //First Section, update the mass matrix used for the explicit convective step
    for (int nc = 0; nc < nCells; ++nc){
      if(mesh_.mapVector()[nc].fineCoarsenTag){
        double massMatrix[sps][sps];
        memset(massMatrix,0,sizeof(double)*sps*sps);
        for (int ii = 0; ii < sps; ++ii){
          for (int jj = 0; jj < sps; ++jj){
            for (int kk = 0; kk < vips; ++kk){
              massMatrix[ii][jj] +=
               T::refFrame::phiAtVolIntPoint[ii][kk]
                *T::refFrame::phiAtVolIntPoint[jj][kk]
                *mesh_.volElems()[nc].volIntWxJ[kk];
            }
          }
        }
        lapack_int ret;
        ret =  LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
                          sps,
                          sps,
                          *massMatrix,
                          sps,
                          ipiv);

        if (ret !=0){
          std::cerr << "Error happend in inversing matrix" << std::endl;
          exit(-1);
        }

        ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
                           sps,
                           *massMatrix,
                           sps,
                           ipiv);
        if (ret !=0){
          std::cerr << "Error happend in inversing matrix" << std::endl;
          exit(-1);
        }
        invMassMatrix[nc].assign(*massMatrix, sps, sps);

        invVolJacob[nc].resize(vips*3,3);
        //Update invVolJacobMatrix for each volIntegration point
        for (int ii = 0; ii < vips; ++ii){
          real3 dxi  (0,0,0);
          real3 deta (0,0,0);
          real3 dzeta(0,0,0);
          for (unsigned int jj = 0; jj < T::refFrame::nVertices; ++jj){
            dxi      += T::refFrame::dpsiAtVolIntPoint[jj][ii][0]*mesh_.volElems()[nc].vertices[jj];
            deta     += T::refFrame::dpsiAtVolIntPoint[jj][ii][1]*mesh_.volElems()[nc].vertices[jj];
            dzeta    += T::refFrame::dpsiAtVolIntPoint[jj][ii][2]*mesh_.volElems()[nc].vertices[jj];
          }
          real jacvol = dxi*crossMult(deta,dzeta);
          assert(jacvol >= 0);
          invVolJacob[nc][3*ii+0][0] = (deta [1]*dzeta[2] - dzeta[1]*deta [2])/jacvol;
          invVolJacob[nc][3*ii+0][1] = (dzeta[0]*deta [2] - deta [0]*dzeta[2])/jacvol;
          invVolJacob[nc][3*ii+0][2] = (deta [0]*dzeta[1] - dzeta[0]*deta [1])/jacvol;
          invVolJacob[nc][3*ii+1][0] = (dzeta[1]*dxi  [2] - dxi  [1]*dzeta[2])/jacvol;
          invVolJacob[nc][3*ii+1][1] = (dxi  [0]*dzeta[2] - dzeta[0]*dxi  [2])/jacvol;
          invVolJacob[nc][3*ii+1][2] = (dzeta[0]*dxi  [1] - dxi  [0]*dzeta[1])/jacvol;
          invVolJacob[nc][3*ii+2][0] = (dxi  [1]*deta [2] - deta [1]*dxi  [2])/jacvol;
          invVolJacob[nc][3*ii+2][1] = (deta [0]*dxi  [2] - dxi  [0]*deta [2])/jacvol;
          invVolJacob[nc][3*ii+2][2] = (dxi  [0]*deta [1] - deta [0]*dxi  [1])/jacvol;
        }
      }
      else{
        llabel oldIdx = mesh_.mapVector()[nc].fcTag.equal.oldLocalIdx;
        invMassMatrix[nc].swap(matrixTemp[oldIdx]);
        invVolJacob[nc].swap(invVolJacobTemp[oldIdx]);
      }
    }
    //Second Section, update the matrix for the poission equation, including the
    //poissionMatrixDiag and poissionMatrixNeighb
    math::matrix DphiDx(sps, 3);
    math::matrix cellMatrixTemp(sps,sps);

    for (int nc = 0; nc < nCells; ++nc){
      if(mesh_.mapVector()[nc].fineCoarsenTag){
        poissionMatrixCell[nc].assign(0.0, sps, sps);
        const T& volElem = mesh_.volElems()[nc];
        //If this is a new cell, update the matrix
        //The poission matrix is symetric
        //Volume integration(\nabla phi, \nabla p)
        for (int ni = 0; ni < vips; ++ni){
          DphiDx.assign(0.0);
          //For each gaussian point, calculate the Dphi*invJ
          for (int ii = 0; ii < sps; ++ii){
            for (int jj = 0; jj < 3; ++jj){
              for (int kk = 0; kk < 3; ++kk){
                DphiDx[ii][jj] += T::refFrame::dphiAtVolIntPoint[ii][ni][kk]*invVolJacob[nc][3*ni+kk][jj];
              }//End for kk
            }//End for jj
          }//End for sps
          cellMatrixTemp.assign(0.0);
          //Then calculate DphiDx*DphiDx^T;
          for (int ii = 0; ii < sps; ++ii){
            for (int jj = 0; jj < sps; ++jj){
              for (int kk = 0; kk < 3; ++kk){
                cellMatrixTemp[ii][jj] += DphiDx[ii][kk]*DphiDx[jj][kk];
              }//End for kk
              cellMatrixTemp[ii][jj] *= volElem.volIntWxJ[ni];
              poissionMatrixCell[nc][ii][jj] += cellMatrixTemp[ii][jj];
            }//End for jj
          }//End for ii
        }//End for vips
      }else{
        //Else, just swap the memory of the old matrix
        llabel oldIdx = mesh_.mapVector()[nc].fcTag.equal.oldLocalIdx;
        poissionMatrixCell[nc].swap(poissionMatrixCellTemp[oldIdx]);
      }
    }

    //Face integration for poissionMatrixCell and poissionMatrixFace
    for (int nf = 0; nf < nFaces; ++nf){
    
    }
    delete[] ipiv;
  }
}
#endif