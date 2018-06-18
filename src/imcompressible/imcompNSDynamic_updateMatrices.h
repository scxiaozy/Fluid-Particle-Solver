#ifndef FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_UPDATEMATRICES_H
#define FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_UPDATEMATRICES_H
#include <lapacke.h>
#include <imcompNSDynamic.h>
namespace fps
{
  namespace
  {
    inline void obtainInvJacob(
     const real3& dxi,
     const real3& deta,
     const real3& dzeta,
     math::matrix& m
     )
    {
      real jacvol = dxi*crossMult(deta,dzeta);
      assert(m.size().first == 3 && m.size().second == 3 && jacvol >= 0);
      m[0][0] = (deta [1]*dzeta[2] - dzeta[1]*deta [2])/jacvol;
      m[0][1] = (dzeta[0]*deta [2] - deta [0]*dzeta[2])/jacvol;
      m[0][2] = (deta [0]*dzeta[1] - dzeta[0]*deta [1])/jacvol;
      m[1][0] = (dzeta[1]*dxi  [2] - dxi  [1]*dzeta[2])/jacvol;
      m[1][1] = (dxi  [0]*dzeta[2] - dzeta[0]*dxi  [2])/jacvol;
      m[1][2] = (dzeta[0]*dxi  [1] - dxi  [0]*dzeta[1])/jacvol;
      m[2][0] = (dxi  [1]*deta [2] - deta [1]*dxi  [2])/jacvol;
      m[2][1] = (deta [0]*dxi  [2] - dxi  [0]*deta [2])/jacvol;
      m[2][2] = (dxi  [0]*deta [1] - deta [0]*dxi  [1])/jacvol;
    }

  }
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
    poissionMatrixFace.resize(nFaces);

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
        poissionMatrixCell[nc].init(0.0, sps, sps);
        const T& volElem = mesh_.volElems()[nc];
        //If this is a new cell, update the matrix
        //The poission matrix is symetric
        //Volume integration(\nabla phi, \nabla p)
        for (int ni = 0; ni < vips; ++ni){
          DphiDx.init(0.0);
          //For each gaussian point, calculate the Dphi*invJ
          for (int ii = 0; ii < sps; ++ii){
            for (int jj = 0; jj < 3; ++jj){
              for (int kk = 0; kk < 3; ++kk){
                DphiDx[ii][jj] += T::refFrame::dphiAtVolIntPoint[ii][ni][kk]*invVolJacob[nc][3*ni+kk][jj];
              }//End for kk
            }//End for jj
          }//End for sps
          cellMatrixTemp.init(0.0);
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

    const int oo =  T::refFrame::order;
    //Face integration for poissionMatrixCell and poissionMatrixFace
    const std::vector<faceElement>& ffem = mesh_.faceElems1();
    for (int nf = 0; nf < nFaces; ++nf){
      llabel ncl = mesh_.faces()[nf].lcellIdx.first;
      llabel ncr = mesh_.faces()[nf].lcellIdx.second;
      const faceElement& fe = ffem[nf];
      //If this is a new face
      llabel oldFaceIdx = mesh_.faces()[nf].oldFaceIdx;
      if(oldFaceIdx == -1){
        //If is a boundary face
        if(ncr == -1) continue;
        //Otherwise
        poissionMatrixFace[nf].init(0.0,sps,sps);
        real voll = 0;
        real volr = 0;
        real S = 0;
        for (int ni = 0; ni < vips; ++ni){
          voll += mesh_.volElems()[ncl].volIntWxJ[ni]; 
          volr += mesh_.volElems()[ncr].volIntWxJ[ni]; 
        }
        for (int np = 0; np < fe.nFaceIntPoints; ++np){
          S += fe.WxJ[np];
        }
        real tau_ip =  S*0.5/std::min(voll,volr)*(oo+1)*(oo+1);
        //First calculate the possionMatrixCell
        std::cout << this->mpiRank() << " in " << __FILE__ << " : " << __LINE__ << std::endl;
        for (int np = 0; np < fe.nFaceIntPoints; ++np){
          int idxl = fe.leftFaceIntIdx[np];
          int idxr = fe.rightFaceIntIdx[np];
          //Calculate the jacobian matrix of this point;
          math::matrix invJacobl(3,3);
          math::matrix invJacobr(3,3);
          real3 dxil  (0,0,0);
          real3 detal (0,0,0);
          real3 dzetal(0,0,0);
          real3 dxir  (0,0,0);
          real3 detar (0,0,0);
          real3 dzetar(0,0,0);
          for (unsigned int jj = 0; jj < T::refFrame::nVertices; ++jj){
            dxil      += T::refFrame::dpsiAtFaceIntPoint[jj][idxl][0]*mesh_.volElems()[ncl].vertices[jj];
            detal     += T::refFrame::dpsiAtFaceIntPoint[jj][idxl][1]*mesh_.volElems()[ncl].vertices[jj];
            dzetal    += T::refFrame::dpsiAtFaceIntPoint[jj][idxl][2]*mesh_.volElems()[ncl].vertices[jj];
            dxir      += T::refFrame::dpsiAtFaceIntPoint[jj][idxr][0]*mesh_.volElems()[ncr].vertices[jj];
            detar     += T::refFrame::dpsiAtFaceIntPoint[jj][idxr][1]*mesh_.volElems()[ncr].vertices[jj];
            dzetar    += T::refFrame::dpsiAtFaceIntPoint[jj][idxr][2]*mesh_.volElems()[ncr].vertices[jj];
          }
          obtainInvJacob(dxil, detal, dzetal, invJacobl);
          obtainInvJacob(dxir, detar, dzetar, invJacobr);
          const real3& norm = fe.intNormVector[np];
          real3 templ;
          real3 tempr;
          for (int ii = 0; ii < 3; ++ii){
            templ[ii] = 
             invJacobl[ii][0] * norm[0]
            +invJacobl[ii][1] * norm[1]
            +invJacobl[ii][2] * norm[2];

            tempr[ii] = 
            -invJacobr[ii][0] * norm[0]
            -invJacobr[ii][1] * norm[1]
            -invJacobr[ii][2] * norm[2];
          }
          real dphidnl[sps];
          real dphidnr[sps];
          for (int ni = 0; ni < sps; ++ni){
            dphidnl[ni] = T::refFrame::dphiAtFaceIntPoint[ni][idxl]*templ;
            dphidnr[ni] = T::refFrame::dphiAtFaceIntPoint[ni][idxr]*tempr;
          }
          real phi_l[sps];
          real phi_r[sps];
          for (int ni = 0; ni < sps; ++ni) phi_l[ni] = T::refFrame::phiAtFaceIntPoint[ni][idxl]*fe.WxJ[np];
          for (int ni = 0; ni < sps; ++ni) phi_r[ni] = T::refFrame::phiAtFaceIntPoint[ni][idxr]*fe.WxJ[np];

          if(idxl < nCells)
          for (int ni = 0; ni < sps; ++ni){
            for (int nj = 0; nj < sps; ++nj){
              poissionMatrixCell[idxl][ni][nj] +=
              (
               -0.5*dphidnl[ni]*phi_l[nj]
               -0.5*dphidnl[nj]*phi_l[ni]
               +tau_ip*T::refFrame::phiAtFaceIntPoint[ni][idxl]*phi_l[nj]
              );
            }
          }

          if(idxr < nCells)
          for (int ni = 0; ni < sps; ++ni){
            for (int nj = 0; nj < sps; ++nj){
              poissionMatrixCell[idxr][ni][nj] +=
              (
               -0.5*dphidnr[ni]*phi_r[nj]
               -0.5*dphidnr[nj]*phi_r[ni]
               +tau_ip*T::refFrame::phiAtFaceIntPoint[ni][idxr]*phi_r[nj]
              );
            }
          }

          //Then calculate the face matrix;
          for (int ni = 0; ni < sps; ++ni){
            for (int nj = 0; nj < sps; ++nj){
              //Note: dphidnr use -normvector;
              poissionMatrixFace[nf][ni][nj] +=
              0.5*(
                dphidnl[ni]*phi_r[nj] + phi_l[ni]*dphidnr[nj]
              )
              -tau_ip*(T::refFrame::phiAtFaceIntPoint[ni][idxl]*phi_r[nj]);
            }
          }
        }
      }
      else{
        //Swap the memory of this face and the new face;
        poissionMatrixFace[nf].swap(poissionMatrixFaceTemp[oldFaceIdx]);
      }
    }
    delete[] ipiv;
  }
}
#endif