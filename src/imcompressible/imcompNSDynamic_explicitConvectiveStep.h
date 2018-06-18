#ifndef FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_EXPLICITCONVECTIVESTEP_H
#define FPS_IMCOMPRESSIBLE_IMCOMPNSDYNAMIC_EXPLICITCONVECTIVESTEP_H
#include <imcompNSDynamic.h>
namespace fps
{
  namespace{
    void imcompInvisicidFlux(
      //Input parameters
      const real3& vell, const real3& velr, const real3& norm,
      //Output parameters
       real3& flux)
    {
      real vnl = vell*norm;
      real vnr = velr*norm;
      real3 fluxl = vell*vnl;
      real3 fluxr = velr*vnr;
      real3 diff = vell-velr;
      //The eigenvalues are 2*vn, vn, vn
      //The lax flux is used.
      flux = (fluxl+fluxr)/2.0 + std::max(vnl,vnr)*diff;
    }
  }

  template <class T>
  void imcompNSDynamic<T>::explicitConvectiveStep()
  {
    llabel nCells = mesh_.nCells();
    llabel nGhostCells = mesh_.nGhostCells();
    llabel nFaces = mesh_.nFaces();
    int step = std::min(start, timeAccuracy);

    //Backup the solution
    for (int ss = step - 1; ss > 0; --ss){
      for (llabel nc = 0; nc < nCells; ++nc){
        up[ss][nc] = up[ss-1][nc];
        rhsConvectiveP[ss][nc] = rhsConvectiveP[ss-1][nc];
      }
    }
    for (llabel nc = 0; nc < nCells; ++nc){
      up[0][nc] = u[nc];
      rhsConvectiveP[0][nc] = rhsConvective[nc];
    }


    //Firslty calculate the volume flux and the source term
    for (llabel nc = 0; nc < nCells; ++nc){
      rhsConvective[nc].init(0.0);
      rhsSource[nc].init(0.0);
      //Integration over volume
      const size_t vips = T::refFrame::nVolIntPoints;
      //Firstly calculate u,v,w in each gaussian integration point
      real vel[vips][3];
      memset(*vel, 0, sizeof(real)*vips*3);
      for (int jj = 0; jj < T::refFrame::nSolutionPoints; ++jj){
        for (int ii = 0; ii < vips; ++ii){
          vel[ii][0] += u[nc][jj][0]*T::refFrame::phiAtVolIntPoint[jj][ii];
          vel[ii][1] += u[nc][jj][1]*T::refFrame::phiAtVolIntPoint[jj][ii];
          vel[ii][2] += u[nc][jj][2]*T::refFrame::phiAtVolIntPoint[jj][ii];
        }
      }

      //Calculate the flux in the volintegration points;
      real flux[vips][3][3];
      for (int ii = 0; ii < vips; ++ii){
        flux[ii][0][0] = vel[ii][0]*vel[ii][0];
        flux[ii][0][1] = vel[ii][0]*vel[ii][1];
        flux[ii][0][2] = vel[ii][0]*vel[ii][2];
        flux[ii][1][1] = vel[ii][1]*vel[ii][1];
        flux[ii][1][2] = vel[ii][1]*vel[ii][2];
        flux[ii][2][2] = vel[ii][2]*vel[ii][2];
        flux[ii][1][0] = flux[ii][0][1];
        flux[ii][2][0] = flux[ii][0][2];
        flux[ii][2][1] = flux[ii][1][2];
      }

      for (int ii = 0; ii < vips; ++ii){
        real temp[3][3];
        for (int ni = 0; ni < 3; ++ni){
          for (int nj = 0; nj < 3; ++nj){
            temp[ni][nj] = 0;
            for (int nk = 0; nk < 3; ++nk){
              temp[ni][nj] += invVolJacob[nc][3*ii+ni][nk]*flux[ii][nj][nk];
            }
          }
        }
        for (int ni = 0; ni < 3; ++ni){
          for (int nj = 0; nj < 3; ++nj){
            temp[ni][nj] *= mesh_.volElems()[nc].volIntWxJ[ii];
          }
        }

        for (int ni = 0; ni < T::refFrame::nSolutionPoints; ++ni){
          for (int nj = 0; nj < 3; ++nj){
            for (int nk = 0; nk < 3; ++nk){
              rhsConvective[nc][ni][nj] += T::refFrame::dphiAtVolIntPoint[ni][ii][nk]*temp[nk][nj];
            }
          }
        }
      }
    }
    //Secondly calculate the face flux
    const std::vector<faceElement> &ffem = mesh_.faceElems2();
    for (llabel nf = 0; nf < nFaces; ++nf){
      llabel ncl = mesh_.faces()[nf].lcellIdx.first;
      llabel ncr = mesh_.faces()[nf].lcellIdx.second;
      const faceElement& fe = ffem[nf];
      if(ncr == -1){
        //Boundary treatment;
      }else{
        for (int np = 0; np < fe.nFaceIntPoints; ++np){
          //Obtain the u v w for left and right cells at this point;
          int idxl = fe.leftFaceIntIdx[np];
          int idxr = fe.rightFaceIntIdx[np];
          real3 vell(0,0,0);
          real3 velr(0,0,0);
          //The face interpolation should be improved later!!!!!!!!!!!!!!!!!!
          //---------------------TO BE IMPROVED LATER------------------------
          for (int ii = 0; ii < T::refFrame::nSolutionPoints; ++ii){
            vell[0] += T::refFrame::phiAtFaceIntPoint[ii][idxl]*u[ncl][ii][0];
            vell[1] += T::refFrame::phiAtFaceIntPoint[ii][idxl]*u[ncl][ii][1];
            vell[2] += T::refFrame::phiAtFaceIntPoint[ii][idxl]*u[ncl][ii][2];

            velr[0] += T::refFrame::phiAtFaceIntPoint[ii][idxr]*u[ncr][ii][0];
            velr[1] += T::refFrame::phiAtFaceIntPoint[ii][idxr]*u[ncr][ii][1];
            velr[2] += T::refFrame::phiAtFaceIntPoint[ii][idxr]*u[ncr][ii][2];
          }
          //---------------------END TO BE IMPROVED LATER---------------------
          const real3& norm = fe.intNormVector[np];
          //Currently only the translation periodic boundary conditions are 
          //considered, thus, normal vector is the same for left and right cell.
          //And the velr remain unchanged for the priodic boundary conditions.
          //For the rotational boundary conditions, both should be translated
          //using the rotational matrix.
          real3 flux;
          imcompInvisicidFlux(vell, velr, norm, flux);
          flux*=fe.WxJ[np];
          for (int ii = 0; ii < T::refFrame::nSolutionPoints; ++ii){
              //For the left cell
              rhsConvective[ncl][ii][0] -= T::refFrame::phiAtFaceIntPoint[ii][idxl]*flux[0];
              rhsConvective[ncl][ii][1] -= T::refFrame::phiAtFaceIntPoint[ii][idxl]*flux[1];
              rhsConvective[ncl][ii][2] -= T::refFrame::phiAtFaceIntPoint[ii][idxl]*flux[2];
              //For the right cell
              rhsConvective[ncr][ii][0] += T::refFrame::phiAtFaceIntPoint[ii][idxr]*flux[0];
              rhsConvective[ncr][ii][1] += T::refFrame::phiAtFaceIntPoint[ii][idxr]*flux[1];
              rhsConvective[ncr][ii][2] += T::refFrame::phiAtFaceIntPoint[ii][idxr]*flux[2];
          }
        }
      }
    }

    //Thirdly, update the velocity field
    for (llabel nc = 0; nc < nCells; ++nc){
      const int sps = T::refFrame::nSolutionPoints;
      real3 rhs[sps];
      //Firslty calculate rhs
      for (int ii = 0; ii < sps; ++ii){
        //u[nc][ii][0] = cofAlpha[step-1][0]*up[0][nc][ii][0];
        //u[nc][ii][1] = cofAlpha[step-1][0]*up[0][nc][ii][1]; 
        //u[nc][ii][2] = cofAlpha[step-1][0]*up[0][nc][ii][2];
        rhs[ii][0] = cofBeta[step - 1][0] * rhsConvective[nc][ii][0];
        rhs[ii][1] = cofBeta[step - 1][0] * rhsConvective[nc][ii][1];
        rhs[ii][2] = cofBeta[step - 1][0] * rhsConvective[nc][ii][2];
        for (int ss = 1; ss < step; ++ss){
           rhs[ii][0] += cofBeta[step-1][ss]*rhsConvectiveP[ss][nc][ii][0]; 
           rhs[ii][1] += cofBeta[step-1][ss]*rhsConvectiveP[ss][nc][ii][1]; 
           rhs[ii][2] += cofBeta[step-1][ss]*rhsConvectiveP[ss][nc][ii][2]; 
        }
        rhs[ii] *= physicalDt;
        //for (int ss = 1; ss < step; ++ss){
        //   u[nc][ii][0] += cofAlpha[step-1][ss]*up[ss][nc][ii][0]; 
        //   u[nc][ii][1] += cofAlpha[step-1][ss]*up[ss][nc][ii][1]; 
        //   u[nc][ii][2] += cofAlpha[step-1][ss]*up[ss][nc][ii][2]; 
        //}
      }
      //Then multiply the rhs by invJacobMatrix
      //Automatically assigned 0 in matrix constructor;
      math::matrix temp(sps, 3);
      for (int ni = 0; ni < sps; ++ni){
        for (int nj = 0; nj < 3; ++nj){
          for (int nk = 0; nk < sps; ++nk){
            temp[ni][nj] += invVolJacob[nc][ni][nk]*rhs[nk][nj];
          }
        }
      }
      for (int ii = 0; ii < sps; ++ii){
        u[nc][ii][0] = cofAlpha[step-1][0]*up[0][nc][ii][0];
        u[nc][ii][1] = cofAlpha[step-1][0]*up[0][nc][ii][1]; 
        u[nc][ii][2] = cofAlpha[step-1][0]*up[0][nc][ii][2];
        for (int ss = 1; ss < step; ++ss){
           u[nc][ii][0] += cofAlpha[step-1][ss]*up[ss][nc][ii][0]; 
           u[nc][ii][1] += cofAlpha[step-1][ss]*up[ss][nc][ii][1]; 
           u[nc][ii][2] += cofAlpha[step-1][ss]*up[ss][nc][ii][2]; 
        }
        u[nc][ii][0] += temp[ii][0];
        u[nc][ii][1] += temp[ii][1];
        u[nc][ii][2] += temp[ii][2];
        u[nc][ii][0] += rhsSource[nc][ii][0]*physicalDt;
        u[nc][ii][1] += rhsSource[nc][ii][1]*physicalDt;
        u[nc][ii][2] += rhsSource[nc][ii][2]*physicalDt;
      }
    }
  }
}
#endif