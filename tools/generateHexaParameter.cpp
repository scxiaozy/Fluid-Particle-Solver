#include <iostream>
#include <iomanip>
#include <cmath>

double l1(double x){return (1.0-x)/2.0;}
double dl1(double x){return -0.5;}
double l2(double x){return (1.0+x)/2.0;}
double dl2(double x){return 0.5;}
using _fp = double (*)(double x);
_fp  phi[2] = {l1,l2};
_fp  dphi[2] = {dl1,dl2};

void generateHexa8_LGL2_Martin()
{
  double volIntPoint[8][3];
  double volIntPointW[8];
  double faceIntPoint[120][3];
  double faceIntPointWeight[120];
  double phiAtVolIntPoint[8][8];
  double phiAtFaceIntPoint[8][144];
  double dphiAtVolIntPoint[8][8][3];
  double dphiAtFaceIntPoint[8][144][3];
  double lgl[2]  = {-1.0, 1.0};

  double gl2[2]  = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
  //double gl2[2]  = {-1.0, 1.0};
  double gl2_w[2]= {1.0, 1.0};
  double gl2_half[4] = {(gl2[0]-1)/2.0, (gl2[1]-1.0)/2.0, (gl2[0]+1)/2.0, (gl2[1]+1.0)/2.0};
  double gl2_half_w[2] = {0.5, 0.5};
  int count = 0;
  //Coord
  std::cout << "dofCoord" << std::endl;
  std::cout.setf(std::ios::scientific|std::ios::showpos);
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << std::setprecision(16) << "real3(" << lgl[ii] << ", " << lgl[jj] << "," << lgl[kk] << "), \n";
      }
    }
  }
  count = 0;
  std::cout << "volIntPoints" << std::endl;
  //VolIntPoints
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2[ii] << ", " << gl2[jj] << "," << gl2[kk] << "),\n";
        volIntPoint[count][0] = gl2[ii];
        volIntPoint[count][1] = gl2[jj];
        volIntPoint[count][2] = gl2[kk];
        ++count;
      }
    }
  }
  std::cout << "faceIntPoints" << std::endl;
  count = 0;
  //Face Int Points;
  //I-direction: -1
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      std::cout << "real3(" << std::setprecision(16) << -1.0 << ", " << gl2[jj] << "," << gl2[kk] << "),\n";
      faceIntPoint[count][0] = -1;
      faceIntPoint[count][1] = gl2[jj];
      faceIntPoint[count][2] = gl2[kk];
        ++count;
    }
  }
  //I-direction: 1
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      std::cout << "real3(" << std::setprecision(16) <<  1.0 << ", " << gl2[jj] << "," << gl2[kk] << "),\n";
      faceIntPoint[count][0] = 1;
      faceIntPoint[count][1] = gl2[jj];
      faceIntPoint[count][2] = gl2[kk];
        ++count;
    }
  }
  //J-direction: -1
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 1; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2[ii] << ", " << -1.0 << "," << gl2[kk] << "),\n";
        faceIntPoint[count][0] = gl2[ii];
        faceIntPoint[count][1] = -1;
        faceIntPoint[count][2] = gl2[kk];
        ++count;
      }
    }
  }
  //J-direction: 1
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 1; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2[ii] << ", " << 1.0 << "," << gl2[kk] << "),\n";
        faceIntPoint[count][0] = gl2[ii];
        faceIntPoint[count][1] = 1;
        faceIntPoint[count][2] = gl2[kk];
        ++count;
      }
    }
  }
  //K-direction: -1
  for (int kk = 0; kk < 1; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2[ii] << ", " << gl2[jj] << "," << -1.0 << "),\n";
        faceIntPoint[count][0] = gl2[ii];
        faceIntPoint[count][1] = gl2[jj];
        faceIntPoint[count][2] = -1;
        ++count;
      }
    }
  }
  //K-direction: 1
  for (int kk = 0; kk < 1; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2[ii] << ", " << gl2[jj] << "," << 1.0 << "),\n";
        faceIntPoint[count][0] = gl2[ii];
        faceIntPoint[count][1] = gl2[jj];
        faceIntPoint[count][2] = 1;
        ++count;
      }
    }
  }

  //std::cout << "halfFaceIntPoints" << std::endl;
  //Face Int Points;
  //I-direction: -1
  for (int kk = 0; kk < 4; ++kk){
    for (int jj = 0; jj < 4; ++jj){
      std::cout << "real3(" << std::setprecision(16) << -1.0 << ", " << gl2_half[jj] << "," << gl2_half[kk] << "),\n";
      faceIntPoint[count][0] = -1.0;
      faceIntPoint[count][1] = gl2_half[jj];
      faceIntPoint[count][2] = gl2_half[kk];
        ++count;
    }
  }
  //I-direction: 1
  for (int kk = 0; kk < 4; ++kk){
    for (int jj = 0; jj < 4; ++jj){
      std::cout << "real3(" << std::setprecision(16) << 1.0 << ", " << gl2_half[jj] << "," << gl2_half[kk] << "),\n";
      faceIntPoint[count][0] = 1.0;
      faceIntPoint[count][1] = gl2_half[jj];
      faceIntPoint[count][2] = gl2_half[kk];
        ++count;
    }
  }
  //J-direction: -1
  for (int kk = 0; kk < 4; ++kk){
    for (int jj = 0; jj < 1; ++jj){
      for(int ii = 0; ii < 4; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2_half[ii] << ", " << -1.0 << "," << gl2_half[kk] << "),\n";
      faceIntPoint[count][0] = gl2_half[ii];
      faceIntPoint[count][1] = -1.0;
      faceIntPoint[count][2] = gl2_half[kk];
        ++count;
      }
    }
  }
  //J-direction: 1
  for (int kk = 0; kk < 4; ++kk){
    for (int jj = 0; jj < 1; ++jj){
      for(int ii = 0; ii < 4; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2_half[ii] << ", " << 1.0 << "," << gl2_half[kk] << "),\n";
        faceIntPoint[count][0] = gl2_half[ii];
        faceIntPoint[count][1] = 1.0;
        faceIntPoint[count][2] = gl2_half[kk];
        ++count;
      }
    }
  }
  //K-direction: -1
  for (int kk = 0; kk < 1; ++kk){
    for (int jj = 0; jj < 4; ++jj){
      for(int ii = 0; ii < 4; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2_half[ii] << ", " << gl2_half[jj] << "," << -1.0 << "),\n";
        faceIntPoint[count][0] = gl2_half[ii];
        faceIntPoint[count][1] = gl2_half[jj];
        faceIntPoint[count][2] = -1.0;
        ++count;
      }
    }
  }
  //K-direction: 1
  for (int kk = 0; kk < 1; ++kk){
    for (int jj = 0; jj < 4; ++jj){
      for(int ii = 0; ii < 4; ++ii){
        std::cout << "real3(" << std::setprecision(16) << gl2_half[ii] << ", " << gl2_half[jj] << "," << 1.0 << "),\n";
        faceIntPoint[count][0] = gl2_half[ii];
        faceIntPoint[count][1] = gl2_half[jj];
        faceIntPoint[count][2] = 1.0;
        ++count;
      }
    }
  }
  std::cout << "phi at volIntPoint" << std::endl;
  //phi At volIntPoint
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "realVector{\n" ;
        for (int ll = 0; ll < 8; ++ll){
          if((ll+1)%2 == 1) std::cout << "  ";
          double value = phi[ii](volIntPoint[ll][0])*phi[jj](volIntPoint[ll][1])*phi[kk](volIntPoint[ll][2]);
          std::cout << std::setprecision(16) << value << ", ";
          if((ll+1)%2 == 0) std::cout << std::endl;
        }
        std::cout << "},\n";
      }
    }
  }

  std::cout << "phi at faceIntPoint" << std::endl;
  //phi At faceIntPoint
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "realVector{\n" ;
        for (int ll = 0; ll < 120; ++ll){
          if((ll+1)%2 == 1) std::cout << "  ";
          double value = phi[ii](faceIntPoint[ll][0])*phi[jj](faceIntPoint[ll][1])*phi[kk](faceIntPoint[ll][2]);
          std::cout << std::setprecision(16) << value << ", ";
          if((ll+1)%2 == 0) std::cout << std::endl;
        }
        std::cout << "},\n";
      }
    }
  }

  std::cout << "dphi at volIntPoint" << std::endl;
  //dphi At volIntPoint
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3Vector{\n" ;
        for (int ll = 0; ll < 8; ++ll){
          std::cout << "  ";
          double value1 = dphi[ii](volIntPoint[ll][0])*phi[jj](volIntPoint[ll][1])*phi[kk](volIntPoint[ll][2]);
          double value2 = phi[ii](volIntPoint[ll][0])*dphi[jj](volIntPoint[ll][1])*phi[kk](volIntPoint[ll][2]);
          double value3 = phi[ii](volIntPoint[ll][0])*phi[jj](volIntPoint[ll][1])*dphi[kk](volIntPoint[ll][2]);
          std::cout << "real3(" << std::setprecision(16) << value1 << ", " << value2 << "," << value3 << "),\n";
        }
        std::cout << "},\n";
      }
    }
  }

  std::cout << "dphi at faceIntPoint" << std::endl;
  //dphi At faceIntPoint
  for (int kk = 0; kk < 2; ++kk){
    for (int jj = 0; jj < 2; ++jj){
      for(int ii = 0; ii < 2; ++ii){
        std::cout << "real3Vector{\n" ;
        for (int ll = 0; ll < 120; ++ll){
          std::cout << "  ";
          double value1 = dphi[ii](faceIntPoint[ll][0])*phi[jj](faceIntPoint[ll][1])*phi[kk](faceIntPoint[ll][2]);
          double value2 = phi[ii](faceIntPoint[ll][0])*dphi[jj](faceIntPoint[ll][1])*phi[kk](faceIntPoint[ll][2]);
          double value3 = phi[ii](faceIntPoint[ll][0])*phi[jj](faceIntPoint[ll][1])*dphi[kk](faceIntPoint[ll][2]);
          std::cout << "real3(" << std::setprecision(16) << value1 << ", " << value2 << "," << value3 << "),\n";
        }
        std::cout << "},\n";
      }
    }
  }
}
int main()
{
  generateHexa8_LGL2_Martin();
  return 0;
}
