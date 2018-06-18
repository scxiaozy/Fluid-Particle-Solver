#ifndef FPS_BASIC_REAL3_H
#define FPS_BASIC_REAL3_H
#include <cmath>
#include <iostream>

#include <dataType.h>
namespace fps
{
  struct real3
  {
    public:
      real x[3];

      real3(real xx = 0, real yy = 0, real zz = 0){
          x[0]=xx;
          x[1]=yy;
          x[2]=zz;
      }

    inline real3 operator/(real x) const;
    inline real3 operator*(real x) const;
    inline friend real3 operator*(real x, const real3&);
    inline friend real  operator*(const real3& x, const real3& y);
    inline friend real3 operator+(const real3& x, const real3& y);
    inline friend real3 operator-(const real3& x, const real3& y);
    inline real3& operator+=(const real3& y);
    inline real3& operator-=(const real3& y);
    inline real3& operator/=(real x);
    inline real3& operator*=(real x);
    inline real3& operator =(const real3& x);
    inline real mod()const;
    inline friend std::ostream& operator<< (std::ostream& out, const real3& x);
    inline real& operator[](int ii);
    inline const real& operator[](int ii)const;
    inline friend real3 crossMult(const real3&, const real3&);
    inline friend real  Det(const real3&, const real3&, const real3&);
  };

  inline real3 real3::operator/(real xx) const
  {
    return real3(x[0]/xx,x[1]/xx,x[2]/xx);
  }
  inline real3 real3::operator*(real xx) const
  {
    return real3(x[0]*xx,x[1]*xx,x[2]*xx);
  }
  inline real3 operator*(real x, const real3& r)
  {
    return real3(x*r.x[0],x*r.x[1],x*r.x[2]);
  }
  inline real3& real3::operator/=(real xx)
  {
    x[0]/=xx;
    x[1]/=xx;
    x[2]/=xx;
    return *this;
  }
  inline real3& real3::operator*=(real xx)
  {
    x[0]*=xx;
    x[1]*=xx;
    x[2]*=xx;
    return *this;
  }
  inline real3 operator+(const real3& x, const real3& y)
  {
    return real3(x.x[0]+y.x[0], x.x[1]+y.x[1], x.x[2]+y.x[2]);
  }
  inline real3 operator-(const real3& x, const real3& y)
  {
    return real3(x.x[0]-y.x[0], x.x[1]-y.x[1], x.x[2]-y.x[2]);
  }
  inline real  operator*(const real3& x, const real3& y)
  {
      return x.x[0]*y.x[0] + x.x[1]*y.x[1]+ x.x[2]*y.x[2];
  }
  inline real3& real3::operator+=(const real3& yy)
  {
    x[0]+=yy.x[0];
    x[1]+=yy.x[1];
    x[2]+=yy.x[2];
    return *this;
  }
  inline real3& real3::operator-=(const real3& yy)
  {
    x[0]-=yy.x[0];
    x[1]-=yy.x[1];
    x[2]-=yy.x[2];
    return *this;
  }

  inline real3& real3::operator=(const real3& xx)
  {
    x[0]=xx.x[0];
    x[1]=xx.x[1];
    x[2]=xx.x[2];
    return *this;
  }
  inline real real3::mod()const
  {
      return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  }
  inline std::ostream& operator<< (std::ostream& out, const real3& x)
  {
      out << "(" << x.x[0] << ", " << x.x[1] << ", " << x.x[2] <<")"; return out;
  }
  inline real& real3::operator[](int ii){return x[ii];}
  inline const real& real3::operator[](int ii)const{return x[ii];}
  inline real3 crossMult(const real3& l, const real3& r)
  {
    return real3(l[1]*r[2]-l[2]*r[1], l[2]*r[0]-l[0]*r[2], l[0]*r[1]-l[1]*r[0]);
  }
  inline real Det(const real3& x, const real3& y, const real3 z)
  {
    return x*crossMult(y,z);
  }

  inline real3 fabs(const real3& x)
  {
    return real3(std::fabs(x[0]),std::fabs(x[1]),std::fabs(x[2]));
  }
  using real3Vector = std::vector<real3>;
  using real3VectorVector = std::vector<std::vector<real3> >;
}
#endif
