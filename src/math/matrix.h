#ifndef FPS_MATH_MATRIX_H
#define FPS_MATH_MATRIX_H
#include <iostream>
#include <cmath>
#include <iomanip>

#include <dataType.h>

#include <mpi.h>
namespace fps
{
  namespace math{
    class matrix;
  }

  void unserialize(math::matrix& m, char* ptr);

  void serialize(const math::matrix& m, char* ptr);

  namespace math
  {

    class matrix
    {
      private:
        size_t m_;
        size_t n_;
        size_t capacity_;
        real* matrix_;
      public:
        ~matrix();
        matrix(size_t m = 0, size_t n = 0);
        matrix(const matrix& MC);
        matrix(float* data, size_t m, size_t n);
        matrix(double* data, size_t m, size_t n);
        void assign(float* data, size_t m, size_t n);
        void assign(double* data, size_t m, size_t n);
        inline void init(real data){
          for (int ii = 0; ii < m_*n_; ++ii){
            matrix_[ii] = data;
          }
        };
        inline void init(real data, size_t m, size_t n);
        inline void swap(matrix& m);
        const std::pair<size_t,size_t> size()const{return std::pair<size_t,size_t>(m_,n_);}
        inline real* operator[](int i);
        inline real* operator[](int i) const;
        inline matrix& reshape(size_t m, size_t n);
        inline matrix& resize(size_t m, size_t n);
        matrix& operator = (const matrix& m);
        void print(std::ostream &out)const;
        matrix operator*(const matrix& m2) const;
        friend void fps::unserialize(math::matrix& m, char* ptr);
        friend void fps::serialize(const math::matrix& m, char* ptr);
        size_t capacity(){return capacity_;}
    };

    //Constructor
    
    matrix::matrix(size_t m, size_t n):m_(m),n_(n)
    {
      matrix_ = new real[m_*n_];
      memset(matrix_, 0, sizeof(real)*m_*n_);
      capacity_ = m_*n_;
    }
    matrix::matrix(float* data, size_t m, size_t n):m_(m),n_(n)
    {
      matrix_ = new real[m*n];
      for (int ii = 0; ii < m*n; ++ii){
        matrix_[ii] = data[ii];
      }
      capacity_ = m_*n_;
    }

    matrix::matrix(double* data, size_t m, size_t n):m_(m),n_(n)
    {
      matrix_ = new real[m*n];
      for (int ii = 0; ii < m*n; ++ii){
        matrix_[ii] = data[ii];
      }
      capacity_ = m_*n_;
    }

    void matrix::assign(float* data, size_t m, size_t n)
    {
      if(capacity_ < m*n){
        delete[] matrix_;
        matrix_ = new real[m*n];
        capacity_ = m*n;
      }
      m_ = m;
      n_ = n;
      for (int ii = 0; ii < m*n; ++ii){
        matrix_[ii] = data[ii];
      }
    }
    inline void matrix::swap(matrix& m)
    {
      size_t temp;
      temp = m_;
      m_ = m.m_;
      m.m_ = temp;

      temp = n_;
      n_ = m.n_;
      m.n_ = temp;

      temp = capacity_;
      capacity_ = m.capacity_;
      m.capacity_ = temp;

      real *tempPtr = matrix_;
      matrix_ = m.matrix_;
      m.matrix_ = tempPtr;
    }
    void matrix::assign(double* data, size_t m, size_t n)
    {
      if(capacity_ < m*n){
        delete[] matrix_;
        matrix_ = new real[m*n];
        capacity_ = m*n;
      }
      m_ = m;
      n_ = n;
      for (int ii = 0; ii < m*n; ++ii){
        matrix_[ii] = data[ii];
      }
    }
    inline void matrix::init(real data, size_t m, size_t n){
      if(capacity_ < m*n){
        delete[] matrix_;
        matrix_ = new real[m*n];
        capacity_ = m*n;
      }
      m_ = m;
      n_ = n;
      for (int ii = 0; ii < m*n; ++ii){
        matrix_[ii] = data;
      }
    }
    
    matrix::matrix(const matrix& m):m_(m.m_),n_(m.n_)
    {
      matrix_ = new real[m.m_*m.n_];
      memcpy(matrix_, m.matrix_, sizeof(real)*m.n_*m.m_);
      capacity_ = m_*n_;
    }

    
    matrix::~matrix()
    {
      delete[] matrix_;
    }

    
    inline real* matrix::operator[](const int i)
    {
      return matrix_+i*n_;
    }

    
    inline real* matrix::operator[](const int i) const
    {
      return matrix_+i*n_;
    }

    
    inline matrix& matrix::reshape(size_t m, size_t n)
    {
      assert(m_*n_ == m*n);
      m_ = m;
      n_ = n;
      return *this;
    }

    
    inline matrix& matrix::resize(size_t m, size_t n)
    {
      if(capacity_ < m*n){
        real* temp = new real[m*n];
        memset(temp, 0, sizeof(real)*m*n);
        size_t len = m_*n_;
        memcpy(temp, matrix_, len*sizeof(real));
        delete[] matrix_;
        matrix_ = temp;
        capacity_ = m*n;
      }else{
        memset(matrix_+m*n, 0, (capacity_-m*n)*sizeof(real));
      }
      m_ = m;
      n_ = n;
      return *this;
    }

    
    matrix& matrix::operator=(const matrix& m)
    {
      if(capacity_ >= m.m_*m.n_){
        memcpy(matrix_, m.matrix_, sizeof(real)*m.n_*m.m_);
      }else{
        delete[] matrix_;
        matrix_ = new real[m.m_*m.n_];
        memcpy(matrix_, m.matrix_, sizeof(real)*m.n_*m.m_);
        capacity_ = m.m_*m.n_;
      }
      m_ = m.m_;
      n_ = m.n_;
      return *this;
    }

    
    std::ostream&  operator << (std::ostream &out, const matrix& m)
    {
      size_t m_ = m.size().first;
      size_t n_ = m.size().second;
      out << "\n{\n";
      for (size_t ii = 0; ii < m_; ++ii){
        out << "  {";
        for (size_t jj = 0; jj < n_; ++jj){
          out << m[ii][jj] << " ";
          if(jj!=n_-1) out << ",";
        }
        out << "}";
        if(ii!=m_-1) out << ",\n";
        out << std::endl;
      }
      out << "\n}\n";
      return out;
    }

    
    void matrix::print(std::ostream &out)const
    {
      out << "\n{\n";
      for (int ii = 0; ii < m_; ++ii){
        out << "  {";
        for (int jj = 0; jj < n_; ++jj){
          out << matrix_[ii*n_+jj] << " ";
          if(jj!=n_-1) out << ",";
        }
        out << "}";
        if(ii!=m_-1) out << ",\n";
        out << std::endl;
      }
      out << "\n}\n";
    }
    
    matrix matrix::operator*(const matrix& m2)const
    {
      matrix result(this->m_, m2.n_);
      if(n_!=m2.m_) {
        std::cerr << "Dimension mismatch of matrix multiply" << std::endl;
        exit(-1);
      }
      for (int ii = 0; ii < m_; ++ii){
        for (int jj = 0; jj < m2.n_; ++jj){
          result[ii][jj] = 0.0;
          for (int kk = 0; kk < n_; ++kk){
            result[ii][jj] += (*this)[ii][kk]*m2[kk][jj];
          }
        }
      }
      return result;
    }
  }
  size_t serializedSize(const math::matrix& m)
  {
    std::pair<size_t,size_t> ss = m.size();
    return sizeof(size_t)*2 + sizeof(real)*(ss.first*ss.second);
  }

  void serialize(const math::matrix& m, char* ptr)
  {
    std::pair<size_t, size_t> ss = m.size();
    ((size_t*)ptr)[0] = ss.first;
    ((size_t*)ptr)[1] = ss.second;
    ptr+= 2*sizeof(size_t)/sizeof(char);
    for (size_t ii = 0; ii < ss.first; ++ii){
      for (size_t jj = 0; jj < ss.second; ++jj){
        ((real*)ptr)[ii*ss.second+jj] = m[ii][jj];
      }
    }
  }

  void unserialize(math::matrix& m, char* ptr)
  {
    m.m_ = ((size_t*)ptr)[0];
    m.n_ = ((size_t*)ptr)[1];
    delete[] m.matrix_;
    m.matrix_ = new real[m.m_*m.n_];
    ptr+= 2*sizeof(size_t)/sizeof(char);
    for (size_t ii = 0; ii < m.m_; ++ii){
      for (size_t jj = 0; jj < m.n_; ++jj){
        m[ii][jj] = ((real*)ptr)[0];
        ptr+=sizeof(real)/sizeof(char);
      }
    }
    m.capacity_ = m.m_*m.n_;
  }
}
#endif
