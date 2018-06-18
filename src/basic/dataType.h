#ifndef FPS_BASIC_DATATYPE_H
#define FPS_BASIC_DATATYPE_H
#include <vector>
#include <cstdint>
#include "config.inc"
namespace fps
{
#ifdef LINT64
  using llabel = int64_t;
#elif defined LINT32
  using llabel = int32_t;
#elif
#error "Either LINT32 or LINT64 should be defiend for llabel"
#endif

#ifdef GINT64
  using glabel = int64_t;
#elif defined GINT32
  using glabel = int32_t;
#elif
#error "Either GINT32 or GINT64 should be defiend for glabel"
#endif

#ifdef FLOAT32
  using real = float;
#define REAL_MAX  FLT_MAX;
#define REAL_MIN  FLT_MIN;
#elif defined FLOAT64
  using real = double;
#define REAL_MAX  DBL_MAX;
#define REAL_MIN  DBL_MIN;
#elif
#error "Either FLOAT32 or FLOAT64 should be defiend for real"
#endif

  using llabelVector       = std::vector<llabel>      ;
  using glabelVector       = std::vector<glabel>      ;
  using llabelVectorVector = std::vector<llabelVector>;
  using glabelVectorVector = std::vector<glabelVector>;
  using llabelPair         = std::pair<llabel,llabel> ;
  using glabelPair         = std::pair<glabel,glabel> ;
  using realVector         = std::vector<real>        ;
  using realVectorVector   = std::vector<realVector>  ;
  using intPair            = std::pair<int,int>       ;
}
#endif
