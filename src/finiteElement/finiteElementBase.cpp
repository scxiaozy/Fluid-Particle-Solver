#include <finiteElementBase.h>
namespace fps{
  size_t serializedSize(const finiteElementBase& elem)
  {
    return
      sizeof(int)*2
      +sizeof(size_t)*4
      +sizeof(real)*3
        *(elem.vertices.size()
            + elem.solutionPhyCoord.size()
            + elem.volIntPhyCoord.size())
      +sizeof(real)*elem.volIntWxJ.size();
  }
  void serialize(const finiteElementBase& elem, char* ptr){
    ((int*) ptr)[0] = elem.type;
    ((int*) ptr)[1] = elem.order;
    ptr += 2*sizeof(int)/sizeof(char);
    size_t ss = ((size_t*)ptr)[0] = elem.vertices.size();
    ptr += sizeof(size_t)/sizeof(char);
    for (size_t ii = 0; ii < ss; ++ii){
      ((real3*)ptr)[ii] = elem.vertices[ii];
    }
    ptr += 3*ss*sizeof(real)/sizeof(char);

    ss = ((size_t*)ptr)[0] = elem.solutionPhyCoord.size();
    ptr += sizeof(size_t)/sizeof(char);
    for (size_t ii = 0; ii < ss; ++ii){
      ((real3*)ptr)[ii] = elem.solutionPhyCoord[ii];
    }
    ptr += 3*ss*sizeof(real)/sizeof(char);

    ss = ((size_t*)ptr)[0] = elem.volIntWxJ.size();
    ptr += sizeof(size_t)/sizeof(char);
    for (size_t ii = 0; ii < ss; ++ii){
      ((real*)ptr)[ii] = elem.volIntWxJ[ii];
    }
    ptr += ss*sizeof(real)/sizeof(char);

    ss = ((size_t*)ptr)[0] = elem.volIntPhyCoord.size();
    ptr += sizeof(size_t)/sizeof(char);
    for (size_t ii = 0; ii < ss; ++ii){
      ((real3*)ptr)[ii] = elem.volIntPhyCoord[ii];
    }
    ptr += 3*ss*sizeof(real)/sizeof(char);
  }

  void unserialize(finiteElementBase& elem, char* ptr)
  {
    elem.type  = ((int*) ptr)[0];
    elem.order = ((int*) ptr)[1];
    ptr += 2*sizeof(int)/sizeof(char);

    size_t ss = ((size_t*)ptr)[0];
    ptr += sizeof(size_t)/sizeof(char);
    elem.vertices.reserve(ss);
    elem.vertices.assign((real3*)ptr, ((real3*)ptr)+ss);
    ptr += 3*ss*sizeof(real)/sizeof(char);

    ss = ((size_t*)ptr)[0];
    elem.solutionPhyCoord.reserve(ss);
    ptr += sizeof(size_t)/sizeof(char);
    elem.solutionPhyCoord.assign((real3*)ptr, ((real3*)ptr)+ss);
    ptr += 3*ss*sizeof(real)/sizeof(char);

    ss = ((size_t*)ptr)[0];
    elem.volIntWxJ.reserve(ss);
    ptr += sizeof(size_t)/sizeof(char);
    elem.volIntWxJ.assign((real*)ptr, ((real*)ptr)+ss);
    ptr += ss*sizeof(real)/sizeof(char);

    ss = ((size_t*)ptr)[0];
    elem.volIntPhyCoord.reserve(ss);
    ptr += sizeof(size_t)/sizeof(char);
    elem.volIntPhyCoord.assign((real3*)ptr, ((real3*)ptr)+ss);
    ptr += 3*ss*sizeof(real)/sizeof(char);
  }
}
