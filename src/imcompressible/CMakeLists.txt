cmake_minimum_required (VERSION 3.10)
set (CMAKE_CXX_STANDARD 11)

INCLUDE_DIRECTORIES(
  /usr/local/include
  ../basic/
  ../mesh/
  ../dynamicMesh/
  ../finiteElement
  ../math
  .
  /Users/jianhuapan/Applications/p4est/include
  /Users/jianhuapan/Applications/libconfig/include
  /Users/jianhuapan/Applications/openmpi/include
  /Users/jianhuapan/Applications/lapack/include
  )

#lib directories
LINK_DIRECTORIES(
  /usr/lib
  /usr/local/lib
  /Users/jianhuapan/Applications/p4est/lib
  /Users/jianhuapan/Applications/libconfig/lib
  /Users/jianhuapan/Applications/openmpi/lib
  /Users/jianhuapan/Applications/lapack/lib
  ../finiteElement/
  )

# 查找当前目录下的所有源文件
# 并将名称保存到 DIR_LIB_SRCS 变量
aux_source_directory(. DIR_SRCS)
add_executable(IMCOMP ${DIR_SRCS})
target_link_libraries(IMCOMP z sc p4est config++ FPSFE lapacke)
