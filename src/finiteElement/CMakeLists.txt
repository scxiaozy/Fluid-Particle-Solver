cmake_minimum_required (VERSION 3.10)
set (CMAKE_CXX_STANDARD 11)
INCLUDE_DIRECTORIES(
  ../basic/
  ../mesh/
  .
  )

#lib directories
LINK_DIRECTORIES(
  /usr/local/lib
  )

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# 查找当前目录下的所有源文件
# 并将名称保存到 DIR_LIB_SRCS 变量
aux_source_directory(. DIR_LIB_SRCS)
ADD_Library(FPSFE ${DIR_LIB_SRCS})
