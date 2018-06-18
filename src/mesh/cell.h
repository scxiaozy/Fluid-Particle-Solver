#ifndef FPS_MESH_CELL_H
#define FPS_MESH_CELL_H
#include <vector>
#include <dataType.h>
#include <real3.h>
namespace fps
{
  //Class cell only stores the linear vertices
  struct cell
  {
    public:
      enum cellType{
        hexohedral,
        others
      };
    public:
      cellType type;
    public:
      //Constructor
      cell();
      //Return the type of the cell
    public:
      //Only the index of vertices are stored
      //to save the memory storage.
      //In dynamic mesh system, this information
      //should be updated after each AMR procedure.
      llabelVector vertices;
  };
  using cellVector = std::vector<cell>;
}
#endif
