#ifndef FPS_DYNAMICMESH_DYNAMICMESH_UPDATEMAPVECTOR_H
#define FPS_DYNAMICMESH_DYNAMICMESH_UPDATEMAPVECTOR_H
#include <dynamicMesh.h>
namespace fps
{
  template<class T>
  void dynamicMesh<T>::updateMapVector()
  {
    llabel nCells = p8estPtr_->local_num_quadrants;
    llabel nGhostCells = p8estGhostPtr_->ghosts.elem_count;
    mapVector_.resize(nCells+nGhostCells);
    int count = 0;
    invMapVector_.resize(oldNCells_);
    std::vector<std::vector<glabel> > temp1(oldNCells_);
    std::vector<std::vector<int> > temp2(oldNCells_);
    std::vector<int> tag(oldNCells_);

    for (int ii = p8estPtr_->first_local_tree; ii <= p8estPtr_->last_local_tree; ++ii){
      p8est_tree_t *tree = p8est_tree_array_index(p8estPtr_->trees, ii);
      for (unsigned int jj = 0; jj < tree->quadrants.elem_count; ++jj){
        p8est_quadrant_t *q = p8est_quadrant_array_index(&(tree->quadrants), jj);
        p8estData *data = (p8estData*)q->p.user_data;
        mapVector_[count] = data->p;
        if(mapVector_[count].fineCoarsenTag==0){
          glabel oldLocalIdx = mapVector_[count].fcTag.equal.oldLocalIdx;
          temp1[oldLocalIdx] .push_back(count);
          tag[oldLocalIdx] = 0;
        }else if(mapVector_[count].fineCoarsenTag==-1){
          for (int h = 0; h < P8EST_CHILDREN; ++h){
            glabel oldLocalIdx = mapVector_[count].fcTag.coarser.oldLocalIdx[h];
            temp1[oldLocalIdx].push_back(count);
            tag[oldLocalIdx] = -1;
          }
        }else if(mapVector_[count].fineCoarsenTag==+1){
           glabel oldLocalIdx = mapVector_[count].fcTag.finer.oldLocalIdx;
           int    cid           = mapVector_[count].fcTag.finer.childID;
           temp1[oldLocalIdx].push_back(count);
           temp2[oldLocalIdx].push_back(cid);
           tag[oldLocalIdx] = +1;
        }
        ++count;
      }
    }
    for (int ii = 0; ii < this->oldNCells_; ++ii){
      if(tag[ii] == 0){
        invMapVector_[ii].fineCoarsenTag = 0;
        invMapVector_[ii].fcTag.equal.newLocalIdx = temp1[ii][0];
      }else if(tag[ii] == 1){
        invMapVector_[ii].fineCoarsenTag = 1;
        assert(temp1[ii].size()==P8EST_CHILDREN);
        for (int h = 0; h < P8EST_CHILDREN; ++h){
          invMapVector_[ii].fcTag.finer.newLocalIdx[h] = temp1[ii][h];
          invMapVector_[ii].fcTag.finer.childID[h] = temp2[ii][h];
        }
      }else if(tag[ii] == -1){
        invMapVector_[ii].fineCoarsenTag = -1;
        invMapVector_[ii].fcTag.coarser.newLocalIdx = temp1[ii][0];
      }
    }
    for (int ii = 0; ii < p8estGhostPtr_->ghosts.elem_count; ++ii){
      p8estData *data = ghostData_ + ii;
      mapVector_[ii+nCells] = data->p;
    }
    MPI_Barrier(primitiveMesh::mpiComm_);
  }

}
#endif
