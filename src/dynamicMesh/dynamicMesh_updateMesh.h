#ifndef FPS_DYNAMICMESH_DYNAMICMESH_UPDATEMESH_H
#define FPS_DYNAMICMESH_DYNAMICMESH_UPDATEMESH_H
//Header files for std library
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <utility>
//Header file for p8est
#include <p8est.h>
//Header files for FPS
#include <dynamicMesh.h>
#include <sync.h>

namespace fps
{
  struct facePtrEqualUpdateMesh{
    bool operator() (const face* lhs, const face* rhs) const
    {
      return
        (lhs->gcellIdx  == rhs->gcellIdx)
        //`firstSide` used to eliminate the
        //ambiguity of self periodic cell
      &&(lhs->firstSide == rhs->firstSide)
        ;
    }
  };
  struct facePtrHashUpdateMesh{
    std::size_t operator()(const face* key)const
    {
      using std::size_t;
      using std::hash;
      return ((hash<glabel>()(key->gcellIdx.first)
            ^ (hash<glabel>()(key->gcellIdx.second) << 1)) >> 1)
            ^ (hash<int>()(key->firstSide) << 1);
    }
  };
  static unsigned long hashpjw(char *arKey, unsigned int nKeyLength)
  {
    unsigned long h = 0, g;
    char *arEnd=arKey+nKeyLength;

    while (arKey < arEnd) {
      h = (h << 4) + *arKey++;
      if ((g = (h & 0xF0000000))){
        h = h ^ (g >> 24);
        h = h ^ g;
      }
    }
    return h;
  }

  template <class T>
  struct meshUpdateContex{
    dynamicMesh<T>* meshPtr_;
    std::unordered_map<face*, llabel, facePtrHashUpdateMesh, facePtrEqualUpdateMesh>* mapPtr_;
    std::vector<llabel> *cellTag;
    std::vector<llabel> *faceTag;
  };

  template<class T>
  void updateMeshIterVol(p8est_iter_volume_info_t *info, void *user_data)
  {
    //Init the
    dynamicMesh<T>*mesh
      = ((meshUpdateContex<T>*) user_data)->meshPtr_;

    p4est_locidx_t      jl;
    p8est_tree_t       *tree;
    p8est_quadrant_t *q;
    tree = p8est_tree_array_index (info->p4est->trees, info->treeid);
    jl = info->quadid + tree->quadrants_offset;
    q = info->quad;
    typename dynamicMesh<T>::p8estData* p8estdataPtr
      = (typename dynamicMesh<T>::p8estData*)q->p.user_data;

    if(!p8estdataPtr->p.fineCoarsenTag){
      mesh->volElems_[jl] = mesh->oldvolElems_[p8estdataPtr->p.fcTag.equal.oldLocalIdx];
    }
    else{
      const std::vector<llabel>& vIdx = mesh->cellToHiOrderNodes[info->treeid];
      //The coordinates in the reference frame
      //Should be [-1,1]\times[-1,1]\times[-1,1]
      std::vector<real3> vvrefFrame;
      vvrefFrame.reserve(8);
      for (int kk = 0; kk < 2; ++kk){
        for (int jj = 0; jj < 2; ++jj){
          for (int ii = 0; ii < 2; ++ii){
            vvrefFrame.push_back(
              real3(
                2.0*real(q->x + ii*P8EST_QUADRANT_LEN(q->level))/P8EST_ROOT_LEN-1.0,
                2.0*real(q->y + jj*P8EST_QUADRANT_LEN(q->level))/P8EST_ROOT_LEN-1.0,
                2.0*real(q->z + kk*P8EST_QUADRANT_LEN(q->level))/P8EST_ROOT_LEN-1.0
                )
              );
          }
        }
      }
      std::vector<real3> vvphyFrame(vvrefFrame);
      T::obtainCoord(mesh->highOrderNodes, vIdx, vvrefFrame, vvphyFrame);
      mesh->volElems_[jl].initiate(vvphyFrame);
    }
  }

  template<class T>
  void updateMeshIterFace(p8est_iter_face_info_t *info, void *user_data)
  {
    int                 h;
    int                 swapsides;
    dynamicMesh<T>*mesh
      = ((meshUpdateContex<T>*) user_data)->meshPtr_;
    std::unordered_map<face*, llabel, facePtrHashUpdateMesh, facePtrEqualUpdateMesh>* mapPtr_
      = ((meshUpdateContex<T>*) user_data)->mapPtr_;
    std::vector<llabel>& cellTag
      = (*((meshUpdateContex<T>*) user_data)->cellTag);
    std::vector<llabel>& faceTag
      = *(((meshUpdateContex<T>*) user_data)->faceTag);

    p4est_locidx_t      jl, jl2, jls[P8EST_HALF];
    //p4est_locidx_t      in_qtoq, halfindex;
    //p4est_locidx_t     *halfentries;
    p8est_tree_t       *tree;
    p8est_iter_face_side_t *side, *side2, *tempside;
    p4est_locidx_t     local_num_quadrants = info->p4est->local_num_quadrants;
    p8est_quadrant_t *q;
    p8est_quadrant_t *q2;
    p8est_quadrant_t *qs[P8EST_HALF];

    if (info->sides.elem_count == 1) {
      /* this face is on an outside boundary of the forest */
      P4EST_ASSERT (info->orientation == 0);
      P4EST_ASSERT (info->tree_boundary);
      side = (p8est_iter_face_side_t *) sc_array_index (&info->sides, 0);
      P4EST_ASSERT (0 <= side->treeid &&
                    side->treeid < info->p4est->connectivity->num_trees);
      P4EST_ASSERT (0 <= side->face && side->face < P4EST_FACES);
      P4EST_ASSERT (!side->is_hanging && !side->is.full.is_ghost);
      tree = p8est_tree_array_index (info->p4est->trees, side->treeid);

      jl = side->is.full.quadid + tree->quadrants_offset;
      P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);

      q = side->is.full.quad;
      typename dynamicMesh<T>::p8estData* p8estdataPtr
        = (typename dynamicMesh<T>::p8estData*)q->p.user_data;


      //face information
      face f;
      f.lcellIdx.first  = jl;
      f.lcellIdx.second = -1;
      f.gcellIdx.first  = jl + mesh->p8estPtr_->global_first_quadrant[mesh->mpiRank()];
      f.gcellIdx.second = -1;
      f.firstSide  = side->face;
      f.secondSide = side->face;
      f.orientation = 0;
      f.firstFaceChildIdx = -1;
      f.secondFaceChildIdx = -1;
      f.firstCellType = T::refFrame::type;
      f.secondCellType = T::refFrame::type;
      f.isHanging = 0;
      f.patchIdx.first = p8estdataPtr->patchIdx[side->face];
      f.patchIdx.second = -1;

      mesh->p8estCece_[jl].ne[side->face].finer = 0;
      mesh->p8estCece_[jl].ne[side->face].index[0] = -1;
      mesh->p8estCece_[jl].ne[side->face].pidx[0] = mesh->nFaces_;

      if(!p8estdataPtr->p.fineCoarsenTag){
        int oldGlobalIdx = p8estdataPtr->p.fcTag.equal.oldGlobalIdx;
        face tempf;
        tempf.gcellIdx.first  = oldGlobalIdx;
        tempf.gcellIdx.second = -1;
        tempf.firstSide = side->face;
        auto it  = mapPtr_->find(&(tempf));
        assert(it!=mapPtr_->end());
        llabel oldfaceIdx = it->second;
        face ttf = mesh->oldfaces_[oldfaceIdx];
        ttf.lcellIdx  = f.lcellIdx;
        ttf.gcellIdx  = f.gcellIdx;
        f = ttf;

        //For later initialization of faceElement
        faceTag.push_back(oldfaceIdx);
        //mesh->faceElem_.push_back(mesh->oldfaceElem_[oldfaceIdx]);
      }else{
        faceTag.push_back(-1);
        //mesh->faceElem_.push_back(
            //faceElement(T::refFrame::nSolutionPoints, f, mesh->volElem_[jl], mesh->volElem_[jl])
            //);
      }
      if(f.patchIdx.first > 0)
        mesh->boundaryFaceIdxList_[f.patchIdx.first].push_back(mesh->nFaces_);
      mesh->faces_.push_back(f);
      ++mesh->nFaces_;
    }
    else {
      /* this face is between two quadrants */
      P4EST_ASSERT (info->orientation == 0 || info->tree_boundary);
      P4EST_ASSERT (info->sides.elem_count == 2);
      side = (p8est_iter_face_side_t *) sc_array_index (&info->sides, 0);
      side2 = (p8est_iter_face_side_t *) sc_array_index (&info->sides, 1);
      P4EST_ASSERT (info->tree_boundary || side->treeid == side2->treeid);
      P4EST_ASSERT (!side->is_hanging || !side2->is_hanging);
      if (!side->is_hanging && !side2->is_hanging){
        /* same-size face neighbors */
        P4EST_ASSERT (!side->is.full.is_ghost || !side2->is.full.is_ghost);
        p4est_gloidx_t jlg;
        p4est_gloidx_t jlg2;

        /* determine both quadrant numbers */
        if (!side->is.full.is_ghost) {
          tree = p8est_tree_array_index (info->p4est->trees, side->treeid);
          jl  = side->is.full.quadid + tree->quadrants_offset;
          P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
        }
        else {
          P4EST_ASSERT (side->is.full.quad !  = NULL);
          P4EST_ASSERT (side->is.full.quadid >= 0);
          jl = local_num_quadrants + side->is.full.quadid;
        }

        if (!side2->is.full.is_ghost) {
          tree = p8est_tree_array_index (info->p4est->trees, side2->treeid);
          jl2 = side2->is.full.quadid + tree->quadrants_offset;
          P4EST_ASSERT (0 <= jl2 && jl2 < local_num_quadrants);
        }
        else {
          P4EST_ASSERT (side2->is.full.quad != NULL);
          P4EST_ASSERT (side2->is.full.quadid >= 0);
          jl2 = local_num_quadrants + side2->is.full.quadid;
        }

        //Global index for the quadrant
        jlg = mesh->cellGlobalIdx()[jl];
        jlg2 = mesh->cellGlobalIdx()[jl2];

        q = side->is.full.quad;
        q2 = side2->is.full.quad;
        typename dynamicMesh<T>::p8estData* p8estdataPtr;

        if(jl<local_num_quadrants){
          p8estdataPtr = (typename dynamicMesh<T>::p8estData*)q->p.user_data;
        }
        else
          p8estdataPtr = mesh->ghostData_ + jl - local_num_quadrants;

        typename dynamicMesh<T>::p8estData* p8estdataPtr2;
        if(jl2<local_num_quadrants)
          p8estdataPtr2 = (typename dynamicMesh<T>::p8estData*)q2->p.user_data;
        else
          p8estdataPtr2 = mesh->ghostData_ + jl2 - local_num_quadrants;

        //Primitive mesh information
        face f;

        int lor;
        //The f.lcellIdx.first should always not greater than the second
        //for nonhaning face
        if(jlg <= jlg2
            && p8estdataPtr->patchIdx[side->face]
            <= p8estdataPtr2->patchIdx[side2->face]){
          lor = 0;
          f.lcellIdx.first  = jl;
          f.lcellIdx.second = jl2;
          f.gcellIdx.first  = jlg;
          f.gcellIdx.second = jlg2;
          f.firstSide  = side->face;
          f.secondSide = side2->face;
          f.orientation = info->orientation;
          f.firstFaceChildIdx = -1;
          f.secondFaceChildIdx = -1;
          f.firstCellType  = T::refFrame::type;
          f.secondCellType = T::refFrame::type;
          f.isHanging = 0;
          f.patchIdx.first = p8estdataPtr->patchIdx[side->face];
          f.patchIdx.second = p8estdataPtr2->patchIdx[side2->face];
        }else{
          lor = 1;
          f.lcellIdx.first  = jl2;
          f.lcellIdx.second = jl;
          f.gcellIdx.first  = jlg2;
          f.gcellIdx.second = jlg;
          f.firstSide  = side2->face;
          f.secondSide = side->face;
          f.orientation = info->orientation;
          f.firstFaceChildIdx = -1;
          f.secondFaceChildIdx = -1;
          f.firstCellType = T::refFrame::type;
          f.secondCellType = T::refFrame::type;
          f.isHanging = 0;
          f.patchIdx.first = p8estdataPtr2->patchIdx[side2->face];
          f.patchIdx.second = p8estdataPtr->patchIdx[side->face];
        }

        if(jl<local_num_quadrants){
          mesh->p8estCece_[jl].ne[side->face].finer = 0;
          mesh->p8estCece_[jl].ne[side->face].index[0] = jl2;
          mesh->p8estCece_[jl].ne[side->face].pidx[0] = mesh->nFaces_;
        }
        if(jl2<local_num_quadrants){
          mesh->p8estCece_[jl2].ne[side2->face].finer = 0;
          mesh->p8estCece_[jl2].ne[side2->face].index[0] = jl;
          mesh->p8estCece_[jl2].ne[side2->face].pidx[0] = mesh->nFaces_;
        }

        if(!p8estdataPtr->p.fineCoarsenTag
            &&!p8estdataPtr2->p.fineCoarsenTag){

          p4est_gloidx_t oldGlobalIdx;
          p4est_gloidx_t oldGlobalIdx2;

          //oldLocalIdx used later should be a non ghost cell
          //if(jl < local_num_quadrants)
            //oldLocalIdx  = p8estdataPtr->p.fcTag.equal.oldLocalIdx;
          //else
            //oldLocalIdx = p8estdataPtr2->p.fcTag.equal.oldLocalIdx;
          if(!lor){
            oldGlobalIdx   = p8estdataPtr->p.fcTag.equal.oldGlobalIdx;
            oldGlobalIdx2  = p8estdataPtr2->p.fcTag.equal.oldGlobalIdx;
          }else{
            oldGlobalIdx2  = p8estdataPtr->p.fcTag.equal.oldGlobalIdx;
            oldGlobalIdx   = p8estdataPtr2->p.fcTag.equal.oldGlobalIdx;
          }

          face tempf;
          tempf.gcellIdx.first  = oldGlobalIdx;
          tempf.gcellIdx.second = oldGlobalIdx2;
          tempf.firstSide = f.firstSide;
          auto it  = mapPtr_->find(&(tempf));
          assert(it!=mapPtr_->end());
          llabel oldfaceIdx = it->second;

          face ttf = mesh->oldfaces_[oldfaceIdx];
          ttf.lcellIdx  = f.lcellIdx;
          ttf.gcellIdx  = f.gcellIdx;
          f = ttf;
          faceTag.push_back(oldfaceIdx);
          //mesh->faceElem_.push_back(mesh->oldfaceElem_[oldfaceIdx]);
        }else{
          faceTag.push_back(-1);
            //mesh->faceElem_.push_back(
                //faceElement(T::refFrame::nSolutionPoints, f,
                  //mesh->volElem_[f.lcellIdx.first], mesh->volElem_[f.lcellIdx.second])
            //);
        }
        mesh->faces_.push_back(f);
        ++mesh->nFaces_;
      }
      else {
        /* one of the faces is hanging, rename so it's always side2 */
        swapsides = side->is_hanging;
        if (swapsides) {
          tempside = side;
          side = side2;
          side2 = tempside;
        }
        P4EST_ASSERT (!side->is_hanging && side2->is_hanging);
        p4est_gloidx_t jlg;
        p4est_gloidx_t jlg2[P8EST_HALF];
        /* determine quadrant number for non-hanging large face */
        if (!side->is.full.is_ghost) {
          tree = p8est_tree_array_index (info->p4est->trees, side->treeid);
          jl = side->is.full.quadid + tree->quadrants_offset;
          P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
        }
        else {
          P4EST_ASSERT (side->is.full.quad != NULL);
          P4EST_ASSERT (side->is.full.quadid >= 0);
          jl = local_num_quadrants + side->is.full.quadid;
        }
        jlg = mesh->cellGlobalIdx()[jl];

        /* determine quadrant numbers for all hanging faces */
        for (h = 0; h < P8EST_HALF; ++h) {
          if (!side2->is.hanging.is_ghost[h]) {
            tree = p8est_tree_array_index (info->p4est->trees, side2->treeid);
            jls[h] = side2->is.hanging.quadid[h] + tree->quadrants_offset;
            P4EST_ASSERT (0 <= jls[h] && jls[h] < local_num_quadrants);
          }
          else {
            P4EST_ASSERT (side2->is.hanging.quad[h] != NULL);
            P4EST_ASSERT (side2->is.hanging.quadid[h] >= 0);
            jls[h] = local_num_quadrants + side2->is.hanging.quadid[h];
          }
          jlg2[h] = mesh->cellGlobalIdx()[jls[h]];
        }
        face f[4];
        q = side->is.full.quad;
        typename dynamicMesh<T>::p8estData* p8estdataPtr;

        if(jl < local_num_quadrants){
          p8estdataPtr = (typename dynamicMesh<T>::p8estData*)q->p.user_data;
        }
        else{
          p8estdataPtr = mesh->ghostData_ + jl - local_num_quadrants;
        }

        const int R[6][6] = {
          {0,1,1,0,0,1},
          {2,0,0,1,1,0},
          {2,0,0,1,1,0},
          {0,2,2,0,0,1},
          {0,2,2,0,0,1},
          {2,0,0,2,2,0},
        };
        const int Q[3][4] = {
          {1,2,5,6},
          {0,3,4,7},
          {0,4,3,7},
        };
        const int P[8][4] = {
          {0,1,2,3},
          {0,2,1,3},
          {1,0,3,2},
          {1,3,0,2},
          {2,0,3,1},
          {2,3,0,1},
          {3,1,2,0},
          {3,2,1,0},
        };
        for (h = 0; h < P8EST_HALF; ++h){
          if(side2->is.hanging.is_ghost[h] && side->is.full.is_ghost) continue;
          qs[h] = side2->is.hanging.quad[h];
          //The left side is always the hanging side
          typename dynamicMesh<T>::p8estData* p8estdataPtr2;
          if(jls[h] < local_num_quadrants)
            p8estdataPtr2 = (typename dynamicMesh<T>::p8estData*)qs[h]->p.user_data;
          else
            p8estdataPtr2 = mesh->ghostData_ + jls[h] - local_num_quadrants;

          f[h].lcellIdx.first  = jl;
          f[h].lcellIdx.second = jls[h];
          f[h].gcellIdx.first  = jlg;
          f[h].gcellIdx.second = jlg2[h];
          f[h].firstSide  = side->face;
          f[h].secondSide = side2->face;
          f[h].orientation = info->orientation;
          //Here since the right face is always hanging, the left face is a sub face of a complete face
          //firstFaceChildIdx denotes the sub-id of this face;
          f[h].firstFaceChildIdx = P[Q[R[side2->face][side->face]][info->orientation]][h];
          f[h].secondFaceChildIdx = -1;
          f[h].firstCellType  = T::refFrame::type;
          f[h].secondCellType = T::refFrame::type;
          f[h].isHanging = 0;
          f[h].patchIdx.first  = p8estdataPtr->patchIdx[side->face];
          f[h].patchIdx.second = p8estdataPtr2->patchIdx[side2->face];
          if(!p8estdataPtr->p.fineCoarsenTag
              && !p8estdataPtr2->p.fineCoarsenTag){
            llabel oldGlobalIdx  = p8estdataPtr->p.fcTag.equal.oldGlobalIdx;
            llabel oldGlobalIdx2 = p8estdataPtr2->p.fcTag.equal.oldGlobalIdx;
            face tempf;
            tempf.gcellIdx.first  = oldGlobalIdx;
            tempf.gcellIdx.second = oldGlobalIdx2;
            tempf.firstSide = f[h].firstSide;
            auto it  = mapPtr_->find(&(tempf));
            assert(it!=mapPtr_->end());
            llabel oldfaceIdx = it->second;
            face ttf = mesh->oldfaces_[oldfaceIdx];
            ttf.lcellIdx  = f[h].lcellIdx;
            ttf.gcellIdx  = f[h].gcellIdx;
            f[h] = ttf;
            faceTag.push_back(oldfaceIdx);
          }else{
            faceTag.push_back(-1);
          }

          if(!side->is.full.is_ghost){
            mesh->p8estCece_[jl].ne[side->face].finer =  1;
            mesh->p8estCece_[jl].ne[side->face].index[h] = jls[h];
            mesh->p8estCece_[jl].ne[side->face].pidx[h] = mesh->nFaces_;
          }
          if(!side2->is.hanging.is_ghost[h]){
            mesh->p8estCece_[jls[h]].ne[side2->face].finer =  0;
            mesh->p8estCece_[jls[h]].ne[side2->face].index[0] = jl;
            mesh->p8estCece_[jls[h]].ne[side2->face].pidx[0] = mesh->nFaces_;
          }
          mesh->faces_.push_back(f[h]);
          ++mesh->nFaces_;
        }
      }
    }
  }

  template<class T>
  void dynamicMesh<T>::updateMesh()
  {
    P4EST_ASSERT(p8est_is_balanced (p8estPtr_, P8EST_CONNECT_FULL));
    llabel nCells      = p8estPtr_->local_num_quadrants;
    llabel nGhostCells = p8estGhostPtr_ ->ghosts.elem_count;


    oldfaceElems1_.swap(faceElems1_);
    if(oversetFaceIntPt)
      oldfaceElems2_.swap(faceElems2_);
    oldvolElems_.swap(volElems_);
    oldfaces_            .swap(this->faces_);
    oldcefa_             .swap(this->cefa_);


    std::unordered_map<face*, llabel, facePtrHashUpdateMesh, facePtrEqualUpdateMesh> facePtrHashMap;
    llabel nFaces = oldfaces_.size();
    for (llabel ii = 0; ii < nFaces; ++ii){
      facePtrHashMap.insert(std::make_pair(&(oldfaces_[ii]),ii));
    }

    //Face number and point number will be determined later
    this->nCells_ = nCells;
    this->nGhostCells_ = nGhostCells;
    this->nFaces_ = 0;
    this->nPoints_ = 0;
    this->cellGlobalIdx_     .resize(nCells+nGhostCells);
    this->ghostCellMPIRank_  .resize(nGhostCells);

    //Local cells' global index
    for (int jl = 0; jl < nCells; ++jl){
      this->cellGlobalIdx_[jl] = p8estPtr_->global_first_quadrant[p8estPtr_->mpirank] + jl;
    }
    //Ghost cells' rank and global index
    int rank = 0;
    for (int jl = 0; jl < nGhostCells; ++jl) {
      p8est_quadrant_t *q = p8est_quadrant_array_index(&(p8estGhostPtr_->ghosts), jl);
      while (p8estGhostPtr_->proc_offsets[rank + 1] <= jl) {
        ++rank;
        P4EST_ASSERT (rank < p4est->mpisize);
      }
      this->ghostCellMPIRank_[jl] = rank;
      this->cellGlobalIdx_[jl+nCells]
        = p8estPtr_->global_first_quadrant[rank] + q->p.piggy3.local_num;
    }

    this->points_.clear();
    this->faces_.clear();
    this->cells_.clear();

    this->volElems_.resize(nCells+nGhostCells);

    this->cece_.resize(nCells);
    this->cefa_.resize(nCells);

    for (int ii = 0; ii < nCells; ++ii){
      this->cece_[ii].clear();
      this->cefa_[ii].clear();
    }

    //Empty the boundary list information
    for(int ii = 0; ii < this->boundaryFaceIdxList_.size(); ++ii){
      this->boundaryFaceIdxList_[ii].clear();
    }

    p8estCece_.resize(nCells);

    //Update the connection information, face- & cell-
    //information and boundary Patchlist information
    //Before the iteration, *p8estPtr_ and *p8estGhostPtr_ data should be up-to-date
    //including the local cell and ghost cell information


    //int vector used to denote the old llabel of the cell
    std::vector<llabel>      cellTag(nCells+nGhostCells, -1);
    //int vector used to denote the old llabel of the face
    std::vector<llabel>      faceTag;  faceTag.clear();

    meshUpdateContex<T> ctx;
    ctx.meshPtr_ = this;
    ctx.mapPtr_  = &facePtrHashMap;
    ctx.cellTag  = &cellTag;
    ctx.faceTag  = &faceTag;

    //ONLY generate the faces information;
    //The volume information are also generated in volElems_
    //The oldfaceIdx information is stored in faceTag;
    //Meanwhile, the cece information are store in p8estCeCe_ structure
    p8est_iterate(p8estPtr_, p8estGhostPtr_, &ctx,
        updateMeshIterVol<T>, updateMeshIterFace<T>, NULL,  NULL);

    //Read the cece and cefa data from p8estCece_ structure
    for (int ii = 0; ii < this->nCells_; ++ii){
      for (int f = 0; f < P8EST_FACES; ++f){
        if(!p8estCece_[ii].ne[f].finer){
          llabel cidx = p8estCece_[ii].ne[f].index[0];
          llabel pidx = p8estCece_[ii].ne[f].pidx[0];
          this->cece_[ii].push_back(cidx);
          this->cefa_[ii].push_back(pidx);
        }else{
          for(int h = 0; h < P8EST_HALF; ++h){
            llabel cidx = p8estCece_[ii].ne[f].index[h];
            llabel pidx = p8estCece_[ii].ne[f].pidx[h];
            this->cece_[ii].push_back(cidx);
            this->cefa_[ii].push_back(pidx);
          }
        }
      }
    }

    //Update the send and receive pattern;
    this->sendPattern_   .resize(this->mpiSize());
    this->receivePattern_.resize(this->mpiSize());

    for (int ii = 0; ii < this->mpiSize(); ++ii){
      //Send
      llabel end    = p8estGhostPtr_->mirror_proc_offsets[ii+1];
      llabel start  = p8estGhostPtr_->mirror_proc_offsets[ii];

      this->sendPattern_[ii].clear();
      if(end-start)this->sendPattern_[ii].resize(end-start);
      for (int jj = 0; jj < end-start; ++jj){
        p4est_locidx_t mirrorIdx = p8estGhostPtr_->mirror_proc_mirrors[start+jj];
        p8est_quadrant_t *q  = p8est_quadrant_array_index(&(p8estGhostPtr_->mirrors),mirrorIdx);
        llabel localIdx = q->p.piggy3.local_num;
        this->sendPattern_[ii][jj] = localIdx;
      }
      //Receive
      end    = p8estGhostPtr_->proc_offsets[ii+1];
      start  = p8estGhostPtr_->proc_offsets[ii];
      this->receivePattern_[ii].clear();
      if(end-start) this->receivePattern_[ii].resize(end-start);
      for (int jj = 0; jj < end-start; ++jj){
        p4est_locidx_t ghostIdx = start + jj;
        p8est_quadrant_t *q  = p8est_quadrant_array_index(&(p8estGhostPtr_->ghosts),ghostIdx);
        llabel localIdx = q->p.piggy3.local_num;
        this->receivePattern_[ii][jj] = localIdx + p8estPtr_->global_first_quadrant[ii];
      }
    }
    //Initiate the ghost cell of the volElems_ elements;
    //!!!!!!!!!!!!To be completed later!!!!!!!!!!!!!!//
    ghostDataSynchronizeDeepCopy<T>(*this, 1, &volElems_);
    ///////////////////////////////////////////////////

    faceElems1_.resize(this->nFaces_);
    if(oversetFaceIntPt) faceElems2_.resize(this->nFaces_);

    int o1 = T::refFrame::order;
    int o2 = 3*o1/2;
    o1 = 2*o1+1;
    o2 = 2*o2+1;
    //Then initiate the faceElement_
    for (unsigned int nf = 0; nf < this->nFaces_; ++nf){
      if(faceTag[nf]!=-1){
        faceElems1_[nf] = oldfaceElems1_[faceTag[nf]];
        if(oversetFaceIntPt) faceElems2_[nf] = oldfaceElems2_[faceTag[nf]];
      }else{
        llabel cIdx1 = this->faces_[nf].lcellIdx.first;
        llabel cIdx2 = this->faces_[nf].lcellIdx.second;

        faceElems1_[nf].initiate(o1, this->faces_[nf], volElems_[cIdx1]);
        if(oversetFaceIntPt) faceElems2_[nf].initiate(
            o2, this->faces_[nf], volElems_[cIdx1]
            );
      }
    }

    //Update the boundaryFaceIdxList_
    //Firstly, clean all the boundaryFacIdxList_ information
    for (std::map<int, std::vector<int> >::iterator
        it = boundaryFaceIdxList_.begin();
        it != boundaryFaceIdxList_.end(); ++it){
      it->second.clear();
    }
    for (unsigned int nf = 0; nf < this->nFaces_; ++nf){
      int l = this->faces()[nf].patchIdx.first;
      int r = this->faces()[nf].patchIdx.second;
      if(!boundaryPropertyList_[l].isPeriodic){
        assert(l == r || r == -1);
        boundaryFaceIdxList_[l].push_back(nf);
      }
    }


    faceVector                           temp5;
    llabelVectorVector                   temp6;
    std::vector<faceElement>             tempFaceElem1;
    std::vector<faceElement>             tempFaceElem2;
    //std::vector<T>                       tempVolElem;

    oldfaceElems1_.swap(tempFaceElem1);
    oldfaceElems2_.swap(tempFaceElem2);
    oldfaces_.swap(temp5);
    oldcefa_.swap(temp6);
    //oldvolElems_.swap(tempVolElem);
    //Do not clear the oldvolElems_ data because it will be used in data mapping
  }

  //template<class T>
  //void dynamicMesh<T>::checkGeom()
  //{
    //int nCells = this->nCells();
    //for (int ii = 0; ii < nCells; ++ii){
      //int nfaces = this->cefa()[ii].size();
      //real3 sumnorm = real3(0,0,0);
      //for(int nf = 0; nf < nfaces; ++nf){
        //int fidx = this->cefa()[ii][nf];
        //const face& f = this->faces()[fidx];
        //int sign = -1;
        //if(f.firstCell==f.secondCell) continue;
        //if(f.firstCell==ii) sign = 1;
        //int npoints = this->faceGaussPairs_[fidx].size();
        //for(int np = 0; np < npoints; ++np){
          //sumnorm += sign*this->faceNormalVectors_[fidx][np]*this->faceGaussPairs_[fidx][np].Jaxw;
        //}
      //}
      //if(sumnorm.length() > 1e-3)
      //printf("Sum of norm of cell %d is %f %f %f\n",ii,sumnorm.x,sumnorm.y,sumnorm.z);
    //}
    //int nFaces = this->nFaces();
    //for (int ii = 0; ii < nFaces; ++ii){
      //const face& f = this->faces()[ii];
      //if(f.norm.length() < 0.5) printf("Fatal error!!\n");
    //}
  //}
}
#endif
