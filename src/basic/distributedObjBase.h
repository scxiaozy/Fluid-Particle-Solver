#ifndef FPS_BASIC_DISTRIBUTED_OBJ_BASE
#define FPS_BASIC_DISTRIBUTED_OBJ_BASE
#include <config.inc>

#include <mpi.h>
namespace fps
{
    class distributedObjBase
    {
        protected:
            int mpiRank_;
            int mpiSize_;
            MPI_Comm mpiComm_;
        protected:
            distributedObjBase(MPI_Comm mpiComm = MPI_COMM_WORLD):mpiComm_(mpiComm)
            {
                MPI_Comm_rank(mpiComm_, &mpiRank_);
                MPI_Comm_size(mpiComm_, &mpiSize_);
            }
        public:
            int mpiRank()const {return mpiRank_;}
            int mpiSize()const {return mpiSize_;}
            MPI_Comm mpiComm()const {return mpiComm_;}
        public:
            inline bool inSameCommDomain(const distributedObjBase& a,
                    const distributedObjBase& b);
    };

    //Definition of inline function
    inline bool distributedObjBase::inSameCommDomain(const distributedObjBase& a,
            const distributedObjBase& b)
    {
        return a.mpiComm_ == b.mpiComm_;
    }
}
#endif
