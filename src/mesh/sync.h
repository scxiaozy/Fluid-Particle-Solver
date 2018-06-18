#ifndef FPS_MESH_SYNC_H
#define FPS_MESH_SYNC_H
#include <vector>
#include <map>
#include <cassert>
#include <mpi.h>
#include <primitiveMesh.h>
#define PARTITION_USER_DATA_PACKAGE_SIZE 100000000
namespace fps
{
  template<class T>
  int ghostDataSynchronizeDeepCopy(
      const primitiveMesh &mesh, int n, std::vector<T> ** f);
  template<class T>
  int ghostDataSynchronizeDeepCopy(
      const primitiveMesh &mesh, int n, std::vector<T> * f);

  template<class T>
  int ghostDataSynchronizeDeepCopy(
      const primitiveMesh &mesh, int n, std::vector<T> ** f)
  {
    size_t serializedSize(const T&);
    void serialize(const T&, char*);
    void unserialize(T&, char*);
    size_t transfer_package_size = PARTITION_USER_DATA_PACKAGE_SIZE;
    int mpiSize = mesh.mpiSize();
    size_t *nbyteReceive  = new size_t[mpiSize];
    size_t *nbyteSend     = new size_t[mpiSize];
    size_t *nbyteReceive2  = new size_t[mpiSize];
    size_t *nbyteSend2     = new size_t[mpiSize];
    MPI_Request *reqSend   = new MPI_Request[mpiSize];
    MPI_Request *reqReceive  = new MPI_Request[mpiSize];
    MPI_Status *status = new MPI_Status[mpiSize];

    int error = 0;
    //The size of the vector should be not smaller than the mesh size
    for (int ii = 0; ii < n; ++ii){
      if((*f[ii]).size() < mesh.nCells()){
        std::cout << "Unmatched size between " << ii
          << "-th std::vector<T> and the mesh" << std::endl;
        error = 1;
      }
    }
    int nerr;
    MPI_Allreduce(&error, &nerr, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
    if(nerr){
      MPI_Abort(mesh.mpiComm(), -1);
    }else{
      for (int ii = 0; ii < n; ++ii){
        if((*f[ii]).size() < mesh.nCells()+mesh.nGhostCells()){
          (*f[ii]).resize((mesh.nCells()+mesh.nGhostCells()));
        }
      }
    }

    //Calculate the send size of this dynamic class
    size_t buffsize = 0;
    for (int ir = 0; ir < mesh.mpiSize(); ++ir){
      nbyteSend[ir] = 0;
      if(ir == mesh.mpiRank()) continue;
      for (unsigned int ii = 0; ii < mesh.sendPattern()[ir].size(); ++ii){
        int cIdx = mesh.sendPattern()[ir][ii];
        //The size of glabel indicating the global index of the data;
        nbyteSend[ir]+= sizeof(glabel)/sizeof(char);
        for (int jj = 0; jj < n; ++jj){
          //The size of the (unsigned int) indicating the data section's size;
          nbyteSend[ir]+=sizeof(size_t)/sizeof(char);
          //The size of the data section
          nbyteSend[ir]+=serializedSize((*f[jj])[cIdx]);
        }
      }
      buffsize += nbyteSend[ir] < transfer_package_size ? nbyteSend[ir]:transfer_package_size;
    }

    if(buffsize==0){
      MPI_Barrier(mesh.mpiComm());
      delete[] nbyteReceive;
      delete[] nbyteSend   ;
      delete[] nbyteReceive2;
      delete[] nbyteSend2   ;
      delete[] reqSend;
      delete[] reqReceive;
      delete[] status;
      return 0;
    }
    buffsize += mesh.mpiSize()*MPI_BSEND_OVERHEAD;

    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Isend((char*)(nbyteSend+ii), sizeof(size_t)/sizeof(char), MPI_CHAR, ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+ii);
    }
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Irecv((char*)(nbyteReceive+ii), sizeof(size_t)/sizeof(char), MPI_CHAR, ii, ii, mesh.mpiComm(), reqReceive+ii);
    }
    MPI_Waitall(mesh.mpiSize(), reqReceive, status);
    MPI_Waitall(mesh.mpiSize(), reqSend, status);

    char *buff = new char[buffsize];
    char **sendBuffer = new char*[mpiSize];
    char **receiveBuffer = new char*[mpiSize];
    for (int ii = 0; ii < mpiSize; ++ii){
      sendBuffer[ii] = nullptr;
      receiveBuffer[ii] = nullptr;
      nbyteSend2[ii] = nbyteSend[ii];
      nbyteReceive2[ii] = nbyteReceive[ii];
    }
    MPI_Buffer_attach(buff,buffsize);
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if((!nbyteReceive[ii])||(ii==mesh.mpiRank())) continue;
      receiveBuffer[ii] = new char[nbyteReceive[ii]];
      assert(receiveBuffer[ii]);
    }
    //Push data into send buffers;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(!nbyteSend[ii]||(ii==mesh.mpiRank())) continue;
      sendBuffer[ii] = new char[nbyteSend[ii]];
      assert(sendBuffer[ii]);
      char* current_pos = sendBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.sendPattern()[ii].size(); ++jj){
        llabel locidx = mesh.sendPattern()[ii][jj];
        //For every assignment, current_pos should move to the updated position;
        ((glabel*)current_pos)[0] = mesh.cellGlobalIdx()[locidx];
        current_pos += sizeof(glabel)/sizeof(char);
        for (int kk = 0; kk < n; ++kk){
          size_t ss = ((size_t*)current_pos)[0] =  serializedSize((*f[kk])[locidx]);
          current_pos += sizeof(size_t)/sizeof(char);
          serialize((*f[kk])[locidx], current_pos);
          current_pos += ss;
        }
      }
    }

    //To eliminate the MPI timeout problem, the data package size is limited to PARTITION_USER_DATA_PACKAGE_SIZE.
    int send_recv_count = 0;
    for(send_recv_count = 0; ; ++send_recv_count){
      int nreceiveRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteReceive2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Receive from rank " << ii << " : " << nbyteReceive2[ii] << std::endl;
        size_t receiveSize =
            nbyteReceive2[ii] < transfer_package_size?nbyteReceive2[ii]:transfer_package_size;
        MPI_Irecv(receiveBuffer[ii]+transfer_package_size*send_recv_count,
            receiveSize,
            MPI_CHAR, ii, ii, mesh.mpiComm(), reqReceive+nreceiveRankSize);
        ++nreceiveRankSize;
        nbyteReceive2[ii]-=receiveSize;
      }
      int nsendRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteSend2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Send to rank " << ii << " : " << nbyteSend2[ii] << std::endl;
        size_t sendSize =
          nbyteSend2[ii] < transfer_package_size?nbyteSend2[ii]:transfer_package_size;
        MPI_Ibsend(sendBuffer[ii]+transfer_package_size*send_recv_count,
            sendSize,
            MPI_CHAR, ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+nsendRankSize);
        ++nsendRankSize;
        nbyteSend2[ii]-=sendSize;
      }
      MPI_Waitall(nreceiveRankSize, reqReceive, status);
      MPI_Waitall(nsendRankSize,    reqSend, status);
      MPI_Barrier(mesh.mpiComm());

      //Cannot break for a single rank, must synchronized
      int active = 1;
      if(!nreceiveRankSize&&!nsendRankSize) active = 0;
      int active_total;
      MPI_Allreduce(&active, &active_total, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
      if(!active_total){
        break;
      }
    }




    std::map<glabel, char*> ghostDataMap;
    std::map<glabel, char*>::iterator it;
    //Insert all the ghost data into map
    for (int ii = 0; ii < mpiSize; ++ii){
      if(!nbyteReceive[ii]||ii==mesh.mpiRank()) continue;
      char* current_pos = receiveBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.receivePattern()[ii].size(); ++jj){
        glabel gloidx = ((glabel*)current_pos)[0];
        current_pos += sizeof(glabel)/sizeof(char);
        assert(ghostDataMap.find(gloidx)==ghostDataMap.end());
        ghostDataMap.insert(std::pair<glabel,char*>(gloidx, current_pos));
        for(int kk = 0; kk < n; ++kk){
          size_t ss = ((size_t*)current_pos)[0];
          current_pos += (sizeof(size_t)/sizeof(char)+ss);
        }
      }
    }
    //fclose(fp);
    int nCells = mesh.nCells();
    for (int ii = 0; ii < mesh.nGhostCells(); ++ii){
      glabel gloidx = mesh.cellGlobalIdx()[ii+nCells];
      char *current_pos = ghostDataMap[gloidx];
      for (int kk = 0; kk < n; ++kk){
        size_t ss = ((size_t*)current_pos)[0];
        current_pos += sizeof(size_t)/sizeof(char);
        unserialize((*f[kk])[ii+nCells], current_pos);
        current_pos += ss;
      }
    }
    for (int ii = 0; ii < mpiSize; ++ii){
      delete []sendBuffer[ii];
      delete []receiveBuffer[ii];
    }
    delete[] sendBuffer;
    delete[] receiveBuffer;
    delete[] nbyteReceive;
    delete[] nbyteSend;
    delete[] nbyteReceive2;
    delete[] nbyteSend2;
    delete[] reqSend;
    delete[] reqReceive;
    delete[] status;
    int detachsize;
    char *temp_ptr;
    MPI_Buffer_detach(&temp_ptr, &detachsize);
    delete[] buff;
    buff = nullptr;
    MPI_Barrier(mesh.mpiComm());
    return 0;
  }

  template<class T>
  int ghostDataSynchronizeDeepCopy(
      const primitiveMesh &mesh, int n, std::vector<T> * f)
  {
    size_t serializedSize(const T&);
    void serialize(const T&, char*);
    void unserialize(T&, char*);
    size_t transfer_package_size = PARTITION_USER_DATA_PACKAGE_SIZE;
    int mpiSize = mesh.mpiSize();
    size_t *nbyteReceive  = new size_t[mpiSize];
    size_t *nbyteSend     = new size_t[mpiSize];
    size_t *nbyteReceive2  = new size_t[mpiSize];
    size_t *nbyteSend2     = new size_t[mpiSize];
    MPI_Request *reqSend   = new MPI_Request[mpiSize];
    MPI_Request *reqReceive  = new MPI_Request[mpiSize];
    MPI_Status *status = new MPI_Status[mpiSize];

    int error = 0;
    //The size of the vector should be not smaller than the mesh size
    for (int ii = 0; ii < n; ++ii){
      if((f[ii]).size() < mesh.nCells()){
        std::cout << "Unmatched size between " << ii
          << "-th std::vector<T> and the mesh" << std::endl;
        error = 1;
      }
    }
    int nerr;
    MPI_Allreduce(&error, &nerr, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
    if(nerr){
      MPI_Abort(mesh.mpiComm(), -1);
    }else{
      for (int ii = 0; ii < n; ++ii){
        if(f[ii].size() < mesh.nCells()+mesh.nGhostCells()){
          f[ii].resize((mesh.nCells()+mesh.nGhostCells()));
        }
      }
    }

    //Calculate the send size of this dynamic class
    size_t buffsize = 0;
    for (int ir = 0; ir < mesh.mpiSize(); ++ir){
      nbyteSend[ir] = 0;
      if(ir == mesh.mpiRank()) continue;
      for (unsigned int ii = 0; ii < mesh.sendPattern()[ir].size(); ++ii){
        int cIdx = mesh.sendPattern()[ir][ii];
        //The size of glabel indicating the global index of the data;
        nbyteSend[ir]+= sizeof(glabel)/sizeof(char);
        for (int jj = 0; jj < n; ++jj){
          //The size of the (unsigned int) indicating the data section's size;
          nbyteSend[ir]+=sizeof(size_t)/sizeof(char);
          //The size of the data section
          nbyteSend[ir]+=serializedSize((f[jj])[cIdx]);
        }
      }
      buffsize += nbyteSend[ir] < transfer_package_size ? nbyteSend[ir]:transfer_package_size;
    }
    if(buffsize==0){
      MPI_Barrier(mesh.mpiComm());
      delete[] nbyteReceive;
      delete[] nbyteSend   ;
      delete[] nbyteReceive2;
      delete[] nbyteSend2   ;
      delete[] reqSend;
      delete[] reqReceive;
      delete[] status;
      return 0;
    }

    buffsize += mesh.mpiSize()*MPI_BSEND_OVERHEAD;

    int nreqSend = 0;
    int nreqReceive = 0;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Isend((char*)(nbyteSend+ii), sizeof(size_t)/sizeof(char), MPI_CHAR, ii,
          mesh.mpiRank(), mesh.mpiComm(), reqSend+(nreqSend++));
    }
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Irecv((char*)(nbyteReceive+ii), sizeof(size_t)/sizeof(char), MPI_CHAR, ii,
          ii, mesh.mpiComm(), reqReceive+(nreqReceive++));
    }
    MPI_Waitall(nreqReceive, reqReceive, status);
    MPI_Waitall(nreqSend, reqSend, status);

    char *buff = new char[buffsize];
    char **sendBuffer = new char*[mpiSize];
    char **receiveBuffer = new char*[mpiSize];
    for (int ii = 0; ii < mpiSize; ++ii){
      sendBuffer[ii] = nullptr;
      receiveBuffer[ii] = nullptr;
      nbyteSend2[ii] = nbyteSend[ii];
      nbyteReceive2[ii] = nbyteReceive[ii];
    }
    MPI_Buffer_attach(buff,buffsize);
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if((!nbyteReceive[ii])||(ii==mesh.mpiRank())) continue;
      receiveBuffer[ii] = new char[nbyteReceive[ii]];
      assert(receiveBuffer[ii]);
    }
    //Push data into send buffers;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(!nbyteSend[ii]||(ii==mesh.mpiRank())) continue;
      sendBuffer[ii] = new char[nbyteSend[ii]];
      assert(sendBuffer[ii]);
      char* current_pos = sendBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.sendPattern()[ii].size(); ++jj){
        llabel locidx = mesh.sendPattern()[ii][jj];
        //For every assignment, current_pos should move to the updated position;
        size_t ss = ((glabel*)current_pos)[0] = mesh.cellGlobalIdx()[locidx];
        current_pos += sizeof(glabel)/sizeof(char);
        for (int kk = 0; kk < n; ++kk){
          size_t ss = ((size_t*)current_pos)[0] =  serializedSize(f[kk][locidx]);
          current_pos += sizeof(size_t)/sizeof(char);
          serialize((f[kk])[locidx], current_pos);
          current_pos += ss;
        }
      }
    }


    //To eliminate the MPI timeout problem, the data package size is limited to PARTITION_USER_DATA_PACKAGE_SIZE.
    int send_recv_count = 0;
    for(send_recv_count = 0; ; ++send_recv_count){
      int nreceiveRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteReceive2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Receive from rank " << ii << " : " << nbyteReceive2[ii] << std::endl;
        size_t receiveSize =
            nbyteReceive2[ii] < transfer_package_size ? nbyteReceive2[ii] : transfer_package_size;
        MPI_Irecv(receiveBuffer[ii]+transfer_package_size*send_recv_count,
            receiveSize,
            MPI_CHAR, ii, ii, mesh.mpiComm(), reqReceive+nreceiveRankSize);
        ++nreceiveRankSize;
        nbyteReceive2[ii]-=receiveSize;
      }
      int nsendRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteSend2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Send to rank " << ii << " : " << nbyteSend2[ii] << std::endl;
        size_t sendSize =
          nbyteSend2[ii] < transfer_package_size?nbyteSend2[ii]:transfer_package_size;
        MPI_Ibsend(sendBuffer[ii]+transfer_package_size*send_recv_count,
            sendSize,
            MPI_CHAR, ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+nsendRankSize);
        ++nsendRankSize;
        nbyteSend2[ii]-=sendSize;
      }

      MPI_Waitall(nreceiveRankSize, reqReceive, status);
      MPI_Waitall(nsendRankSize,    reqSend, status);
      MPI_Barrier(mesh.mpiComm());

      //Cannot break for a single rank, must synchronized
      int active = 1;
      if(!nreceiveRankSize&&!nsendRankSize) active = 0;

      int active_total;
      MPI_Allreduce(&active, &active_total, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
      if(!active_total){
        break;
      }
    }


    std::map<glabel, char*> ghostDataMap;
    std::map<glabel, char*>::iterator it;
    //Insert all the ghost data into map
    for (int ii = 0; ii < mpiSize; ++ii){
      if(!nbyteReceive[ii]||ii==mesh.mpiRank()) continue;
      char* current_pos = receiveBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.receivePattern()[ii].size(); ++jj){
        glabel gloidx = ((glabel*)current_pos)[0];

        current_pos += sizeof(glabel)/sizeof(char);
        assert(ghostDataMap.find(gloidx)==ghostDataMap.end());
        ghostDataMap.insert(std::pair<glabel,char*>(gloidx, current_pos));
        for(int kk = 0; kk < n; ++kk){
          size_t ss = ((size_t*)current_pos)[0];
          current_pos += (sizeof(size_t)/sizeof(char)+ss);
        }
      }
    }
    llabel nCells = mesh.nCells();
    for (int ii = 0; ii < mesh.nGhostCells(); ++ii){
      glabel gloidx = mesh.cellGlobalIdx()[ii+nCells];
      char *current_pos = ghostDataMap[gloidx];
      for (int kk = 0; kk < n; ++kk){
        size_t ss = ((size_t*)current_pos)[0];
        current_pos += sizeof(size_t)/sizeof(char);
        unserialize((f[kk])[ii+nCells], current_pos);
        current_pos += ss;
      }
    }

    for (int ii = 0; ii < mpiSize; ++ii){
      delete []sendBuffer[ii];
      delete []receiveBuffer[ii];
    }
    delete[] sendBuffer;
    delete[] receiveBuffer;
    delete[] nbyteReceive;
    delete[] nbyteSend;
    delete[] nbyteReceive2;
    delete[] nbyteSend2;
    delete[] reqSend;
    delete[] reqReceive;
    delete[] status;
    int detachsize;
    char *temp_ptr;
    MPI_Buffer_detach(&temp_ptr, &detachsize);
    delete[] buff;
    buff = nullptr;
    MPI_Barrier(mesh.mpiComm());

    return 0;
  }

  template<class T>
  int ghostDataSynchronizeShallowCopy(
      const primitiveMesh &mesh, int n, std::vector<T> ** f, int dataPack = 1)
  {
    size_t transfer_package_size = PARTITION_USER_DATA_PACKAGE_SIZE;
    int mpiSize = mesh.mpiSize();
    size_t *nbyteReceive  = new size_t[mpiSize];
    size_t *nbyteSend     = new size_t[mpiSize];
    size_t *nbyteReceive2  = new size_t[mpiSize];
    size_t *nbyteSend2     = new size_t[mpiSize];
    MPI_Request *reqSend   = new MPI_Request[mpiSize];
    MPI_Request *reqReceive  = new MPI_Request[mpiSize];
    MPI_Status *status = new MPI_Status[mpiSize];

    int error = 0;
    //The size of the vector should be not smaller than the mesh size
    for (int ii = 0; ii < n; ++ii){
      if((*f[ii]).size() < mesh.nCells()*dataPack){
        std::cout << "Unmatched size between " << ii
          << "-th std::vector<T> and the mesh in " << __FILE__ << " : " << __LINE__ << std::endl;
        error = 1;
      }
    }
    int nerr;
    MPI_Allreduce(&error, &nerr, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
    if(nerr){
      MPI_Abort(mesh.mpiComm(), -1);
    }else{
      for (int ii = 0; ii < n; ++ii){
        if((*f[ii]).size() < dataPack*(mesh.nCells()+mesh.nGhostCells())){
          (*f[ii]).resize((mesh.nCells()+mesh.nGhostCells()));
        }
      }
    }


    //Calculate the send size of this dynamic class
    size_t buffsize = 0;
    for (int ir = 0; ir < mesh.mpiSize(); ++ir){
      nbyteSend[ir] = 0;
      if(ir == mesh.mpiRank()) continue;
      //The size of glabel indicating the global index of the data;
      nbyteSend[ir]+= mesh.sendPattern()[ir].size()*
        (sizeof(glabel) + sizeof(T)*n*dataPack)/sizeof(char);
      buffsize += nbyteSend[ir] < transfer_package_size ? nbyteSend[ir]:transfer_package_size;
    }
    if(buffsize==0){
      MPI_Barrier(mesh.mpiComm());
      delete[] nbyteReceive;
      delete[] nbyteSend   ;
      delete[] nbyteReceive2;
      delete[] nbyteSend2   ;
      delete[] reqSend;
      delete[] reqReceive;
      delete[] status;
      return 0;
    }
    buffsize += mesh.mpiSize()*MPI_BSEND_OVERHEAD;

    int nreqSend = 0, nreqReceive = 0;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Isend((char*)(nbyteSend+ii), sizeof(size_t)/sizeof(char), MPI_CHAR,
          ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+nreqSend++);
    }
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Irecv((char*)(nbyteReceive+ii), sizeof(size_t)/sizeof(char), MPI_CHAR,
          ii, ii, mesh.mpiComm(), reqReceive+nreqReceive++);
    }
    MPI_Waitall(nreqReceive, reqReceive, status);
    MPI_Waitall(nreqSend,    reqSend,    status);

    char *buff = new char[buffsize];
    char **sendBuffer = new char*[mpiSize];
    char **receiveBuffer = new char*[mpiSize];
    for (int ii = 0; ii < mpiSize; ++ii){
      sendBuffer[ii] = nullptr;
      receiveBuffer[ii] = nullptr;
      nbyteSend2[ii] = nbyteSend[ii];
      nbyteReceive2[ii] = nbyteReceive[ii];
    }
    MPI_Buffer_attach(buff,buffsize);
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if((!nbyteReceive[ii])||(ii==mesh.mpiRank())) continue;
      receiveBuffer[ii] = new char[nbyteReceive[ii]];
      assert(receiveBuffer[ii]);
    }
    //Push data into send buffers;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(!nbyteSend[ii]||(ii==mesh.mpiRank())) continue;
      sendBuffer[ii] = new char[nbyteSend[ii]];
      assert(sendBuffer[ii]);
      char* current_pos = sendBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.sendPattern()[ii].size(); ++jj){
        llabel locidx = mesh.sendPattern()[ii][jj];
        //For every assignment, current_pos should move to the updated position;
        ((glabel*)current_pos)[0] = mesh.cellGlobalIdx()[locidx];
        current_pos += sizeof(glabel)/sizeof(char);
        for (int kk = 0; kk < n; ++kk){
          for (int dd = 0; dd < dataPack; ++dd){
            ((T*)current_pos)[0] = (*f[kk])[locidx*dataPack+dd];
            current_pos += sizeof(T)/sizeof(char);
          }
        }
      }
    }

    //To eliminate the MPI timeout problem, the data package size is limited to PARTITION_USER_DATA_PACKAGE_SIZE.
    int send_recv_count = 0;
    for(send_recv_count = 0; ; ++send_recv_count){
      int nreceiveRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteReceive2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Receive from rank " << ii << " : " << nbyteReceive2[ii] << std::endl;
        size_t receiveSize =
            nbyteReceive2[ii] < transfer_package_size?nbyteReceive2[ii]:transfer_package_size;
        MPI_Irecv(receiveBuffer[ii]+transfer_package_size*send_recv_count,
            receiveSize,
            MPI_CHAR, ii, ii, mesh.mpiComm(), reqReceive+nreceiveRankSize);
        ++nreceiveRankSize;
        nbyteReceive2[ii]-=receiveSize;
      }
      int nsendRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteSend2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Send to rank " << ii << " : " << nbyteSend2[ii] << std::endl;
        size_t sendSize =
          nbyteSend2[ii] < transfer_package_size?nbyteSend2[ii]:transfer_package_size;
        MPI_Ibsend(sendBuffer[ii]+transfer_package_size*send_recv_count,
            sendSize,
            MPI_CHAR, ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+nsendRankSize);
        ++nsendRankSize;
        nbyteSend2[ii]-=sendSize;
      }
      MPI_Waitall(nreceiveRankSize, reqReceive, status);
      MPI_Waitall(nsendRankSize,    reqSend, status);
      MPI_Barrier(mesh.mpiComm());
      //Cannot break for a single rank, must synchronized
      int active = 1;
      if(!nreceiveRankSize&&!nsendRankSize) active = 0;
      int active_total;
      MPI_Allreduce(&active, &active_total, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
      if(!active_total){
        break;
      }
    }




    std::map<glabel, char*> ghostDataMap;
    std::map<glabel, char*>::iterator it;
    //Insert all the ghost data into map
    for (int ii = 0; ii < mpiSize; ++ii){
      if(!nbyteReceive[ii]||ii==mesh.mpiRank()) continue;
      char* current_pos = receiveBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.receivePattern()[ii].size(); ++jj){
        glabel gloidx = ((glabel*)current_pos)[0];
        current_pos += sizeof(glabel)/sizeof(char);
        //assert(ghostDataMap.find(gloidx)==ghostDataMap.end());
        ghostDataMap.insert(std::pair<glabel,char*>(gloidx, current_pos));
        current_pos += (sizeof(T)/sizeof(char))*n*dataPack;
      }
    }

    int nCells = mesh.nCells();
    for (int ii = 0; ii < mesh.nGhostCells(); ++ii){
      glabel gloidx = mesh.cellGlobalIdx()[ii+nCells];
      char *current_pos = ghostDataMap[gloidx];
      for (int kk = 0; kk < n; ++kk){
        for (int dd = 0; dd < dataPack; ++dd){
          (*f[kk])[(ii+nCells)*dataPack + dd] = ((T*)(current_pos))[0];
          current_pos += sizeof(T)/sizeof(char);
        }
      }
    }
    for (int ii = 0; ii < mpiSize; ++ii){
      delete []sendBuffer[ii];
      delete []receiveBuffer[ii];
    }
    delete[] sendBuffer;
    delete[] receiveBuffer;
    delete[] nbyteReceive;
    delete[] nbyteSend;
    delete[] nbyteReceive2;
    delete[] nbyteSend2;
    delete[] reqSend;
    delete[] reqReceive;
    delete[] status;
    int detachsize;
    char *temp_ptr;
    MPI_Buffer_detach(&temp_ptr, &detachsize);
    delete[] buff;
    MPI_Barrier(mesh.mpiComm());
    return 0;
  }

  template<class T>
  int ghostDataSynchronizeShallowCopy(
      const primitiveMesh &mesh, int n, std::vector<T> * f, int dataPack = 1)
  {
    size_t transfer_package_size = PARTITION_USER_DATA_PACKAGE_SIZE;
    int mpiSize = mesh.mpiSize();
    size_t *nbyteReceive  = new size_t[mpiSize];
    size_t *nbyteSend     = new size_t[mpiSize];
    size_t *nbyteReceive2  = new size_t[mpiSize];
    size_t *nbyteSend2     = new size_t[mpiSize];
    MPI_Request *reqSend   = new MPI_Request[mpiSize];
    MPI_Request *reqReceive  = new MPI_Request[mpiSize];
    MPI_Status *status = new MPI_Status[mpiSize];

    int error = 0;
    //The size of the vector should be not smaller than the mesh size
    for (int ii = 0; ii < n; ++ii){
      if((f[ii]).size() < mesh.nCells()*dataPack){
        std::cout << "Unmatched size between " << ii
          << "-th std::vector<T> and the mesh in " << __FILE__ << " : " << __LINE__ << std::endl;
        error = 1;
      }
    }
    int nerr;
    MPI_Allreduce(&error, &nerr, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
    if(nerr){
      MPI_Abort(mesh.mpiComm(), -1);
    }else{
      for (int ii = 0; ii < n; ++ii){
        if(f[ii].size() < dataPack*(mesh.nCells()+mesh.nGhostCells())){
          f[ii].resize((mesh.nCells()+mesh.nGhostCells()));
        }
      }
    }


    //Calculate the send size of this dynamic class
    size_t buffsize = 0;
    for (int ir = 0; ir < mesh.mpiSize(); ++ir){
      nbyteSend[ir] = 0;
      if(ir == mesh.mpiRank()) continue;
      //The size of glabel indicating the global index of the data;
      nbyteSend[ir]+= mesh.sendPattern()[ir].size()*
        (sizeof(glabel) + sizeof(T)*n*dataPack)/sizeof(char);
      buffsize += nbyteSend[ir] < transfer_package_size ? nbyteSend[ir]:transfer_package_size;
    }
    if(buffsize==0){
      MPI_Barrier(mesh.mpiComm());
      delete[] nbyteReceive;
      delete[] nbyteSend   ;
      delete[] nbyteReceive2;
      delete[] nbyteSend2   ;
      delete[] reqSend;
      delete[] reqReceive;
      delete[] status;
      return 0;
    }
    buffsize += mesh.mpiSize()*MPI_BSEND_OVERHEAD;

    int nreqSend = 0, nreqReceive = 0;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Isend((char*)(nbyteSend+ii), sizeof(size_t)/sizeof(char), MPI_CHAR,
          ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+nreqSend++);
    }
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(ii == mesh.mpiRank()) continue;
      MPI_Irecv((char*)(nbyteReceive+ii), sizeof(size_t)/sizeof(char), MPI_CHAR,
          ii, ii, mesh.mpiComm(), reqReceive+nreqReceive++);
    }
    MPI_Waitall(nreqReceive, reqReceive, status);
    MPI_Waitall(nreqSend,    reqSend,    status);

    char *buff = new char[buffsize];
    char **sendBuffer = new char*[mpiSize];
    char **receiveBuffer = new char*[mpiSize];
    for (int ii = 0; ii < mpiSize; ++ii){
      sendBuffer[ii] = nullptr;
      receiveBuffer[ii] = nullptr;
      nbyteSend2[ii] = nbyteSend[ii];
      nbyteReceive2[ii] = nbyteReceive[ii];
    }
    MPI_Buffer_attach(buff,buffsize);
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if((!nbyteReceive[ii])||(ii==mesh.mpiRank())) continue;
      receiveBuffer[ii] = new char[nbyteReceive[ii]];
      assert(receiveBuffer[ii]);
    }
    //Push data into send buffers;
    for (int ii = 0; ii < mesh.mpiSize(); ++ii){
      if(!nbyteSend[ii]||(ii==mesh.mpiRank())) continue;
      sendBuffer[ii] = new char[nbyteSend[ii]];
      assert(sendBuffer[ii]);
      char* current_pos = sendBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.sendPattern()[ii].size(); ++jj){
        llabel locidx = mesh.sendPattern()[ii][jj];
        //For every assignment, current_pos should move to the updated position;
        ((glabel*)current_pos)[0] = mesh.cellGlobalIdx()[locidx];
        current_pos += sizeof(glabel)/sizeof(char);
        for (int kk = 0; kk < n; ++kk){
          for (int dd = 0; dd < dataPack; ++dd){
            ((T*)current_pos)[0] = (f[kk])[locidx*dataPack+dd];
            current_pos += sizeof(T)/sizeof(char);
          }
        }
      }
    }

    //To eliminate the MPI timeout problem, the data package size is limited to PARTITION_USER_DATA_PACKAGE_SIZE.
    int send_recv_count = 0;
    for(send_recv_count = 0; ; ++send_recv_count){
      int nreceiveRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteReceive2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Receive from rank " << ii << " : " << nbyteReceive2[ii] << std::endl;
        size_t receiveSize =
            nbyteReceive2[ii] < transfer_package_size?nbyteReceive2[ii]:transfer_package_size;
        MPI_Irecv(receiveBuffer[ii]+transfer_package_size*send_recv_count,
            receiveSize,
            MPI_CHAR, ii, ii, mesh.mpiComm(), reqReceive+nreceiveRankSize);
        ++nreceiveRankSize;
        nbyteReceive2[ii]-=receiveSize;
      }
      int nsendRankSize = 0;
      for (int ii = 0; ii < mpiSize; ++ii){
        if((!nbyteSend2[ii])||(ii==mesh.mpiRank())) continue;
        //test << send_recv_count << "  Send to rank " << ii << " : " << nbyteSend2[ii] << std::endl;
        size_t sendSize =
          nbyteSend2[ii] < transfer_package_size?nbyteSend2[ii]:transfer_package_size;
        MPI_Ibsend(sendBuffer[ii]+transfer_package_size*send_recv_count,
            sendSize,
            MPI_CHAR, ii, mesh.mpiRank(), mesh.mpiComm(), reqSend+nsendRankSize);
        ++nsendRankSize;
        nbyteSend2[ii]-=sendSize;
      }
      MPI_Waitall(nreceiveRankSize, reqReceive, status);
      MPI_Waitall(nsendRankSize,    reqSend, status);
      MPI_Barrier(mesh.mpiComm());
      //Cannot break for a single rank, must synchronized
      int active = 1;
      if(!nreceiveRankSize&&!nsendRankSize) active = 0;
      int active_total;
      MPI_Allreduce(&active, &active_total, 1, MPI_INT, MPI_SUM, mesh.mpiComm());
      if(!active_total){
        break;
      }
    }




    std::map<glabel, char*> ghostDataMap;
    std::map<glabel, char*>::iterator it;
    //Insert all the ghost data into map
    for (int ii = 0; ii < mpiSize; ++ii){
      if(!nbyteReceive[ii]||ii==mesh.mpiRank()) continue;
      char* current_pos = receiveBuffer[ii];
      for (unsigned int jj = 0; jj < mesh.receivePattern()[ii].size(); ++jj){
        glabel gloidx = ((glabel*)current_pos)[0];
        current_pos += sizeof(glabel)/sizeof(char);
        //assert(ghostDataMap.find(gloidx)==ghostDataMap.end());
        ghostDataMap.insert(std::pair<glabel,char*>(gloidx, current_pos));
        current_pos += (sizeof(T)/sizeof(char))*n*dataPack;
      }
    }

    int nCells = mesh.nCells();
    for (int ii = 0; ii < mesh.nGhostCells(); ++ii){
      glabel gloidx = mesh.cellGlobalIdx()[ii+nCells];
      char *current_pos = ghostDataMap[gloidx];
      for (int kk = 0; kk < n; ++kk){
        for (int dd = 0; dd < dataPack; ++dd){
          (f[kk])[(ii+nCells)*dataPack + dd] = ((T*)(current_pos))[0];
          current_pos += sizeof(T)/sizeof(char);
        }
      }
    }
    for (int ii = 0; ii < mpiSize; ++ii){
      delete []sendBuffer[ii];
      delete []receiveBuffer[ii];
    }
    delete[] sendBuffer;
    delete[] receiveBuffer;
    delete[] nbyteReceive;
    delete[] nbyteSend;
    delete[] nbyteReceive2;
    delete[] nbyteSend2;
    delete[] reqSend;
    delete[] reqReceive;
    delete[] status;
    int detachsize;
    char *temp_ptr;
    MPI_Buffer_detach(&temp_ptr, &detachsize);
    delete[] buff;
    MPI_Barrier(mesh.mpiComm());
    return 0;
  }
}
#endif
