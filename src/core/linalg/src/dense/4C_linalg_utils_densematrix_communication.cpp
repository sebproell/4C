// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_utils_densematrix_communication.hpp"

#include "4C_utils_exceptions.hpp"

#include <numeric>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Core::LinAlg::find_my_pos(int nummyelements, MPI_Comm comm)
{
  const int myrank = Core::Communication::my_mpi_rank(comm);
  const int numproc = Core::Communication::num_mpi_ranks(comm);

  std::vector<int> snum(numproc, 0);
  std::vector<int> rnum(numproc);
  snum[myrank] = nummyelements;

  Core::Communication::sum_all(snum.data(), rnum.data(), numproc, comm);

  return std::accumulate(rnum.data(), rnum.data() + myrank, 0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::allreduce_vector(
    const std::vector<int>& src, std::vector<int>& dest, MPI_Comm comm)
{
  // communicate size
  int localsize = static_cast<int>(src.size());
  int globalsize;
  Core::Communication::sum_all(&localsize, &globalsize, 1, comm);

  // communicate values
  int pos = find_my_pos(localsize, comm);
  std::vector<int> sendglobal(globalsize, 0);
  dest.resize(globalsize);
  std::copy(src.begin(), src.end(), sendglobal.data() + pos);
  Core::Communication::sum_all(sendglobal.data(), dest.data(), globalsize, comm);

  // sort & unique
  std::sort(dest.begin(), dest.end());
  std::vector<int>::iterator i = std::unique(dest.begin(), dest.end());
  dest.erase(i, dest.end());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::allreduce_e_map(std::vector<int>& rredundant, const Epetra_Map& emap)
{
  const int mynodepos =
      find_my_pos(emap.NumMyElements(), Core::Communication::unpack_epetra_comm(emap.Comm()));

  std::vector<int> sredundant(emap.NumGlobalElements(), 0);

  int* gids = emap.MyGlobalElements();
  std::copy(gids, gids + emap.NumMyElements(), sredundant.data() + mynodepos);

  rredundant.resize(emap.NumGlobalElements());
  Core::Communication::sum_all(sredundant.data(), rredundant.data(), emap.NumGlobalElements(),
      Core::Communication::unpack_epetra_comm(emap.Comm()));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::allreduce_e_map(std::map<int, int>& idxmap, const Epetra_Map& emap)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not emap.UniqueGIDs()) FOUR_C_THROW("works only for unique Epetra_Maps");
#endif

  idxmap.clear();

  std::vector<int> rredundant;
  allreduce_e_map(rredundant, emap);

  for (std::size_t i = 0; i < rredundant.size(); ++i)
  {
    idxmap[rredundant[i]] = i;
  }
}

/*----------------------------------------------------------------------*
 |  create an allreduced map on a distinct processor (public)  gjb 12/07|
 *----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::allreduce_e_map(const Epetra_Map& emap, const int pid)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not emap.UniqueGIDs()) FOUR_C_THROW("works only for unique Epetra_Maps");
#endif
  std::vector<int> rv;
  allreduce_e_map(rv, emap);
  std::shared_ptr<Epetra_Map> rmap;

  if (Core::Communication::my_mpi_rank(Core::Communication::unpack_epetra_comm(emap.Comm())) == pid)
  {
    rmap = std::make_shared<Epetra_Map>(-1, rv.size(), rv.data(), 0, emap.Comm());
    // check the map
    FOUR_C_ASSERT(rmap->NumMyElements() == rmap->NumGlobalElements(),
        "Processor with pid does not get all map elements");
  }
  else
  {
    rv.clear();
    rmap = std::make_shared<Epetra_Map>(-1, 0, nullptr, 0, emap.Comm());
    // check the map
    FOUR_C_ASSERT(rmap->NumMyElements() == 0, "At least one proc will keep a map element");
  }
  return rmap;
}

/*----------------------------------------------------------------------*
 |  create an allreduced map on EVERY processor (public)        tk 12/07|
 *----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::allreduce_e_map(const Epetra_Map& emap)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not emap.UniqueGIDs()) FOUR_C_THROW("works only for unique Epetra_Maps");
#endif
  std::vector<int> rv;
  allreduce_e_map(rv, emap);
  std::shared_ptr<Epetra_Map> rmap;

  rmap = std::make_shared<Epetra_Map>(-1, rv.size(), rv.data(), 0, emap.Comm());

  return rmap;
}

/*----------------------------------------------------------------------*
|  create an allreduced map on EVERY processor (public)                 |
 *----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::allreduce_overlapping_e_map(const Epetra_Map& emap)
{
  std::vector<int> rv;
  allreduce_e_map(rv, emap);

  // remove duplicates
  std::set<int> rs(rv.begin(), rv.end());
  rv.assign(rs.begin(), rs.end());

  return std::make_shared<Epetra_Map>(-1, rv.size(), rv.data(), 0, emap.Comm());
}

/*----------------------------------------------------------------------*
| create an allreduced map on a distinct processor (public)  ghamm 10/14|
 *----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::allreduce_overlapping_e_map(
    const Epetra_Map& emap, const int pid)
{
  std::vector<int> rv;
  allreduce_e_map(rv, emap);
  std::shared_ptr<Epetra_Map> rmap;

  if (Core::Communication::my_mpi_rank(Core::Communication::unpack_epetra_comm(emap.Comm())) == pid)
  {
    // remove duplicates only on proc pid
    std::set<int> rs(rv.begin(), rv.end());
    rv.assign(rs.begin(), rs.end());

    rmap = std::make_shared<Epetra_Map>(-1, rv.size(), rv.data(), 0, emap.Comm());
    // check the map
    FOUR_C_ASSERT(rmap->NumMyElements() == rmap->NumGlobalElements(),
        "Processor with pid does not get all map elements");
  }
  else
  {
    rv.clear();
    rmap = std::make_shared<Epetra_Map>(-1, 0, nullptr, 0, emap.Comm());
    // check the map
    FOUR_C_ASSERT(rmap->NumMyElements() == 0, "At least one proc will keep a map element");
  }
  return rmap;
}

/*----------------------------------------------------------------------*
 |  Send and receive lists of ints.  (heiner 09/07)                     |
 *----------------------------------------------------------------------*/
void Core::LinAlg::all_to_all_communication(
    MPI_Comm comm, const std::vector<std::vector<int>>& send, std::vector<std::vector<int>>& recv)
{
  if (Core::Communication::num_mpi_ranks(comm) == 1)
  {
    FOUR_C_ASSERT(send.size() == 1, "there has to be just one entry for sending");

    // make a copy
    recv.clear();
    recv.push_back(send[0]);
  }
  else
  {
    std::vector<int> sendbuf;
    std::vector<int> sendcounts;
    sendcounts.reserve(Core::Communication::num_mpi_ranks(comm));
    std::vector<int> sdispls;
    sdispls.reserve(Core::Communication::num_mpi_ranks(comm));

    int displacement = 0;
    sdispls.push_back(0);
    for (std::vector<std::vector<int>>::const_iterator iter = send.begin(); iter != send.end();
        ++iter)
    {
      sendbuf.insert(sendbuf.end(), iter->begin(), iter->end());
      sendcounts.push_back(iter->size());
      displacement += iter->size();
      sdispls.push_back(displacement);
    }

    std::vector<int> recvcounts(Core::Communication::num_mpi_ranks(comm));

    // initial communication: Request. Send and receive the number of
    // ints we communicate with each process.

    int status = MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    if (status != MPI_SUCCESS) FOUR_C_THROW("MPI_Alltoall returned status={}", status);

    std::vector<int> rdispls;
    rdispls.reserve(Core::Communication::num_mpi_ranks(comm));

    displacement = 0;
    rdispls.push_back(0);
    for (std::vector<int>::const_iterator iter = recvcounts.begin(); iter != recvcounts.end();
        ++iter)
    {
      displacement += *iter;
      rdispls.push_back(displacement);
    }

    std::vector<int> recvbuf(rdispls.back());

    // transmit communication: Send and get the data.

    status = MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT,
        recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_INT, comm);
    if (status != MPI_SUCCESS) FOUR_C_THROW("MPI_Alltoallv returned status={}", status);

    recv.clear();
    for (int proc = 0; proc < Core::Communication::num_mpi_ranks(comm); ++proc)
    {
      recv.push_back(
          std::vector<int>(recvbuf.data() + rdispls[proc], recvbuf.data() + rdispls[proc + 1]));
    }
  }
}

/*----------------------------------------------------------------------*
 |  Send and receive lists of ints.                                     |
 *----------------------------------------------------------------------*/
void Core::LinAlg::all_to_all_communication(
    MPI_Comm comm, const std::vector<std::vector<int>>& send, std::vector<int>& recv)
{
  if (Core::Communication::num_mpi_ranks(comm) == 1)
  {
    FOUR_C_ASSERT(send.size() == 1, "there has to be just one entry for sending");

    // make a copy
    recv.clear();
    recv = send[0];
  }
  else
  {
    std::vector<int> sendbuf;
    std::vector<int> sendcounts;
    sendcounts.reserve(Core::Communication::num_mpi_ranks(comm));
    std::vector<int> sdispls;
    sdispls.reserve(Core::Communication::num_mpi_ranks(comm));

    int displacement = 0;
    sdispls.push_back(0);
    for (std::vector<std::vector<int>>::const_iterator iter = send.begin(); iter != send.end();
        ++iter)
    {
      sendbuf.insert(sendbuf.end(), iter->begin(), iter->end());
      sendcounts.push_back(iter->size());
      displacement += iter->size();
      sdispls.push_back(displacement);
    }

    std::vector<int> recvcounts(Core::Communication::num_mpi_ranks(comm));

    // initial communication: Request. Send and receive the number of
    // ints we communicate with each process.

    int status = MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    if (status != MPI_SUCCESS) FOUR_C_THROW("MPI_Alltoall returned status={}", status);

    std::vector<int> rdispls;
    rdispls.reserve(Core::Communication::num_mpi_ranks(comm));

    displacement = 0;
    rdispls.push_back(0);
    for (std::vector<int>::const_iterator iter = recvcounts.begin(); iter != recvcounts.end();
        ++iter)
    {
      displacement += *iter;
      rdispls.push_back(displacement);
    }

    std::vector<int> recvbuf(rdispls.back());

    // transmit communication: Send and get the data.

    recv.clear();
    recv.resize(rdispls.back());

    status = MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, recv.data(),
        recvcounts.data(), rdispls.data(), MPI_INT, comm);
    if (status != MPI_SUCCESS) FOUR_C_THROW("MPI_Alltoallv returned status={}", status);
  }
}

FOUR_C_NAMESPACE_CLOSE
