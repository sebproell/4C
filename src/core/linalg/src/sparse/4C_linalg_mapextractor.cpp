// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_mapextractor.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"

#include <Teuchos_getConst.hpp>

#include <cmath>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MultiMapExtractor::MultiMapExtractor() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MultiMapExtractor::MultiMapExtractor(
    const Epetra_Map& fullmap, const std::vector<std::shared_ptr<const Epetra_Map>>& maps)
{
  setup(fullmap, maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::setup(
    const Epetra_Map& fullmap, const std::vector<std::shared_ptr<const Epetra_Map>>& maps)
{
  fullmap_ = std::make_shared<Epetra_Map>(fullmap);
  maps_ = maps;

  importer_.resize(maps_.size());
  for (unsigned i = 0; i < importer_.size(); ++i)
  {
    if (maps_[i] != nullptr)
    {
      importer_[i] = std::make_shared<Epetra_Import>(*maps_[i], *fullmap_);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::check_for_valid_map_extractor() const
{
  if (maps_.size() == 0)
  {
    FOUR_C_THROW("no maps_ available");
  }

  for (unsigned i = 0; i < maps_.size(); ++i)
  {
    if (maps_[i] != nullptr)
    {
      if (maps_[i]->DataPtr() == nullptr)
      {
        FOUR_C_THROW("Got zero data pointer on setup of block {} of maps_\n", i);
      }
      if (not maps_[i]->UniqueGIDs())
      {
        FOUR_C_THROW("map {} not unique", i);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::MultiMapExtractor::merge_maps(
    const std::vector<std::shared_ptr<const Epetra_Map>>& maps)
{
  if (maps.size() == 0) FOUR_C_THROW("no maps to merge");
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == nullptr) FOUR_C_THROW("can not merge extractor with null maps");
    if (not maps[i]->UniqueGIDs()) FOUR_C_THROW("map {} not unique", i);
  }
  std::set<int> mapentries;
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    std::copy(map.MyGlobalElements(), map.MyGlobalElements() + map.NumMyElements(),
        std::inserter(mapentries, mapentries.begin()));
  }
  return Core::LinAlg::create_map(
      mapentries, Core::Communication::unpack_epetra_comm(maps[0]->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::MultiMapExtractor::merge_maps_keep_order(
    const std::vector<std::shared_ptr<const Epetra_Map>>& maps)
{
  if (maps.empty()) FOUR_C_THROW("no maps to merge");

  // sanity checks
  for (std::size_t i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == nullptr) FOUR_C_THROW("can not merge extractor with null maps");
    if (not maps[i]->UniqueGIDs()) FOUR_C_THROW("map {} not unique", i);
  }

  // collect gids
  std::vector<int> gids;
  for (std::size_t i = 0; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    for (int j = 0; j < map.NumMyElements(); ++j) gids.push_back(map.GID(j));
  }

  // build combined map
  std::shared_ptr<Epetra_Map> fullmap =
      std::make_shared<Epetra_Map>(-1, gids.size(), gids.data(), 0, maps[0]->Comm());
  return fullmap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Core::LinAlg::MultiMapExtractor::intersect_maps(
    const std::vector<std::shared_ptr<const Epetra_Map>>& maps)
{
  if (maps.size() == 0) FOUR_C_THROW("no maps to intersect");
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == nullptr) FOUR_C_THROW("can not intersect extractor with null maps");
    if (not maps[i]->UniqueGIDs()) FOUR_C_THROW("map {} not unique", i);
  }
  std::set<int> mapentries(
      maps[0]->MyGlobalElements(), maps[0]->MyGlobalElements() + maps[0]->NumMyElements());
  for (unsigned i = 1; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    std::set<int> newset;
    int numele = map.NumMyElements();
    int* ele = map.MyGlobalElements();
    for (int j = 0; j < numele; ++j)
    {
      if (mapentries.find(ele[j]) != mapentries.end())
      {
        newset.insert(ele[j]);
      }
    }
    std::swap(mapentries, newset);
  }
  return Core::LinAlg::create_map(
      mapentries, Core::Communication::unpack_epetra_comm(maps[0]->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Core::LinAlg::MultiMapExtractor::extract_vector(
    const Core::LinAlg::Vector<double>& full, int block) const
{
  if (maps_[block] == nullptr) FOUR_C_THROW("null map at block {}", block);
  std::shared_ptr<Core::LinAlg::Vector<double>> vec =
      std::make_shared<Core::LinAlg::Vector<double>>(*maps_[block]);
  extract_vector(full, block, *vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Core::LinAlg::MultiMapExtractor::extract_vector(
    const Core::LinAlg::MultiVector<double>& full, int block) const
{
  if (maps_[block] == nullptr) FOUR_C_THROW("null map at block {}", block);
  std::shared_ptr<Core::LinAlg::MultiVector<double>> vec =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*maps_[block], full.NumVectors());
  extract_vector(full, block, *vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::extract_vector(const Core::LinAlg::MultiVector<double>& full,
    int block, Core::LinAlg::MultiVector<double>& partial) const
{
  if (maps_[block] == nullptr) FOUR_C_THROW("null map at block {}", block);
  int err = partial.Import(full, *importer_[block], Insert);
  if (err) FOUR_C_THROW("Import using importer returned err={}", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Core::LinAlg::MultiMapExtractor::insert_vector(
    const Core::LinAlg::Vector<double>& partial, int block) const
{
  std::shared_ptr<Core::LinAlg::Vector<double>> full =
      std::make_shared<Core::LinAlg::Vector<double>>(*fullmap_);
  insert_vector(partial, block, *full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Core::LinAlg::MultiMapExtractor::insert_vector(
    const Core::LinAlg::MultiVector<double>& partial, int block) const
{
  std::shared_ptr<Core::LinAlg::MultiVector<double>> full =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*fullmap_, partial.NumVectors());
  insert_vector(partial, block, *full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::insert_vector(
    const Core::LinAlg::MultiVector<double>& partial, int block,
    Core::LinAlg::MultiVector<double>& full) const
{
  if (maps_[block] == nullptr) FOUR_C_THROW("null map at block {}", block);
  int err = full.Export(partial, *importer_[block], Insert);
  if (err) FOUR_C_THROW("Export using importer returned err={}", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::add_vector(const Core::LinAlg::MultiVector<double>& partial,
    int block, Core::LinAlg::MultiVector<double>& full, double scale) const
{
  std::shared_ptr<Core::LinAlg::MultiVector<double>> v = extract_vector(full, block);
  if (not v->Map().SameAs(partial.Map())) FOUR_C_THROW("The maps of the vectors must be the same!");
  v->Update(scale, partial, 1.0);
  insert_vector(*v, block, full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::put_scalar(
    Core::LinAlg::Vector<double>& full, int block, double scalar) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *full_map();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) FOUR_C_THROW("maps do not match");
    full[lid] = scalar;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Core::LinAlg::MultiMapExtractor::norm2(
    const Core::LinAlg::Vector<double>& full, int block) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *full_map();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  double local_norm = 0;

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) FOUR_C_THROW("maps do not match");
    double value = full[lid];
    local_norm += value * value;
  }

  double global_norm = 0;
  Core::Communication::sum_all(
      &local_norm, &global_norm, 1, Core::Communication::unpack_epetra_comm(fm.Comm()));
  return std::sqrt(global_norm);
}


/*----------------------------------------------------------------------*
 | Scale one block only                                      fang 08/16 |
 *----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::scale(
    Core::LinAlg::Vector<double>& full, int block, double scalar) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *full_map();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) FOUR_C_THROW("maps do not match");
    full[lid] *= scalar;
  }
}


/*----------------------------------------------------------------------*
 | Scale one block only                                      fang 08/16 |
 *----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::scale(
    Core::LinAlg::MultiVector<double>& full, int block, double scalar) const
{
  for (int i = 0; i < full.NumVectors(); ++i) scale(full(i), block, scalar);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MapExtractor::MapExtractor() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MapExtractor::MapExtractor(const Epetra_Map& fullmap,
    std::shared_ptr<const Epetra_Map> condmap, std::shared_ptr<const Epetra_Map> othermap)
{
  setup(fullmap, condmap, othermap);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MapExtractor::MapExtractor(
    const Epetra_Map& fullmap, std::shared_ptr<const Epetra_Map> partialmap, bool iscondmap)
{
  // initialise other DOFs by inserting all DOFs of full map
  std::set<int> othergids;
  const int* fullgids = fullmap.MyGlobalElements();
  copy(fullgids, fullgids + fullmap.NumMyElements(), inserter(othergids, othergids.begin()));

  // throw away all DOFs which are in condmap
  if (partialmap->NumMyElements() > 0)
  {
    const int* condgids = partialmap->MyGlobalElements();
    for (int lid = 0; lid < partialmap->NumMyElements(); ++lid) othergids.erase(condgids[lid]);
  }

  // create (non-overlapping) othermap for non-condmap DOFs
  std::shared_ptr<Epetra_Map> othermap =
      Core::LinAlg::create_map(othergids, Core::Communication::unpack_epetra_comm(fullmap.Comm()));

  // create the extractor based on choice 'iscondmap'
  if (iscondmap)
    setup(fullmap, partialmap, othermap);
  else
    setup(fullmap, othermap, partialmap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MapExtractor::setup(const Epetra_Map& fullmap,
    const std::shared_ptr<const Epetra_Map>& condmap,
    const std::shared_ptr<const Epetra_Map>& othermap)
{
  std::vector<std::shared_ptr<const Epetra_Map>> maps;
  maps.push_back(othermap);
  maps.push_back(condmap);
  MultiMapExtractor::setup(fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MapExtractor::setup(
    const Epetra_Map& fullmap, const std::shared_ptr<const Epetra_Map>& partialmap, bool iscondmap)
{
  // initialise other DOFs by inserting all DOFs of full map
  std::set<int> othergids;
  const int* fullgids = fullmap.MyGlobalElements();
  copy(fullgids, fullgids + fullmap.NumMyElements(), inserter(othergids, othergids.begin()));

  // throw away all DOFs which are in condmap
  if (partialmap->NumMyElements() > 0)
  {
    const int* condgids = partialmap->MyGlobalElements();
    for (int lid = 0; lid < partialmap->NumMyElements(); ++lid) othergids.erase(condgids[lid]);
  }

  // create (non-overlapping) othermap for non-condmap DOFs
  std::shared_ptr<Epetra_Map> othermap =
      Core::LinAlg::create_map(othergids, Core::Communication::unpack_epetra_comm(fullmap.Comm()));

  // create the extractor based on choice 'iscondmap'
  if (iscondmap)
    setup(fullmap, partialmap, othermap);
  else
    setup(fullmap, othermap, partialmap);
}

FOUR_C_NAMESPACE_CLOSE
