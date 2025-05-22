// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_MAP_HPP
#define FOUR_C_LINALG_MAP_HPP

#include "4C_config.hpp"

#include "4C_linalg_view.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_owner_or_view.hpp"

#include <Epetra_Map.h>
#include <mpi.h>

#include <memory>
#include <variant>

FOUR_C_NAMESPACE_OPEN

// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Epetra_Map

// NOLINTBEGIN(readability-identifier-naming)



namespace Core::LinAlg
{
  class Map
  {
    // This wrapper may have two different variants
    // Either it holds an Epetra_Map or Epetra_BlockMap
    using MapVariant =
        std::variant<Utils::OwnerOrView<Epetra_Map>, Utils::OwnerOrView<Epetra_BlockMap>>;

    // helper function to access each variant safely
    template <typename Func>
    decltype(auto) visit_variant(Func&& func) const
    {
      return std::visit(
          [&](const auto& wrapped) -> decltype(auto)
          {
            return func(*wrapped);  // Unwrap OwnerOrView
          },
          map_);
    }

   public:
    Map(int NumGlobalElements, int IndexBase, const MPI_Comm& Comm);

    Map(int NumGlobalElements, int NumMyElements, int IndexBase, const MPI_Comm& Comm);

    Map(int NumGlobalElements, int NumMyElements, const int* MyGlobalElements, int IndexBase,
        const MPI_Comm& Comm);

    Map(const Map& Source);

    /// Copy constructor from Epetra_Map
    explicit Map(const Epetra_Map& Source);

    /// Copy constructor from Epetra_BlockMap
    explicit Map(const Epetra_BlockMap& Source);

    ~Map() = default;

    Map& operator=(const Map& other);

    //! Print object to the output stream
    void Print(std::ostream& os) const
    {
      visit_variant([&](const auto& map) -> void { map.Print(os); });
    }

    //! Returns a reference of the Epetra_Map if available.
    const Epetra_Map& get_epetra_map() const
    {
      auto* map = std::get_if<Utils::OwnerOrView<Epetra_Map>>(&map_);
      if (map == nullptr)
      {
        FOUR_C_THROW(
            "This Map is based on an Epetra_BlockMap, not an Epetra_Map. This cast is not "
            "possible.");
      }
      return **map;
    }

    //! Returns a reference of the Epetra_Map if available.
    Epetra_Map& get_epetra_map()
    {
      auto* map = std::get_if<Utils::OwnerOrView<Epetra_Map>>(&map_);
      if (map == nullptr)
      {
        FOUR_C_THROW(
            "This Map is based on an Epetra_BlockMap, not an Epetra_Map. This cast is not "
            "possible.");
      }
      return **map;
    }


    //! Returns a const reference to the underlying Epetra_BlockMap.
    const Epetra_BlockMap& get_epetra_block_map() const
    {
      // If the map_ holds directly the Epetra_BlockMap, just return it.
      if (auto map = std::get_if<Utils::OwnerOrView<Epetra_BlockMap>>(&map_))
      {
        return **map;
      }

      // If the Epetra_Map is stored try to cast it.
      const auto& map = *std::get<Utils::OwnerOrView<Epetra_Map>>(map_);

      return map;
    }


    //! Returns a reference to the underlying Epetra_BlockMap.
    Epetra_BlockMap& get_epetra_block_map()
    {
      // If the map_ holds directly the Epetra_BlockMap, just return it.
      if (auto map = std::get_if<Utils::OwnerOrView<Epetra_BlockMap>>(&map_))
      {
        return **map;
      }

      // If the Epetra_Map is stored try to cast it.
      auto& map = *std::get<Utils::OwnerOrView<Epetra_Map>>(map_);
      return map;
    }

    //! Returns true if this and Map are identical maps
    bool SameAs(const Map& other) const
    {
      return std::visit(
          [&](const auto& this_wrapped)
          {
            return std::visit([&](const auto& other_wrapped)
                { return (*this_wrapped).SameAs(*other_wrapped); }, other.map_);
          },
          map_);
    }

    //! Returns true if this and Map have identical point-wise structure
    bool PointSameAs(const Map& Map) const
    {
      return visit_variant(
          [&](const auto& map) { return map.PointSameAs(Map.get_epetra_block_map()); });
    }

    //! Number of elements across all processors.
    int NumGlobalElements() const
    {
      return visit_variant([](const auto& map) { return map.NumGlobalElements(); });
    }

    //! Number of elements on the calling processor.
    int NumMyElements() const
    {
      return visit_variant([](const auto& map) { return map.NumMyElements(); });
    }

    //! returns the index base for this map.
    int IndexBase() const
    {
      return visit_variant([](const auto& map) { return map.IndexBase(); });
    }

    //! Pointer to internal array containing a mapping between the local elements and the first
    //! local point number in each element.
    int FirstPointInElementList(int* LID) const
    {
      return visit_variant([&](const auto& map) { return map.FirstPointInElementList(LID); });
    }

    //! Returns true if map is defined across more than one processor.
    bool DistributedGlobal() const
    {
      return visit_variant([](const auto& map) { return map.DistributedGlobal(); });
    }

    //! Returns true if this and Map are identical maps
    bool SameAs(const Epetra_Map& other) const
    {
      return visit_variant([&](const auto& map) { return map.SameAs(other); });
    }

    //! Returns true if this and Map are identical maps
    bool SameAs(const Epetra_BlockMap& other) const
    {
      return visit_variant([&](const auto& map) { return map.SameAs(other); });
    }

    //! Returns true if the GID passed in belongs to the calling processor in this map, otherwise
    //! returns false.
    bool MyGID(int GID_in) const
    {
      return visit_variant([&](const auto& map) { return map.MyGID(GID_in); });
    }

    //! Returns global ID of local ID, return IndexBase-1 if not found on this processor.
    int GID(int LID) const
    {
      return visit_variant([&](const auto& map) { return map.GID(LID); });
    }

    //! Returns the size of elements in the map; only valid if map has constant element size.
    int ElementSize() const
    {
      return visit_variant([](const auto& map) { return map.ElementSize(); });
    }

    //! Returns the size of elements in the map; only valid if map has constant element size.
    int ElementSize(int LID) const
    {
      return visit_variant([&](const auto& map) { return map.ElementSize(LID); });
    }

    //! Returns the maximum global ID across the entire map.
    int MaxAllGID() const
    {
      return visit_variant([](const auto& map) { return map.MaxAllGID(); });
    }

    //! Returns the minimum global ID across the entire map.
    int MinAllGID() const
    {
      return visit_variant([](const auto& map) { return map.MinAllGID(); });
    }

    //! Returns local ID of global ID, return -1 if not found on this processor.
    int LID(int GID) const
    {
      return visit_variant([&](const auto& map) { return map.LID(GID); });
    }

    //! Returns true if this and Map have identical point-wise structure
    bool PointSameAs(const Epetra_Map& Map) const
    {
      return visit_variant([&](const auto& map) { return map.PointSameAs(Map); });
    }

    //! Returns true if this and Map have identical point-wise structure
    bool PointSameAs(const Epetra_BlockMap& Map) const
    {
      return visit_variant([&](const auto& map) { return map.PointSameAs(Map); });
    }

    //! Returns the processor IDs and corresponding local index value for a given list of global
    //! indices
    int RemoteIDList(int NumIDs, int* GIDList, int* PIDList, int* LIDList) const
    {
      return visit_variant(
          [&](const auto& map) { return map.RemoteIDList(NumIDs, GIDList, PIDList, LIDList); });
    }

    //! Returns the processor IDs and corresponding local index value for a given list of global
    //! indices
    int RemoteIDList(
        int NumIDs, const int* GIDList, int* PIDList, int* LIDList, int* SizeList) const
    {
      return visit_variant([&](const auto& map)
          { return map.RemoteIDList(NumIDs, GIDList, PIDList, LIDList, SizeList); });
    }

    //! Returns the minimum global ID owned by this processor.
    int MinMyGID(void) const
    {
      return visit_variant([](const auto& map) { return map.MinMyGID(); });
    }

    //! Access function for Epetra_Comm communicator.
    MPI_Comm Comm() const;

    //! Returns true if map GIDs are 1-to-1.
    bool UniqueGIDs(void) const
    {
      return visit_variant([](const auto& map) { return map.UniqueGIDs(); });
    }

    //! Pointer to internal array containing list of global IDs assigned to the calling processor.
    int* MyGlobalElements(void) const
    {
      return visit_variant([](const auto& map) { return map.MyGlobalElements(); });
    }

    //! Maximum element size across all processors.
    int MaxElementSize(void) const
    {
      return visit_variant([](const auto& map) { return map.MaxElementSize(); });
    }

    //! Puts list of global elements on this processor into the user-provided array.
    int MyGlobalElements(int* MyGlobalElementList) const
    {
      return visit_variant(
          [&](const auto& map) { return map.MyGlobalElements(MyGlobalElementList); });
    }
    //! Number of local points for this map; equals the sum of all element sizes on the calling
    //! processor.
    int NumMyPoints() const
    {
      return visit_variant([](const auto& map) { return map.NumMyPoints(); });
    }

    //! Returns a pointer to the BlockMapData instance this BlockMap uses.
    const Epetra_BlockMapData* DataPtr() const
    {
      return visit_variant([](const auto& map) { return map.DataPtr(); });
    }


    [[nodiscard]] static std::unique_ptr<Map> create_view(Epetra_Map& view);
    [[nodiscard]] static std::unique_ptr<const Map> create_view(const Epetra_Map& view);
    [[nodiscard]] static std::unique_ptr<Map> create_view(Epetra_BlockMap& view);
    [[nodiscard]] static std::unique_ptr<const Map> create_view(const Epetra_BlockMap& view);

   private:
    Map() = default;

    //! stores an Epetra_BlockMap or Epetra_Map
    MapVariant map_;
  };

  inline std::ostream& operator<<(std::ostream& os, const Map& m)
  {
    os << m.get_epetra_block_map();
    return os;
  }

  template <>
  struct EnableViewFor<Epetra_Map>
  {
    using type = Map;
  };

  template <>
  struct EnableViewFor<Epetra_BlockMap>
  {
    using type = Map;
  };

}  // namespace Core::LinAlg

// NOLINTEND(readability-identifier-naming)

FOUR_C_NAMESPACE_CLOSE

#endif
