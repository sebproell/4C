// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_PARAMETER_CONTAINER_HPP
#define FOUR_C_IO_INPUT_PARAMETER_CONTAINER_HPP


#include "4C_config.hpp"

#include "4C_utils_demangle.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

#include <any>
#include <functional>
#include <map>
#include <optional>
#include <ostream>
#include <string>
#include <typeindex>
#include <vector>


FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * A type that either contains a value of type T or nothing. This type is useful to
   * represent the absence of a value marked by "none" in the input file. The type is just an alias
   * for std::optional<T> to more precisely convey the intended meaning.
   */
  template <typename T>
  using Noneable = std::optional<T>;

  /**
   * A constant to represent the absence of a value in the input file. The constant is templated
   * on the type of the value inside the Noneable to properly handle cases where a Noneable is
   * itself stored inside a std::optional type.
   */
  template <typename T>
  inline constexpr auto none = Noneable<T>{std::nullopt};

  /**
   * A container to store dynamic input parameters. The container can store arbitrary types of
   * parameters. Parameters can be grouped in sub-containers.
   *
   * This class is a core part of the input mechanism, as it contains the parsed data from the input
   * file and grants access to it.
   */
  class InputParameterContainer
  {
   public:
    /**
     * A type to store a list of InputParameterContainers. This type is used to represent what is
     * often called a *list*, *array*, or *sequence* of data in the input file.
     *
     * @note This type is called `List` to more clearly distinguish it from a simple `std::vector`
     * entry in the container. A `List` contains nested InputParameterContainers and encodes
     * rather complex data structures. Nevertheless, it is implemented as a `std::vector`.
     */
    using List = std::vector<InputParameterContainer>;

    /**
     * \brief Add @data to the container at the given key @name.
     *
     * If an entry with given @p name already exists, it will be overwritten.
     */
    template <typename T>
    void add(const std::string& name, const T& data);

    /**
     * Access group @p name. If the group does not exist, it will be created.
     */
    InputParameterContainer& group(const std::string& name);

    /**
     * Access group @p name. This function throws an error if the group does not exist.
     */
    [[nodiscard]] const InputParameterContainer& group(const std::string& name) const;

    /**
     * Check if a group with the name @p name exists.
     */
    [[nodiscard]] bool has_group(const std::string& name) const;

    /**
     * Add the list @p data at the given key @p name.
     *
     * @note This functions is a more obvious way to add a list to the container compared to
     * the add() function with a List template argument, although this is precisely what happens
     * internally.
     */
    void add_list(const std::string& name, List&& list);

    /**
     * Access the list @p name. This function throws an error if the list does not exist.
     *
     * @note This functions is a more obvious way to get a list from the container compared to
     * the get() function with a List template argument, although this is precisely
     * what happens internally.
     */
    [[nodiscard]] const List& get_list(const std::string& name) const;

    /**
     * Combine the data from another container with this one. Conflicting data will throw an
     * error.
     */
    void merge(const InputParameterContainer& other);

    /*!
     * Get a const reference to the data stored at the key @p name from the container. An error
     * is thrown in case no value of specified type is stored under @p name in the container.
     */
    template <typename T>
    const T& get(const std::string& name) const;

    /*!
     * Get the data stored at the key @p name from the container. Return the @p default_value if
     * no value of specified type is stored under @p name in the container.
     *
     * @note This function returns the value as a copy.
     */
    template <typename T>
    T get_or(const std::string& name, T default_value) const;

    /*!
     * Get the const pointer to the data stored at the key @p name from the container. Return a
     * `nullptr` if no value of specified type is stored under @p name in the container.
     */
    template <typename T>
    const T* get_if(const std::string& name) const;

    /**
     * Print the data in the container to the output stream @p os.
     */
    void print(std::ostream& os) const;

    /**
     * Clear the container.
     */
    void clear();

    /**
     * Convert the data in this container to a Teuchos::ParameterList. All groups are converted to
     * sublists.
     */
    void to_teuchos_parameter_list(Teuchos::ParameterList& list) const;

   private:
    //! Entry stored in the container.
    struct Entry
    {
      //! The actual data.
      std::any data;
    };

    //! Data stored in this container.
    std::map<std::string, Entry> entries_;

    //! Groups present in this container. Groups are InputParameterContainers themselves.
    std::map<std::string, InputParameterContainer> groups_;

    /**
     * Gather different actions that can be performed on a type.
     */
    struct TypeActions
    {
      //! The function to print the data.
      std::function<void(std::ostream&, const std::any&)> print;

      //! Function to write the data into a Teuchos::ParameterList.
      std::function<void(Teuchos::ParameterList&, const std::string& name, const std::any&)>
          write_to_pl;
    };

    /**
     * Add the type actions if not already present.
     */
    template <typename T>
    void ensure_type_action_registered();

    /**
     * Access the shared storage for all the type actions.
     */
    static std::map<std::type_index, TypeActions>& get_type_actions();
  };  // class InputParameterContainer

}  // namespace Core::IO


// --- template and inline functions ---//


namespace Core::IO::Internal::InputParameterContainerImplementation
{
  template <typename T>
  const T* try_get_any_data(const std::string& name, const std::any& data)
  {
    if (typeid(T) == data.type())
    {
      const T* any_ptr = std::any_cast<T>(&data);
      FOUR_C_ASSERT(any_ptr != nullptr, "Implementation error.");
      return any_ptr;
    }
    else
    {
      FOUR_C_THROW(
          "You tried to get the data named %s from the container as type '%s'.\n"
          "Actually, it has type '%s'.",
          name.c_str(), Core::Utils::get_type_name<T>().c_str(),
          Core::Utils::try_demangle(data.type().name()).c_str());
    }
  }

  // Default printer if not printable.
  template <typename T>
  struct PrintHelper
  {
    void operator()(std::ostream& os, const std::any& data) { os << "<not printable> "; }
  };

  template <typename T>
  concept StreamInsertable = requires(std::ostream& os, const T& t) { os << t; };

  // Specialization for stream insert.
  template <StreamInsertable T>
  struct PrintHelper<T>
  {
    void operator()(std::ostream& os, const std::any& data) { os << std::any_cast<T>(data) << " "; }
  };

  template <StreamInsertable T>
  struct PrintHelper<Noneable<T>>
  {
    void operator()(std::ostream& os, const std::any& data)
    {
      auto val = std::any_cast<Noneable<T>>(data);
      if (val.has_value())
        PrintHelper<T>{}(os, *val);
      else
        os << "none ";
    }
  };

  // Specialization for vectors.
  template <typename T>
  struct PrintHelper<std::vector<T>>
  {
    void operator()(std::ostream& os, const std::any& data)
    {
      FOUR_C_ASSERT(typeid(std::vector<T>) == data.type(), "Implementation error.");
      const auto& vec = std::any_cast<std::vector<T>>(data);
      for (const auto& v : vec)
      {
        PrintHelper<T>{}(os, v);
      }
    }
  };

  // Specialization for maps.
  template <typename Key, typename Value>
  struct PrintHelper<std::map<Key, Value>>
  {
    void operator()(std::ostream& os, const std::any& data)
    {
      FOUR_C_ASSERT(typeid(std::map<Key, Value>) == data.type(), "Implementation error.");
      const auto& map = std::any_cast<std::map<Key, Value>>(data);
      for (const auto& [key, value] : map)
      {
        os << key << " : ";
        PrintHelper<Value>{}(os, value);
      }
    }
  };

}  // namespace Core::IO::Internal::InputParameterContainerImplementation


template <typename T>
void Core::IO::InputParameterContainer::add(const std::string& name, const T& data)
{
  entries_[name] = {.data = std::any{data}};
  ensure_type_action_registered<T>();
}

template <typename T>
const T& Core::IO::InputParameterContainer::get(const std::string& name) const
{
  if (const T* p = get_if<T>(name))
    return *p;
  else
    FOUR_C_THROW("Key '%s' cannot be found in the container.", name.c_str());
}

template <typename T>
T Core::IO::InputParameterContainer::get_or(const std::string& name, T default_value) const
{
  if (const T* p = get_if<T>(name))
    return *p;
  else
    return default_value;
}

template <typename T>
const T* Core::IO::InputParameterContainer::get_if(const std::string& name) const
{
  const auto it = entries_.find(name);
  if (it != entries_.end())
  {
    return Internal::InputParameterContainerImplementation::try_get_any_data<T>(
        name, it->second.data);
  }
  else
  {
    return nullptr;
  }
}

template <typename T>
void Core::IO::InputParameterContainer::ensure_type_action_registered()
{
  auto& type_actions = get_type_actions();
  if (!type_actions.contains(typeid(T)))
  {
    type_actions[typeid(T)] = {
        .print = Internal::InputParameterContainerImplementation::PrintHelper<T>{},
        .write_to_pl = [](Teuchos::ParameterList& pl, const std::string& name, const std::any& data)
        { pl.set<T>(name, std::any_cast<T>(data)); },
    };
  }
}

FOUR_C_NAMESPACE_CLOSE

#endif
