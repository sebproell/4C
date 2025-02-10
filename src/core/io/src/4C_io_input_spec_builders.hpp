// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_SPEC_BUILDERS_HPP
#define FOUR_C_IO_INPUT_SPEC_BUILDERS_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_io_yaml.hpp"
#include "4C_utils_string.hpp"

#include <algorithm>
#include <functional>
#include <optional>
#include <ostream>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace Internal
  {
    /**
     * Helper to print values in the format of a dat file.
     */
    struct DatPrinter
    {
      template <typename T>
      void operator()(std::ostream& out, const T& val) const
      {
        out << val;
      }

      void operator()(std::ostream& out, bool val) const { out << " " << (val ? "true" : "false"); }

      template <typename T>
      void operator()(std::ostream& out, const std::optional<T>& val) const
      {
        if (val.has_value())
          (*this)(out, *val);
        else
          out << "none";
      }

      template <typename T>
      void operator()(std::ostream& out, const std::vector<T>& val) const
      {
        for (const auto& v : val)
        {
          (*this)(out, v);
          out << " ";
        }
      }

      template <typename T, typename U>
      void operator()(std::ostream& out, const std::pair<T, U>& val) const
      {
        (*this)(out, val.first);
        out << " ";
        (*this)(out, val.second);
      }
    };

    template <typename T>
    concept CustomDatPrintable = requires(const T& t, std::ostream& stream, std::size_t indent) {
      { t.print(stream, indent) } -> std::same_as<void>;
    };

    template <typename T, typename AlwaysVoid = void>
    constexpr bool has_set_default_value = false;

    template <typename T>
    constexpr bool has_set_default_value<T,
        std::void_t<decltype(std::declval<const std::decay_t<T>>().set_default_value(
            std::declval<Core::IO::InputParameterContainer&>()))>> = true;


    template <typename T, typename AlwaysVoid = void>
    struct PrettyTypeName
    {
      std::string operator()() { return Utils::try_demangle(typeid(T).name()); }
    };

    template <>
    struct PrettyTypeName<std::string>
    {
      std::string operator()() { return "string"; }
    };

    template <>
    struct PrettyTypeName<std::filesystem::path>
    {
      std::string operator()() { return "path"; }
    };

    template <typename T>
    struct PrettyTypeName<std::vector<T>>
    {
      std::string operator()() { return "vector<" + PrettyTypeName<T>{}() + ">"; }
    };

    template <typename T, typename U>
    struct PrettyTypeName<std::pair<T, U>>
    {
      std::string operator()()
      {
        return "pair<" + PrettyTypeName<T>{}() + ", " + PrettyTypeName<U>{}() + ">";
      }
    };

    template <typename T>
    struct PrettyTypeName<Noneable<T>>
    {
      std::string operator()() { return "Noneable<" + PrettyTypeName<T>{}() + ">"; }
    };

    template <typename T>
    std::string get_pretty_type_name()
    {
      return PrettyTypeName<T>{}();
    }

    class MatchTree;

    /**
     * Entries in the MatchTree.
     */
    struct MatchEntry
    {
      MatchTree* tree;
      const InputSpec* spec;
      std::vector<MatchEntry*> children;

      /**
       * A MatchEntry can only match a single node. Logical specs like all_of and one_of are not
       * considered to match nodes themselves, as this is done by their children.
       */
      ryml::id_type matched_node{ryml::npos};

      enum class State : std::uint8_t
      {
        unmatched,
        matched,
        partial,
        defaulted,
      };

      enum class Type : std::uint8_t
      {
        unknown,
        entry,
        group,
        all_of,
        one_of,
      };

      State state{State::unmatched};
      Type type{Type::unknown};


      /**
       * Append a child for the @p in_spec to the current entry and return a reference to it. This
       * child is passed on to the match function of @p in_spec.
       */
      MatchEntry& append_child(const InputSpec* in_spec);
    };

    /**
     * A tree that tracks how well a spec matches a yaml tree.
     */
    class MatchTree
    {
     public:
      MatchTree(const InputSpec& root);

      MatchEntry& root() { return entries_.front(); }

      MatchEntry& append_child(const InputSpec* spec);

      void dump(std::ostream& stream) const;

      /**
       * Throw an exception that contains the input and match tree in a user-friendly format, if
       * the match was not successful.
       */
      void assert_match(ryml::ConstNodeRef node) const;

     private:
      std::vector<MatchEntry> entries_;
    };

    class InputSpecTypeErasedBase
    {
     public:
      /**
       * Store a reduced version of the common data entries.
       */
      struct CommonData
      {
        /**
         * The name or key under which the value is found in the input. This name is also used to
         * store the value in the container.
         */
        std::string name;

        /**
         * An optional description of the value.
         */
        std::string description;

        /**
         * Whether the value is required or optional.
         */
        bool required;

        /**
         * Whether the spec has a default value.
         */
        bool has_default_value;

        /**
         * The total number of specs that make up this spec. This includes the spec itself, meaning,
         * that the minimum value is 1. This value can used to reserve memory ahead of time.
         */
        std::size_t n_specs;
      };

      virtual ~InputSpecTypeErasedBase() = default;

      /**
       * @param data The common data of the spec.
       * @param n_specs The number of specs (itself included) that make up this spec.
       */
      InputSpecTypeErasedBase(CommonData data) : data(std::move(data)) {}

      virtual void parse(ValueParser& parser, InputParameterContainer& container) const = 0;

      /**
       * Returns true if the node matches the spec and stores the value in the container. The passed
       * @p node is the parent node which might contain data matching the spec. Every spec needs
       * to check if it can find required data in this node. If yes, the InputSpec should report
       * itself as matched in the @p match_entry. Note that the @p match_entry already refers to the
       * spec that is being matched. When matching more specs internally, the spec needs to append
       * children to the match_entry.
       */
      virtual bool match(ConstYamlNodeRef node, InputParameterContainer& container,
          MatchEntry& match_entry) const = 0;

      virtual void set_default_value(InputParameterContainer& container) const = 0;

      //! Emit metadata. This function always emits into a map, i.e., the implementation must
      //! insert keys and values into the yaml emitter.
      virtual void emit_metadata(ryml::NodeRef node) const = 0;

      [[nodiscard]] virtual std::string pretty_type_name() const = 0;

      [[nodiscard]] virtual std::unique_ptr<InputSpecTypeErasedBase> clone() const = 0;

      void print(std::ostream& stream, std::size_t indent) const { do_print(stream, indent); }

      [[nodiscard]] const std::string& name() const { return data.name; }

      [[nodiscard]] const std::string& description() const { return data.description; }

      [[nodiscard]] std::string description_one_line() const;

      [[nodiscard]] bool required() const { return data.required; }

      [[nodiscard]] bool has_default_value() const { return data.has_default_value; }


      CommonData data;

     protected:
      InputSpecTypeErasedBase(const InputSpecTypeErasedBase&) = default;
      InputSpecTypeErasedBase& operator=(const InputSpecTypeErasedBase&) = default;
      InputSpecTypeErasedBase(InputSpecTypeErasedBase&&) noexcept = default;
      InputSpecTypeErasedBase& operator=(InputSpecTypeErasedBase&&) noexcept = default;

     private:
      virtual void do_print(std::ostream& stream, std::size_t indent) const = 0;
    };

    template <typename T>
    concept StoresType = requires(T t) { typename T::DataType::StoredType; };

    template <typename T>
    struct InputSpecTypeErasedImplementation : public InputSpecTypeErasedBase
    {
      template <typename T2>
      explicit InputSpecTypeErasedImplementation(T2&& wrapped, CommonData data)
          : InputSpecTypeErasedBase(std::move(data)), wrapped(std::forward<T2>(wrapped))
      {
      }

      void parse(ValueParser& parser, InputParameterContainer& container) const override
      {
        wrapped.parse(parser, container);
      }

      bool match(ConstYamlNodeRef node, InputParameterContainer& container,
          MatchEntry& match_entry) const override
      {
        return wrapped.match(node, container, match_entry);
      }

      void set_default_value(InputParameterContainer& container) const override
      {
        if constexpr (Internal::has_set_default_value<T>)
        {
          wrapped.set_default_value(container);
        }
        else
        {
          FOUR_C_ASSERT(wrapped.data.default_value.has_value(),
              "Implementation error: this function should only be called if the wrapped type has "
              "an optional default value.");

          container.add(wrapped.name, *wrapped.data.default_value);
        }
      }

      void emit_metadata(ryml::NodeRef node) const override { wrapped.emit_metadata(node); }

      void do_print(std::ostream& stream, std::size_t indent) const override
      {
        if constexpr (CustomDatPrintable<T>)
        {
          wrapped.print(stream, indent);
        }
        else
        {
          stream << "// " << std::string(indent, ' ') << name();

          // pretty printed type of the parameter
          stream << " <" << Internal::get_pretty_type_name<typename T::DataType::StoredType>()
                 << ">";

          if (!required())
          {
            stream << " (optional)";
          }

          if (has_default_value())
          {
            stream << " (default: ";
            Internal::DatPrinter{}(stream, wrapped.data.default_value.value());
            stream << ")";
          }

          if (!description().empty())
          {
            stream << " " << std::quoted(description_one_line());
          }
          stream << "\n";
        }
      }

      std::string pretty_type_name() const override
      {
        if constexpr (StoresType<T>)
        {
          return Internal::get_pretty_type_name<typename T::DataType::StoredType>();
        }
        else
        {
          FOUR_C_ASSERT(false,
              "Implementation error: this function should only be called for "
              "types that store a type.");
          return "unknown";
        }
      }

      [[nodiscard]] std::unique_ptr<InputSpecTypeErasedBase> clone() const override
      {
        return std::make_unique<InputSpecTypeErasedImplementation<T>>(wrapped, data);
      }

      T wrapped;
    };

    template <typename T>
    InputSpec make_spec(T&& wrapped, InputSpecTypeErasedBase::CommonData data)
    {
      return InputSpec(std::make_unique<InputSpecTypeErasedImplementation<std::decay_t<T>>>(
          std::forward<T>(wrapped), std::move(data)));
    }


    template <typename T>
    struct SupportedTypeHelper : std::false_type
    {
    };

    template <typename T>
    concept SupportedTypePrimitives =
        std::same_as<T, int> || std::same_as<T, double> || std::same_as<T, bool> ||
        std::same_as<T, std::string> || std::same_as<T, std::filesystem::path> || std::is_enum_v<T>;

    template <SupportedTypePrimitives T>
    struct SupportedTypeHelper<T> : std::true_type
    {
    };

    template <typename T>
    struct SupportedTypeHelper<std::vector<T>> : SupportedTypeHelper<T>
    {
    };

    template <typename T, typename U>
    struct SupportedTypeHelper<std::pair<T, U>> : SupportedTypeHelper<T>, SupportedTypeHelper<U>
    {
    };

    template <typename T>
    struct SupportedTypeHelper<Noneable<T>> : SupportedTypeHelper<T>
    {
    };

  }  // namespace Internal

  /**
   * Helper functions to create InputSpec objects. When you want to create an InputSpec, you
   * can "import" the functions from this namespace by using the following line:
   *
   * @code
   * using namespace Core::IO::InputSpecBuilders;
   * @endcode
   *
   * This allows you to create InputSpec objects in a concise notation like this:
   *
   * @code
   * auto input = group("params",
   *   {
   *   entry<int>("a", {.description = "An integer value", .required = true}),
   *   entry<std::vector<double>>("b", {.description = "A vector of doubles", .required = false,
   *     .size = 4}),
   *   }
   * );
   * @endcode
   */
  namespace InputSpecBuilders
  {
    /**
     * We deliberately limit ourselves to a few generally useful types. While it would not be too
     * difficult to support all the fundamental and container types that C++ provides, this would
     * likely lead to more confusion for users than it would provide benefits. After all, when
     * consuming the parsed input, the user will have to use the exact type of the parameter. Also,
     * input file formats are often not able to distinguish fundamental types like `double` and
     * `float` and there is little benefit in supporting both in the input mechanism. Any conversion
     * between types can be done in the user code, which usually entails additional validation and
     * error handling anyway.
     *
     * The supported types are:
     * - `int`
     * - `double`
     * - `bool`
     * - `std::string`
     * - `std::filesystem::path`
     * - `enum` and `enum class` types
     * - `Noneable<T>`, where `T` is one of the supported types
     * - `std::vector<T>`, where `T` is one of the supported types
     * - `std::pair<T, U>`, where `T` and `U` are supported types
     */
    template <typename T>
    concept SupportedType = Internal::SupportedTypeHelper<T>::value;

    // Import the Noneable type into the InputSpecBuilders namespace to make it easier to use.
    using Core::IO::none;
    using Core::IO::Noneable;

    /**
     * Callback function that may be attached to entry().
     */
    using EntryCallback = std::function<void(InputParameterContainer&)>;

    //! Additional parameters for a scalar-valued entry().
    template <typename StoredTypeIn>
    struct ScalarData
    {
      using StoredType = StoredTypeIn;
      /**
       * An optional description of the value.
       */
      std::string description{};

      /**
       * Whether the value is required or optional.
       */
      std::optional<bool> required{};

      /**
       * The default value of the parameter. If this optional fields is set, the parameter is
       * optional. If it is not set, the parameter is required.
       */
      std::optional<StoredType> default_value{};

      /**
       * An optional callback that is called after the value has been parsed. This can be used to
       * set additional values in the container.
       */
      EntryCallback on_parse_callback{nullptr};
    };

    //! Additional parameters for a vector-valued entry().
    template <typename StoredTypeIn, typename ScalarTypeIn>
    struct VectorData
    {
      using StoredType = StoredTypeIn;
      using ScalarType = ScalarTypeIn;

      /**
       * A function that determines the size of the vector. This function is called with the
       * already parsed content of the input line, which may be used to query the size as the
       * value of another parameter.
       */
      using SizeCallback = std::function<int(const InputParameterContainer&)>;

      /**
       * An optional description of the value.
       */
      std::string description{};

      /**
       * Whether the value is required or optional.
       */
      std::optional<bool> required{};

      /**
       * The default value of the parameter. If this optional fields is set, the parameter is
       * optional. If it is not set, the parameter is required.
       */
      std::optional<StoredType> default_value{};

      /**
       * The size of the data. This can either be a fixed size or a callback that determines the
       * size based on the value of another parameter.
       *
       * @note We use `int` instead of `size_t` since this is what most users write. Template
       *      argument deduction will work a lot better with `int` than with `size_t`.
       */
      std::variant<int, SizeCallback> size;

      /**
       * An optional callback that is called after the value has been parsed. This can be used to
       * set additional values in the container.
       */
      EntryCallback on_parse_callback{nullptr};
    };

    //! Additional parameters for a group().
    struct GroupData
    {
      /**
       * An optional description of the Group.
       */
      std::string description{};

      /**
       * Whether the Group is required or optional.
       */
      bool required{true};
    };

    namespace Internal
    {
      template <typename T>
      struct DataForHelper
      {
        using type = ScalarData<T>;
      };

      template <typename T>
      struct DataForHelper<std::vector<T>>
      {
        using type = VectorData<std::vector<T>, T>;
      };

      template <typename T>
      using DataFor = typename DataForHelper<T>::type;


      template <typename DataType>
      static constexpr bool is_sized_data = false;

      template <typename D, typename S>
      static constexpr bool is_sized_data<VectorData<D, S>> = true;

      //! Make .required field consistent with .default_value.
      template <typename DataType>
      void sanitize_required_default(DataType& data)
      {
        if (data.default_value.has_value())
        {
          if (data.required.has_value())
          {
            FOUR_C_ASSERT_ALWAYS(!data.required.value(),
                "A parameter cannot be both required and have a default value.");
          }
          else
          {
            data.required = false;
          }
        }
        else
        {
          if (!data.required.has_value()) data.required = true;
        }
        FOUR_C_ASSERT(data.required.has_value(), "Required field must now be set.");
      }

      template <typename DataTypeIn>
      struct BasicSpec
      {
        std::string name;
        using DataType = std::decay_t<DataTypeIn>;
        using StoredType = typename DataType::StoredType;
        DataType data;
        void parse(ValueParser& parser, InputParameterContainer& container) const;
        bool match(ConstYamlNodeRef node, InputParameterContainer& container,
            IO::Internal::MatchEntry& match_entry) const;
        void emit_metadata(ryml::NodeRef node) const;
      };


      template <typename DataTypeIn>
      struct SelectionSpec
      {
        std::string name;
        using DataType = std::decay_t<DataTypeIn>;
        using StoredType = typename DataType::StoredType;
        DataType data;
        std::vector<std::pair<std::string, StoredType>> choices;
        void parse(ValueParser& parser, InputParameterContainer& container) const;
        bool match(ConstYamlNodeRef node, InputParameterContainer& container,
            IO::Internal::MatchEntry& match_entry) const;
        void print(std::ostream& stream, std::size_t indent) const;
        void emit_metadata(ryml::NodeRef node) const;
      };

      struct GroupSpec
      {
        std::string name;
        GroupData data;
        std::vector<InputSpec> specs;

        void parse(ValueParser& parser, InputParameterContainer& container) const;
        bool match(ConstYamlNodeRef node, InputParameterContainer& container,
            IO::Internal::MatchEntry& match_entry) const;
        void set_default_value(InputParameterContainer& container) const;
        void print(std::ostream& stream, std::size_t indent) const;
        void emit_metadata(ryml::NodeRef node) const;
      };

      struct AllOfSpec
      {
        GroupData data;
        std::vector<InputSpec> specs;

        void parse(ValueParser& parser, InputParameterContainer& container) const;
        bool match(ConstYamlNodeRef node, InputParameterContainer& container,
            IO::Internal::MatchEntry& match_entry) const;
        void set_default_value(InputParameterContainer& container) const;
        void print(std::ostream& stream, std::size_t indent) const;
        void emit_metadata(ryml::NodeRef node) const;
      };

      struct OneOfSpec
      {
        // A one_of spec is essentially an unnamed group with additional logic to ensure that
        // exactly one of the contained specs is present.
        GroupData data;
        std::vector<InputSpec> specs;

        //! This callback may be used to perform additional actions after parsing one of the specs.
        //! The index of the parsed spec as given inside #specs is passed as an argument.
        std::function<void(InputParameterContainer& container, std::size_t index)>
            on_parse_callback;

        void parse(ValueParser& parser, InputParameterContainer& container) const;

        bool match(ConstYamlNodeRef node, InputParameterContainer& container,
            IO::Internal::MatchEntry& match_entry) const;

        void set_default_value(InputParameterContainer& container) const;

        void print(std::ostream& stream, std::size_t indent) const;

        void emit_metadata(ryml::NodeRef node) const;
      };


      //! Helper to create selection() specs.
      template <typename T, typename DataType = Internal::DataFor<T>>
      [[nodiscard]] InputSpec selection_internal(
          std::string name, std::vector<std::pair<std::string, T>> choices, DataType data = {});
    }  // namespace Internal

    /**
     * Create a normal entry with given @p name. All entries are parameterized by a struct which
     * contains the optional `description`, `required` and `default_value` fields. The following
     * examples demonstrate how entries can be created:
     *
     * @code
     * // An entry with name and description. By default, the entry is required.
     * entry<std::string>("my_string", {.description = "A string value."});
     *
     * // An entry with a default value. This entry is implicitly optional because a default value
     * // is given.
     * entry<double>("my_double", {.default_value = 3.14});
     * // This is equivalent to:
     * entry<double>("my_double", {.required = false, .default_value = 3.14});
     *
     * // An optional entry. This value does not have a default value.
     * entry<int>("my_int", {.required = false});
     *
     * // An alternative way to create this optional int entry is achieved with the Noneable type.
     * // This value is optional and by default has an empty value represented by "none" in the
     * // input file.
     * entry<Noneable<int>>("my_int", .default_value = none<int>);
     *
     * // A vector entry with a fixed size of 3.
     * entry<std::vector<double>>("my_vector", {.size = 3});
     *
     * // A vector may also contain Noneable values.
     * entry<std::vector<Noneable<double>>>("my_vector", {.size = 3});
     *
     * // A vector entry with a size that is determined by the value of another parameter. The
     * // size is given as a callback.
     * entry<int>("N");
     * entry<std::vector<double>>("my_vector", {.size = from_parameter<int>("N")});
     *
     * // A vector entry which performs an additional action after parsing.
     * entry<std::filesystem::path>("data_file", {.description = "A path to a file.",
     *   .on_parse_callback = [](InputParameterContainer& container) {
     *     // Perform an action with the parsed path.
     *     std::filesystem::path path = container.get<std::filesystem::path>("my_path");
     *     // e.g. read the file content and also store it in the container.
     *     // auto data_table = ...
     *     container.add("data_table", data_table);
     *   }});
     * @endcode
     *
     * After parsing an InputSpec with fully_parse(), the value of the entry can be retrieved from
     * an InputParameterContainer. The details depend on the `required` and `default_value` fields:
     *
     *   - `required` is used to determine whether the parameter is required in the input file. Note
     *     that by default, if `required` is not set, the parameter is implicitly required. Failing
     *     to read a required parameter will result in an exception.  If successfully parsed, the
     *     container that is filled by the fully_parse() function will contain the parameter and its
     *     value. On the other hand, if a parameter is not required, it is optional and can be left
     *     out in the input file. The container that is filled by the fully_parse() function will
     *     not contain the optional parameter.
     *
     *   - `default_value` is used to provide a value that is used if the parameter is not present
     *     in the input file. If `default_value` is set, the parameter is implicitly optional. If
     *     the parameter is not present in the input file, the default value will be stored in the
     *     container that is filled by the fully_parse() function.
     *
     *   - Setting both `required = true` and a `default_value` is a logical error and will result
     *     in an exception.
     *
     * When you decide how to set the `required` and `default_value` fields, consider the following
     * cases:
     *
     *   - If you always require a parameter and there is no reasonable default value, do not set
     *     `required` or `default_value`. This will make the parameter required by default. Parsing
     *     will fail if the parameter is not present, but after parsing you can be sure that the
     *     parameter can safely be retrieved with InputParameterContainer::get(). A good example is
     *     the time step size in a time integration scheme: this parameter is always required and
     *     taking an arbitrary default value is not a good idea.
     *
     *   - Is there a reasonable default value for the parameter, which works in most situations? If
     *     yes, set `default_value` to this value (`required` implicitly is `false` then). This
     *     guarantees that you can always read the a value from the container with
     *     InputParameterContainer::get(). A good example is a parameter that activates or
     *     deactivates a feature, e.g., defaulting the EAS element technology to off might be
     *     reasonable.
     *
     *   - If the parameter is not required, but there is no reasonable default value, set
     *     `required = false`. This makes the parameter optional and the container might not contain
     *     the parameter if it is not present in the input file. Use
     *     InputParameterContainer::get_if() or InputParameterContainer::get_or() to safely retrieve
     *     the value from the container. Since this complicates retrieval of input data, this case
     *     should be used with caution. Try to formulate your input requirements in a way that one
     *     of the cases above applies. To still give an example, the present case may be useful for
     *     a damping parameter that, when present, activates damping using the provided value. This
     *     example demonstrates that the parameter has a double role: its presence activates
     *     damping and its value determines the damping strength. An often better way to selectively
     *     activate parameters can be achieved with the group() function, especially if a set of
     *     parameters is always required together.
     *
     *   - As an alternative to the last case, you could also use a Noneable type which allows you
     *     to treat the non-existence of a parameter explicitly via the "none" value. In this case,
     *     wrap the type T of the parameter in a Noneable<T> type and specify
     *     `.default_value = none<T>` (see also the example code above). After parsing, the
     *     container will be guaranteed to contain a Noneable<T> value which you can query with
     *     InputParameterContainer::get(). If the parameter is not present in the input file or set
     *     to "none", the Noneable<T> value will be empty.
     *
     * @tparam T The data type of the entry. Must be a SupportedType.
     *
     * @relatedalso InputSpec
     */
    template <SupportedType T, typename DataType = Internal::DataFor<T>>
    [[nodiscard]] InputSpec entry(std::string name, DataType&& data = {});

    /**
     * Create a callback that returns the value of the parameter with the given @p name. Such a
     * callback can, e.g., be used to determine the size of a vector parameter based on the value
     * of another parameter. Of course, the parameter with the given @p name must be read before
     * the parameter that uses this callback. Example:
     *
     * @code
     *   auto input_spec = group({
     *    entry<int>("N"),
     *    entry<std::vector<double>>("data", {.size = read_from_parameter<int>("N")}),
     *    });
     * @endcode
     *
     * @tparam T The type of the parameter of given @p name.
     *
     * @relatedalso InputSpec
     */
    template <typename T>
    [[nodiscard]] auto from_parameter(const std::string& name);

    /**
     * An entry whose value is a a selection from a list of choices. For example:
     *
     * @code
     * selection<int>("my_selection", {{"a", 1}, {"b", 2}, {"c", 3}});
     * @endcode
     *
     * The choices are given as a vector of pairs. The first element of the pair is the string
     * that is expected in the input file. The second element is the value that is stored and may be
     * any type. This function is for convenience, as you do not need to convert parsed string
     * values to another type yourself. A frequent use case is to map strings to enum constants.
     *
     * The remaining parameterization options follow the same rules as for the entry() function.
     *
     * @note If you want to store the choices as strings and not map them to another type, use the
     * other selection() function.
     *
     * @relatedalso InputSpec
     */
    template <typename T, typename DataType = Internal::DataFor<T>>
      requires(!std::same_as<T, std::string>)
    [[nodiscard]] InputSpec selection(
        std::string name, std::vector<std::pair<std::string, T>> choices, DataType data = {});


    /**
     * Like the the other selection() function, but the choices are stored as strings and not mapped
     * to another type.
     *
     * @note Although this function only works with strings, you still need to provide a type for
     * the first template parameter for consistency with the other functions.
     *
     * @relatedalso InputSpec
     */
    template <std::same_as<std::string> T, typename DataType = Internal::DataFor<T>>
    [[nodiscard]] InputSpec selection(
        std::string name, std::vector<std::string> choices, DataType data = {});

    /**
     * A group of InputSpecs. This groups one or more InputSpecs under a name. A group can be
     * required (default) or optional. If the group is optional and not present in the input, all
     * InputSpecs in the group are implicitly optional. Examples:
     *
     * @code
     * // A required group.
     * group("group",
     *    {
     *      entry<double>("a"),
     *      entry<int>("b"),
     *    });
     *
     * // An optional group. If the group is not present in the input, none of the entries are
     * // required. This is useful to require a group of parameters together.
     * group("GenAlpha",
     *   {
     *      entry<double>("alpha_f"),
     *      entry<double>("alpha_m"),
     *      entry<double>("gamma"),
     *   },
     *   {.required = false}
     *   );
     *
     * //Groups may be nested
     * group("outer",
     *  {
     *    entry<int>("a"),
     *    group("inner",
     *    {
     *      entry<double>("b"),
     *      entry<std::string>("c"),
     *    }
     *    ),
     *  });
     *
     * @endcode
     *
     * A group introduces a new scope in the input. This group scope is not only useful for
     * structuring the input file, but also to selectively activate a set of parameters, as
     * demonstrated in the second example. If an optional group is not present in the input, none of
     * the entries are required. If the group is present, all entries are required. This is often
     * exactly what you need: a group activates a feature which requires certain parameters to
     * be present.
     *
     * Whether a group has a default value is determined by the default values of its entries. If
     * all entries have default values, the group implicitly has a default value. In this case, if
     * the group is optional and not present in the input, the default values of the entries are
     * stored under the group name in the container.
     *
     * @note If you want to group multiple InputSpecs without creating a new named scope, use the
     * all_of() function.
     *
     * @relatedalso InputSpec
     */
    [[nodiscard]] InputSpec group(
        std::string name, std::vector<InputSpec> specs, GroupData data = {});

    /**
     * All of the given InputSpecs are expected, e.g.,
     *
     * @code
     * all_of({
     *   entry<int>("a"),
     *   entry<double>("b"),
     *   entry<std::string>("c"),
     *   });
     * @endcode
     *
     * will require all three entries to be present in the input.
     *
     * The main application of this function is to gather multiple InputSpecs on the same level and
     * treat them as a single InputSpec. Nesting multiple all_of() specs is possible but does not
     * have any effect on the structure of the input. The following two examples are equivalent:
     *
     * @code
     * // version 1
     * group("outer",
     * {
     *   all_of({
     *     entry<int>("a"),
     *     all_of({
     *       entry<double>("b"),
     *     }),
     *   }),
     *   entry<std::string>("c"),
     * });
     *
     * // version 2
     * group("outer",
     * {
     *   entry<int>("a"),
     *   entry<double>("b"),
     *   entry<std::string>("c"),
     * });
     * @endcode
     *
     * Note that all_of() does not changed the `required` state of its contained InputSpecs. It
     * simply tries to parse all of them and missing optional InputSpecs are not an error. In
     * practice, all_of() is essentially a group() without a name and without an associated scope in
     * the input file. An all_of() InputSpec will be required if at least one of its contained
     * InputSpecs is required. If none of the contained InputSpecs are required, the all_of()
     * InputSpec is not required.
     *
     * @relatedalso InputSpec
     */
    [[nodiscard]] InputSpec all_of(std::vector<InputSpec> specs);

    /**
     * Exactly one of the given InputSpecs is expected. For example, to require exactly one of two
     * groups:
     *
     * @code
     * one_of({
     *  group("OneStepTheta",
     *  {
     *    entry<double>("theta"),
     *  }),
     *  group("GenAlpha",
     *  {
     *    entry<double>("alpha_f"),
     *    entry<double>("alpha_m"),
     *    entry<double>("gamma"),
     *    entry<bool>("do_logging", {.default_value = false}),
     *  }),
     *  });
     * @endcode
     *
     * Here, one_of() requires either the "OneStepTheta" group or the "GenAlpha" group to be
     * present in the input. If both or none of them are present, an exception is thrown. Note that
     * all InputSpecs handed to one_of() need to be `required = true`. While this could silently be
     * changed internally, you will encounter an error if any InputSpec is not required to avoid
     * confusion and stop you from constructing difficult to understand InputSpecs. You can use
     * entries that are `required = false` nested inside other InputSpecs, see e.g. the `do_logging`
     * entry in the example code. The return one_of() InputSpec is always treated as required.
     *
     * The optional @p on_parse_callback may be used to perform additional actions after parsing one
     * of the specs. The index of the parsed spec inside the @p specs vector is passed as an
     * argument. An exemplary use case is to map the index to an enum value and store it in the
     * container. This let's you perform a switch on the enum value to easily obtain the correct
     * parsed data from the container. The store_index_as() function can be used to create such a
     * callback.
     *
     * @note The one_of() function is not intended to be used for selecting from a fixed set of
     * different values of the same type. Use the selection() function for this purpose.
     *
     * @relatedalso InputSpec
     */
    [[nodiscard]] InputSpec one_of(std::vector<InputSpec> specs,
        std::function<void(InputParameterContainer& container, std::size_t index)>
            on_parse_callback = nullptr);

    /**
     * This function may be used to produce the optional argument of the one_of() function. It
     * returns a callback that stores the index of the parsed spec in the container under @p name.
     * The index is stored on the same level as the parsed spec from the one_of() function.
     * Example:
     *
     * @code
     * enum class TimeIntegrationMethod
     * {
     *   OneStepTheta,
     *   GenAlpha,
     * };
     *
     * one_of({
     *  group("OneStepTheta",
     *  {
     *    entry<double>("theta"),
     *  }),
     *  group("GenAlpha",
     *  {
     *    entry<double>("alpha_f"),
     *    entry<double>("alpha_m"),
     *    entry<double>("gamma"),
     *  }),
     *  },
     *  store_index_as<TimeIntegrationMethod>("index",
     *    {TimeIntegrationMethod::OneStepTheta,
     *     TimeIntegrationMethod::GenAlpha})
     * );
     * @endcode
     *
     * Additionally, you can provide a reindexing vector to map the indices to other values. This is
     * especially useful if you map the often arbitrarily ordered indices to enum constants that
     * document the meaning of the index. This is demonstrated in the example above. If the
     * @p reindexing vector is not provided, the index is stored as is.
     *
     * @relatedalso InputSpec
     */
    template <typename T>
    auto store_index_as(std::string name, std::vector<T> reindexing = {});
  }  // namespace InputSpecBuilders
}  // namespace Core::IO


// --- template definitions --- //

template <typename DataType>
void Core::IO::InputSpecBuilders::Internal::BasicSpec<DataType>::parse(
    ValueParser& parser, InputParameterContainer& container) const
{
  parser.consume(name);

  if constexpr (InputSpecBuilders::Internal::is_sized_data<DataType>)
  {
    struct SizeVisitor
    {
      std::size_t operator()(std::size_t size) const { return size; }
      std::size_t operator()(const typename DataType::SizeCallback& callback) const
      {
        return callback(container);
      }
      InputParameterContainer& container;
    };

    std::size_t size = std::visit(SizeVisitor{container}, data.size);
    auto parsed = parser.read<typename DataType::StoredType>(size);
    container.add(name, parsed);
  }
  else
  {
    container.add(name, parser.read<typename DataType::StoredType>());
  }

  if (data.on_parse_callback) data.on_parse_callback(container);
}


template <typename DataTypeIn>
bool Core::IO::InputSpecBuilders::Internal::BasicSpec<DataTypeIn>::match(ConstYamlNodeRef node,
    InputParameterContainer& container, IO::Internal::MatchEntry& match_entry) const
{
  match_entry.type = IO::Internal::MatchEntry::Type::entry;
  auto spec_name = ryml::to_csubstr(name);
  if (!node.node.has_child(spec_name))
  {
    // It is OK to not encounter an optional entry
    if (data.default_value.has_value())
    {
      container.add(name, data.default_value.value());
      match_entry.state = IO::Internal::MatchEntry::State::defaulted;
      return true;
    }
    return !data.required.value();
  }

  // A child with the name of the spec exists, so this is at least a partial match.
  match_entry.state = IO::Internal::MatchEntry::State::partial;
  auto entry_node = node.wrap(node.node[spec_name]);

  FOUR_C_ASSERT(entry_node.node.key() == name, "Internal error.");

  try
  {
    typename DataTypeIn::StoredType value;
    read_value_from_yaml(entry_node, value);
    container.add(name, value);
    match_entry.state = IO::Internal::MatchEntry::State::matched;
    match_entry.matched_node = entry_node.node.id();
  }
  catch (const Core::Exception& e)
  {
    return false;
  }
  return true;
}


template <typename DataTypeIn>
void Core::IO::InputSpecBuilders::Internal::BasicSpec<DataTypeIn>::emit_metadata(
    ryml::NodeRef node) const
{
  node |= ryml::MAP;
  node["name"] << name;

  node["type"] << IO::Internal::get_pretty_type_name<StoredType>();
  if (!data.description.empty())
  {
    node["description"] << data.description;
  }
  emit_value_as_yaml(node["required"], data.required.value());
  if (data.default_value.has_value())
  {
    emit_value_as_yaml(node["default"], data.default_value.value());
  }

  if constexpr (InputSpecBuilders::Internal::is_sized_data<DataType>)
  {
    // Size can only be emitted if it is fixed.
    if (data.size.index() == 0)
    {
      node["size"] << std::get<0>(data.size);
    }
  }
}


template <typename DataTypeIn>
void Core::IO::InputSpecBuilders::Internal::SelectionSpec<DataTypeIn>::parse(
    ValueParser& parser, InputParameterContainer& container) const
{
  parser.consume(name);
  auto value = parser.read<std::string>();
  for (const auto& choice : choices)
  {
    if (choice.first == value)
    {
      container.add(name, choice.second);
      return;
    }
  }
  FOUR_C_THROW("Invalid value '%s'", value.c_str());
}


template <typename DataTypeIn>
bool Core::IO::InputSpecBuilders::Internal::SelectionSpec<DataTypeIn>::match(ConstYamlNodeRef node,
    InputParameterContainer& container, IO::Internal::MatchEntry& match_entry) const
{
  match_entry.type = IO::Internal::MatchEntry::Type::entry;
  auto spec_name = ryml::to_csubstr(name);
  if (!node.node.has_child(spec_name))
  {
    // It is OK to not encounter an optional entry
    if (data.default_value.has_value())
    {
      container.add(name, data.default_value.value());
      match_entry.state = IO::Internal::MatchEntry::State::defaulted;
      return true;
    }
    return !data.required.value();
  }

  auto entry_node = node.wrap(node.node[spec_name]);

  FOUR_C_ASSERT(entry_node.node.key() == name, "Internal error.");

  try
  {
    std::string value;
    read_value_from_yaml(entry_node, value);

    for (const auto& choice : choices)
    {
      if (choice.first == value)
      {
        container.add(name, choice.second);
        match_entry.state = IO::Internal::MatchEntry::State::matched;
        match_entry.matched_node = entry_node.node.id();
        if (data.on_parse_callback) data.on_parse_callback(container);
        return true;
      }
    }
  }
  catch (const std::exception& e)
  {
    return false;
  }

  return false;
}

template <typename DataTypeIn>
void Core::IO::InputSpecBuilders::Internal::SelectionSpec<DataTypeIn>::print(
    std::ostream& stream, std::size_t indent) const
{
  stream << "// " << std::string(indent, ' ') << name;

  if (!data.required)
  {
    stream << " (optional)";
  }
  if (data.default_value.has_value())
  {
    // Find the choice that corresponds to the default value.
    auto default_value_it = std::find_if(choices.begin(), choices.end(),
        [&](const auto& choice) { return choice.second == data.default_value.value(); });
    FOUR_C_ASSERT(
        default_value_it != choices.end(), "Internal error: default value not found in choices.");

    stream << " (default: ";
    stream << default_value_it->first;
    stream << ")";
  }
  {
    stream << " (choices: ";
    for (const auto& choice : choices)
    {
      stream << choice.first << "|";
    }
    stream << ")";
  }
  if (!data.description.empty())
  {
    stream << " " << std::quoted(Core::Utils::trim(data.description));
  }
  stream << "\n";
}

template <typename DataTypeIn>
void Core::IO::InputSpecBuilders::Internal::SelectionSpec<DataTypeIn>::emit_metadata(
    ryml::NodeRef node) const
{
  node |= ryml::MAP;
  node["name"] << name;

  node["type"] = "selection";
  if (!data.description.empty())
  {
    node["description"] << data.description;
  }
  emit_value_as_yaml(node["required"], data.required.value());
  if (data.default_value.has_value())
  {
    // Find the choice that corresponds to the default value.
    auto default_value_it = std::find_if(choices.begin(), choices.end(),
        [&](const auto& choice) { return choice.second == data.default_value.value(); });
    FOUR_C_ASSERT(
        default_value_it != choices.end(), "Internal error: default value not found in choices.");
    node["default"] << default_value_it->first;
  }
  node["choices"] |= ryml::SEQ;
  for (const auto& choice : choices)
  {
    auto entry = node["choices"].append_child();
    // Write every choice entry as a map to easily extend the information at a later point.
    entry |= ryml::MAP;
    entry["name"] << choice.first;
  }
}

template <typename T>
auto Core::IO::InputSpecBuilders::from_parameter(const std::string& name)
{
  return [name](const InputParameterContainer& container) -> T
  {
    const T* val = container.get_if<T>(name);
    FOUR_C_ASSERT_ALWAYS(val, "Parameter '%s' not found in container.", name.c_str());
    return *val;
  };
}


template <Core::IO::InputSpecBuilders::SupportedType T, typename DataType>
Core::IO::InputSpec Core::IO::InputSpecBuilders::entry(std::string name, DataType&& data)
{
  Internal::sanitize_required_default(data);
  return IO::Internal::make_spec(Internal::BasicSpec<DataType>{.name = name, .data = data},
      {
          .name = name,
          .description = data.description,
          .required = data.required.value(),
          .has_default_value = data.default_value.has_value(),
          .n_specs = 1,
      });
}


template <typename T, typename DataType>
Core::IO::InputSpec Core::IO::InputSpecBuilders::Internal::selection_internal(
    std::string name, std::vector<std::pair<std::string, T>> choices, DataType data)
{
  FOUR_C_ASSERT_ALWAYS(!choices.empty(), "Selection must have at least one choice.");
  Internal::sanitize_required_default(data);

  if (data.default_value.has_value())
  {
    auto default_value_it = std::find_if(choices.begin(), choices.end(),
        [&](const auto& choice) { return choice.second == data.default_value.value(); });

    if (default_value_it == choices.end())
    {
      std::string error_message;
      for (const auto& choice : choices)
      {
        error_message += choice.first + "|";
      }

      std::stringstream default_value_stream;
      if constexpr (Core::IO::Internal::CustomDatPrintable<T>)
      {
        Core::IO::Internal::DatPrinter{}(default_value_stream, data.default_value.value());
      }
      else
      {
        default_value_stream << "<unprintable>";
      }
      FOUR_C_THROW("Default value '%s' of selection not found in choices '%s'.",
          default_value_stream.str().c_str(), error_message.c_str());
    }
  }

  return IO::Internal::make_spec(
      Internal::SelectionSpec<DataType>{.name = name, .data = data, .choices = choices},
      {
          .name = name,
          .description = data.description,
          .required = data.required.value(),
          .has_default_value = data.default_value.has_value(),
          .n_specs = 1,
      });
}


template <typename T, typename DataType>
  requires(!std::same_as<T, std::string>)
Core::IO::InputSpec Core::IO::InputSpecBuilders::selection(
    std::string name, std::vector<std::pair<std::string, T>> choices, DataType data)
{
  return Internal::selection_internal(name, choices, data);
}


template <std::same_as<std::string> T, typename DataType>
Core::IO::InputSpec Core::IO::InputSpecBuilders::selection(
    std::string name, std::vector<std::string> choices, DataType data)
{
  std::vector<std::pair<std::string, std::string>> choices_with_strings;
  for (const auto& choice : choices)
  {
    choices_with_strings.emplace_back(choice, choice);
  }
  return Internal::selection_internal(name, choices_with_strings, data);
}



template <typename T>
auto Core::IO::InputSpecBuilders::store_index_as(std::string name, std::vector<T> reindexing)
{
  return [name, reindexing](InputParameterContainer& container, std::size_t index)
  {
    if (reindexing.empty())
    {
      container.add(name, static_cast<T>(index));
    }
    else
    {
      FOUR_C_ASSERT(index < reindexing.size(), "Index out of bounds.");
      container.add(name, reindexing[index]);
    }
  };
}

FOUR_C_NAMESPACE_CLOSE

#endif
