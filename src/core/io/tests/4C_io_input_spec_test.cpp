// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_spec_builders.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_unittest_utils_assertions_test.hpp"


namespace
{
  using namespace FourC;
  using namespace FourC::Core::IO;
  using namespace FourC::Core::IO::InputSpecBuilders;

  TEST(InputSpecTest, Simple)
  {
    auto line = all_of({
        tag("marker"),
        entry<int>("a", {.description = "An integer", .default_value = 1}),
        entry<double>("b", {.required = true}),
        entry<std::string>("c", {.required = false}),
        entry<bool>("d"),
    });
    InputParameterContainer container;
    std::string stream("marker b 2.0 d true // trailing comment");
    ValueParser parser(stream);
    line.fully_parse(parser, container);
    EXPECT_EQ(container.get<int>("a"), 1);
    EXPECT_EQ(container.get<double>("b"), 2.0);
    EXPECT_EQ(container.get_if<std::string>("c"), nullptr);
    EXPECT_EQ(container.get<bool>("d"), true);
  }

  TEST(InputSpecTest, OutofOrder)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<double>("b"),
        entry<std::string>("c"),
    });
    InputParameterContainer container;
    std::string stream("b 2.0 c string a 1");
    ValueParser parser(stream);
    line.fully_parse(parser, container);
    EXPECT_EQ(container.get<int>("a"), 1);
    EXPECT_EQ(container.get<double>("b"), 2.0);
    EXPECT_EQ(container.get<std::string>("c"), "string");
  }

  TEST(InputSpecTest, OptionalLeftOut)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<double>("b"),
        entry<std::string>("c", {.default_value = "default"}),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0 // c 1");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2.0);
      EXPECT_EQ(container.get<std::string>("c"), "default");
    }
  }

  TEST(InputSpecTest, RequiredLeftOut)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<double>("b"),
        entry<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
          "Required value 'c' not found in input line");
    }
  }

  TEST(InputSpecTest, Vector)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<std::vector<double>>("b", {.size = 3}),
        entry<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 1.0 2.0 3.0 c string");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      const auto& b = container.get<std::vector<double>>("b");
      EXPECT_EQ(b.size(), 3);
      EXPECT_EQ(b[0], 1.0);
      EXPECT_EQ(b[1], 2.0);
      EXPECT_EQ(b[2], 3.0);
      EXPECT_EQ(container.get<std::string>("c"), "string");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b 1.0 2.0 c string");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
          "Could not parse 'c' as a double value");
    }
  }

  TEST(InputSpecTest, Noneable)
  {
    auto line = all_of({
        entry<int>("size"),
        entry<std::vector<Noneable<int>>>("a", {.size = from_parameter<int>("size")}),
        entry<Noneable<std::string>>("b", {.description = "b"}),
        entry<Noneable<double>>("c", {.default_value = 1.0}),
        entry<Noneable<bool>>("d", {.required = false}),
        entry<Noneable<int>>("e", {.default_value = none<int>}),
    });

    {
      InputParameterContainer container;
      std::string stream("size 3 a 1 2 3 b none");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      const auto& a = container.get<std::vector<Noneable<int>>>("a");
      EXPECT_EQ(a.size(), 3);
      EXPECT_EQ(a[0].has_value(), true);
      EXPECT_EQ(a[0].value(), 1);
      EXPECT_EQ(a[1].has_value(), true);
      EXPECT_EQ(a[1].value(), 2);
      EXPECT_EQ(a[2].has_value(), true);
      EXPECT_EQ(a[2].value(), 3);

      const auto& b = container.get<Noneable<std::string>>("b");
      EXPECT_EQ(b.has_value(), false);

      const auto& c = container.get<Noneable<double>>("c");
      EXPECT_EQ(c.has_value(), true);
      EXPECT_EQ(c.value(), 1.0);

      const auto& e = container.get<Noneable<int>>("e");
      EXPECT_EQ(e.has_value(), false);
    }

    {
      InputParameterContainer container;
      std::string stream("size 3 a 1 none 3 b none c none d none e none");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      const auto& a = container.get<std::vector<Noneable<int>>>("a");
      EXPECT_EQ(a.size(), 3);
      EXPECT_EQ(a[0].has_value(), true);
      EXPECT_EQ(a[0].value(), 1);
      EXPECT_EQ(a[1].has_value(), false);
      EXPECT_EQ(a[2].has_value(), true);
      EXPECT_EQ(a[2].value(), 3);

      const auto& b = container.get<Noneable<std::string>>("b");
      EXPECT_EQ(b.has_value(), false);

      const auto& c = container.get<Noneable<double>>("c");
      EXPECT_EQ(c.has_value(), false);

      const auto& d = container.get<Noneable<bool>>("d");
      EXPECT_EQ(d.has_value(), false);

      const auto& e = container.get<Noneable<int>>("e");
      EXPECT_EQ(e.has_value(), false);
    }

    {
      InputParameterContainer container;
      std::string stream("size 3 a 1 2 3 b string c 2.0 d true e 42");
      ValueParser parser(stream);
      line.fully_parse(parser, container);

      const auto& b = container.get<Noneable<std::string>>("b");
      EXPECT_EQ(b.has_value(), true);
      EXPECT_EQ(b.value(), "string");

      const auto& c = container.get<Noneable<double>>("c");
      EXPECT_EQ(c.has_value(), true);
      EXPECT_EQ(c.value(), 2.0);

      const auto& d = container.get<Noneable<bool>>("d");
      EXPECT_EQ(d.has_value(), true);
      EXPECT_EQ(d.value(), true);

      const auto& e = container.get<Noneable<int>>("e");
      EXPECT_EQ(e.has_value(), true);
      EXPECT_EQ(e.value(), 42);
    }
  }

  TEST(InputSpecTest, VectorWithParsedLength)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<std::vector<double>>("b", {.size = from_parameter<int>("a")}),
        entry<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 3 b 1.0 2.0 3.0 c string");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 3);
      const auto& b = container.get<std::vector<double>>("b");
      EXPECT_EQ(b.size(), 3);
      EXPECT_EQ(b[0], 1.0);
      EXPECT_EQ(b[1], 2.0);
      EXPECT_EQ(b[2], 3.0);
      EXPECT_EQ(container.get<std::string>("c"), "string");
    }

    {
      InputParameterContainer container;
      std::string stream("a 3 b 1.0 2.0 c string");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
          "Could not parse 'c' as a double value");
    }
  }

  TEST(InputSpecTest, UserDefined)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<double>("b"),
        user_defined<std::string>(
            "c", {.description = "A string", .default_value = "Not found"},
            [name = "c"](ValueParser& parser, InputParameterContainer& container)
            {
              parser.consume(name);
              // Some special parsing logic
              parser.consume("_");
              parser.consume("_");

              container.add<std::string>("c", "I found c");
            },
            [](std::ostream& out, std::size_t) { out << "_ _ \n"; }),
        entry<std::string>("s"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0 c _ _ s hello");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2.0);
      EXPECT_EQ(container.get<std::string>("c"), "I found c");
      EXPECT_EQ(container.get<std::string>("s"), "hello");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0 c _ s hello");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          line.fully_parse(parser, container), Core::Exception, "expected string '_'");
    }
  }

  TEST(InputSpecTest, Selection)
  {
    auto line = all_of({
        entry<int>("a"),
        selection<int>("b", {{"b1", 1}, {"b2", 2}}, {.default_value = 1}),
        selection<std::string>("c", {"c1", "c2"}, {.default_value = "c2"}),
        entry<std::string>("d"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b b2 d string c c1");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<int>("b"), 2);
      EXPECT_EQ(container.get<std::string>("c"), "c1");
      EXPECT_EQ(container.get<std::string>("d"), "string");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b b4 c string");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          line.fully_parse(parser, container), Core::Exception, "Invalid value 'b4'");
    }
  }

  TEST(InputSpecTest, Unparsed)
  {
    auto line = all_of({
        entry<int>("a"),
        entry<int>("optional", {.default_value = 42}),
        entry<double>("b"),
        entry<std::string>("c"),
    });
    InputParameterContainer container;
    std::string stream("a 1 b 2.0 c string unparsed unparsed");
    ValueParser parser(stream);
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
        "line still contains 'unparsed unparsed'");
  }


  TEST(InputSpecTest, Groups)
  {
    auto line = all_of({
        entry<int>("a"),
        group("group1",
            {
                entry<double>("b"),
            }),
        group("group2",
            {
                entry<double>("b", {.default_value = 3.0}),
                entry<std::string>("c"),
            },
            {
                .required = false,
            }),
        group("group3",
            {
                entry<std::string>("c", {.default_value = "default"}),
            },
            {
                .required = false,
            }),
        entry<std::string>("c"),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 group1 b 2.0 c string");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_EQ(const_container.group("group1").get<double>("b"), 2.0);
      EXPECT_ANY_THROW([[maybe_unused]] const auto& c = const_container.group("group2"));
      // Group 3 only contains entries that have default values, so it implicitly has a default
      // value.
      EXPECT_EQ(const_container.group("group3").get<std::string>("c"), "default");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group2 b 2.0 c string group1 b 4.0 c string");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_EQ(const_container.group("group2").get<double>("b"), 2.0);
      EXPECT_EQ(const_container.group("group1").get<double>("b"), 4.0);
      EXPECT_EQ(const_container.get<std::string>("c"), "string");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group1 b 4.0 c string");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_ANY_THROW([[maybe_unused]] const auto& c = const_container.group("group2"));
      EXPECT_EQ(const_container.group("group1").get<double>("b"), 4.0);
      EXPECT_EQ(const_container.get<std::string>("c"), "string");
    }
  }

  TEST(InputSpecTest, NestedAllOf)
  {
    auto line = all_of({
        entry<int>("a"),
        all_of({
            all_of({
                entry<double>("b"),
            }),
        }),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2.0");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      const auto& const_container = container;
      EXPECT_EQ(const_container.get<int>("a"), 1);
      EXPECT_EQ(const_container.get<double>("b"), 2.0);
    }
  }

  TEST(InputSpecTest, OneOf)
  {
    auto line = all_of({
        entry<int>("a", {.default_value = 42}),
        one_of({
            entry<double>("b"),
            group("group",
                {
                    entry<std::string>("c"),
                    entry<double>("d"),
                }),
        }),
    });

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2);
    }

    {
      InputParameterContainer container;
      std::string stream("group c string d 2.0");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 42);
      EXPECT_EQ(container.group("group").get<std::string>("c"), "string");
      EXPECT_EQ(container.group("group").get<double>("d"), 2);
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group c string d 2.0 b 3.0");
      ValueParser parser(stream);
      // More than one of the one_of entries is present. Refuse to parse any of them.
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
          "still contains 'group c string d 2.0 b 3.0'");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 group c string");
      // Note: we start to parse the group, but the entries are not complete, so we backtrack.
      // The result is that the parts of the group remain unparsed.
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          line.fully_parse(parser, container), Core::Exception, "still contains 'group c string'");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(
          line.fully_parse(parser, container), Core::Exception, "Required 'one_of {b, group}'");
    }
  }

  TEST(InputSpecTest, OneOfTopLevel)
  {
    auto line = one_of(
        {
            all_of({
                entry<int>("a"),
                entry<double>("b"),
            }),
            all_of({
                entry<std::string>("c"),
                entry<double>("d"),
            }),
        },
        // Additionally store the index of the parsed group but map it to a different value.
        store_index_as<int>("index", /*reindex*/ {1, 10}));

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<double>("b"), 2);
      EXPECT_EQ(container.get<int>("index"), 1);
    }

    {
      InputParameterContainer container;
      std::string stream("c string d 2.0");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<std::string>("c"), "string");
      EXPECT_EQ(container.get<double>("d"), 2);
      EXPECT_EQ(container.get<int>("index"), 10);
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b 2 c string d 2.0");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
          "both 'all_of {a, b}' and 'all_of {c, d}'");
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 c string");
      ValueParser parser(stream);
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(line.fully_parse(parser, container), Core::Exception,
          "one_of {all_of {a, b}, all_of {c, d}}");
    }
  }

  TEST(InputSpecTest, PrintAsDat)
  {
    auto line =
        group("g", {
                       // Note: the all_of entries will be pulled into the parent group.
                       all_of({
                           entry<int>("a", {.description = "An integer"}),
                           tag("b", {.description = "A tag", .required = false}),
                           selection<int>("c", {{"c1", 1}, {"c2", 2}},
                               {.description = "Selection", .default_value = 1}),
                       }),
                       entry<int>("d", {.description = "Another integer", .default_value = 42}),
                   });

    {
      std::ostringstream out;
      line.print_as_dat(out);
      EXPECT_EQ(out.str(), R"(// g:
//   a <int> "An integer"
//   b <tag> (optional) "A tag"
//   c (default: c1) (choices: c1|c2|) "Selection"
//   d <int> (optional) (default: 42) "Another integer"
)");
    }
  }

  TEST(InputSpecTest, EmitMetadata)
  {
    auto line = all_of({
        entry<int>("a", {.default_value = 42}),
        entry<std::vector<double>>("b", {.default_value = std::vector{1., 2., 3.}, .size = 3}),
        one_of({
            all_of({
                entry<std::pair<int, std::string>>("b", {.default_value = std::pair{1, "abc"}}),
                entry<std::string>("c"),
            }),
            group("group",
                {
                    entry<std::string>("c"),
                    entry<double>("d"),
                }),
        }),
        selection<int>("e", {{"e1", 1}, {"e2", 2}}, {.default_value = 1}),
    });


    {
      std::ostringstream out;
      ryml::Tree tree;
      ryml::NodeRef root = tree.rootref();
      YamlEmitter yaml{root};
      line.emit_metadata(yaml);
      out << tree;

      std::string expected = R"('all_of {a, b, one_of {all_of {b, c}, group}, e}':
  type: all_of
  description: 'all_of {a, b, one_of {all_of {b, c}, group}, e}'
  required: true
  specs:
    a:
      type: int
      description: ''
      required: false
      default: 42
    b:
      type: vector<double>
      description: ''
      required: false
      default: [1,2,3]
    'one_of {all_of {b, c}, group}':
      type: one_of
      description: 'one_of {all_of {b, c}, group}'
      required: true
      specs:
        'all_of {b, c}':
          type: all_of
          description: 'all_of {b, c}'
          required: true
          specs:
            b:
              type: 'pair<int, string>'
              description: ''
              required: false
              default: [1,abc]
            c:
              type: string
              description: ''
              required: true
        group:
          type: group
          description: ''
          required: true
          specs:
            c:
              type: string
              description: ''
              required: true
            d:
              type: double
              description: ''
              required: true
    e:
      type: selection
      description: ''
      required: false
      default: 1
      choices:
        e1: 1
        e2: 2
)";
      EXPECT_EQ(out.str(), expected);
    }
  }

  TEST(InputSpecTest, Copyable)
  {
    InputSpec line;
    {
      auto tmp = all_of({
          entry<int>("a"),
          entry<std::string>("b"),
          selection<int>("c", {{"c1", 1}, {"c2", 2}}, {.default_value = 1}),
      });

      line = all_of({
          tmp,
          entry<int>("d"),
      });
    }

    {
      InputParameterContainer container;
      std::string stream("a 1 b string c c2 d 42");
      ValueParser parser(stream);
      line.fully_parse(parser, container);
      EXPECT_EQ(container.get<int>("a"), 1);
      EXPECT_EQ(container.get<std::string>("b"), "string");
      EXPECT_EQ(container.get<int>("c"), 2);
      EXPECT_EQ(container.get<int>("d"), 42);
    }
  }
}  // namespace
