// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_symbolic_expression.hpp"

#include "4C_utils_enum.hpp"
#include "4C_utils_exceptions.hpp"

#include <Sacado.hpp>

#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <std::size_t n>
  [[nodiscard]] bool contains(const std::array<std::string, n>& array, const std::string& element)
  {
    return std::find(array.begin(), array.end(), element) != array.end();
  }
}  // namespace

namespace Core::Utils::SymbolicExpressionDetails
{
  using IndexType = std::size_t;
  constexpr IndexType invalid_index = -1uL;

  enum class NodeType : unsigned char
  {
    // independent variable: 't', 'x', ...
    variable,
    // numeric value
    number,
    // mathematical function e.g. 'sin', 'cos', 'exp', ...
    function,
    // (binary) operator e.g. '+', '-', '*', '/', ...
    op,
  };

  /*!
  \brief Syntax tree node holding binary, unary operator or literals
  */
  struct SyntaxTreeNode
  {
    /**
     * Node either holds a number, or string operator, variable name or a function name.
     */
    union
    {
      double number{};
      std::string_view str;  // this views into the expression string or a string literal
    } as;

    //! Index of the node in the syntax tree
    IndexType my_index;

    //! Lhs and rhs nodes
    IndexType lhs_index{-1u};
    IndexType rhs_index{-1u};

    //! type of the node, ie operator, literal, ...
    NodeType type;

    [[nodiscard]] const double& number() const
    {
      FOUR_C_ASSERT(type == NodeType::number, "Node is not a number.");
      return as.number;
    }

    [[nodiscard]] const std::string_view& str() const
    {
      FOUR_C_ASSERT(
          type == NodeType::variable || type == NodeType::function || type == NodeType::op,
          "Node is not a variable or function.");
      return as.str;
    }
  };


  /*----------------------------------------------------------------------*/
  /*!
  \brief Class holds auxiliary variables for Lexan method which steps through
         the string destilling the function tokens
  */
  class Lexer
  {
   public:
    //! constructor
    Lexer(std::string_view funct) : funct_(funct), pos_(0) {}

    //! delivers funct_ character at position pos_++
    int get_next();

    //! identifies a token (value and kind) in the funct_ string
    void lexan();

    //! type of identifiable tokens in string funct_
    enum TokenType
    {
      tok_none,
      tok_done,
      tok_name,  // operator name, e.g. 'sin'
      tok_int,   // integer number
      tok_real,  // reals number
      tok_add,   // addition '+'
      tok_sub,   // subtraction and negation '-'
      tok_mul,   // multiplication '*'
      tok_div,   // division '/'
      tok_pow,   // power '^'
      tok_lpar,  // left parenthesis '('
      tok_rpar,  // right parenthesis ')'
      tok_comma  // comma ',' (used to separate function arguments)
    };

    std::string_view funct_;  // function description as string, i.e. "t^2", "sin(t)", etc
    unsigned pos_;            // current position in string funct_
    TokenType tok_;           // current token of string funct_
    const char* str_;         // pointer to begin of current string token in funct_
    int integer_;             // translated integer number or length of operator word
    double real_;             // translated real number
  };

  //! helper structs
  template <typename T>
  struct IsFAD : public std::false_type
  {
  };

  template <typename T>
  struct IsFAD<Sacado::Fad::DFad<T>> : public std::true_type
  {
  };


  /*----------------------------------------------------------------------*/
  /*!
  \brief Parser
  */
  template <class T>
  class Parser
  {
   public:
    using NodePtr = SyntaxTreeNode*;

    //! constructor
    Parser(std::string funct);

    //! destructor
    ~Parser() = default;

    //! copy constructor
    Parser(const Parser& other) : symbolicexpression_(other.symbolicexpression_) { parse(); }

    //! copy assignment operator
    Parser& operator=(const Parser& other)
    {
      if (&other == this) return *this;
      symbolicexpression_ = other.symbolicexpression_;
      node_arena_.clear();
      parsed_variable_constant_names_.clear();
      variable_values_ = nullptr;
      constant_values_ = nullptr;

      parse();
      return *this;
    };

    //! move constructor
    Parser(Parser&& other) noexcept = default;

    //! move assignment operator
    Parser& operator=(Parser&& other) noexcept = default;

    /*!
     * @brief evaluates the parsed expression for a given set of variables
     *
     * @param[in] variable_values A map containing all variables (variablename, value) necessary to
     * evaluate the parsed expression
     * @return Value of the parsed expression
     */
    T evaluate_expression(const std::map<std::string, T>& variable_values) const;

    /*!
     * @brief evaluates the derivative of the parsed expression with respect to a given set of
     * variables
     *
     * @param[in] variable_values A map containing all variables (variablename, value) necessary to
     * evaluate the parsed expression. Since the derivative of the parsed expression is evaluated,
     * only  Sacado::Fad::DFad<T> types are allowed
     * @param[in] constants A map containing all constants (constantname, value) necessary
     * to evaluate the parsed expression
     * @return  Derivative of the parsed expression with respect to the variables
     */
    template <typename T2>
      requires IsFAD<T2>::value
    T evaluate_derivative(const std::map<std::string, T2>& variable_values,
        const std::map<std::string, double>& constants = {}) const;

    //! Check if a variable with name 'varname' exists
    [[nodiscard]] bool is_variable(const std::string& varname) const;

   private:
    void parse();
    NodePtr parse_primary(Lexer& lexer);
    NodePtr parse_pow(Lexer& lexer);
    NodePtr parse_term(Lexer& lexer);
    NodePtr parse_expr(Lexer& lexer);
    NodePtr parse(Lexer& lexer);

    //! Create a new node in the node arena and return a pointer to it.
    //! Optionally, a left-hand side and right-hand side node can be provided.
    NodePtr create_node(NodeType type, NodePtr lhs = nullptr, NodePtr rhs = nullptr);

    //! given symbolic expression
    std::string symbolicexpression_;

    //! evaluates the parsed expression
    T evaluate(const std::map<std::string, T>& variable_values,
        const std::map<std::string, double>& constants = {}) const;

    //! recursively extract corresponding number out of a syntax tree node
    T interpret(const SyntaxTreeNode& node) const;

    const SyntaxTreeNode& at(IndexType index) const
    {
      FOUR_C_ASSERT(index < node_arena_.size(), "Index out of bounds: {}", index);
      return node_arena_[index];
    }

    /**
     * Actual storage for all parsed syntax tree nodes. The first element is the root node of
     * the tree.
     */
    std::vector<SyntaxTreeNode> node_arena_;

    //! root node of the syntax tree
    SyntaxTreeNode* root_{nullptr};

    //! set of all parsed variables
    std::set<std::string_view> parsed_variable_constant_names_;

    //! necessary variable values for evaluating the parsed expression
    mutable const std::map<std::string, T>* variable_values_{nullptr};

    //! necessary constant values for evaluating the parsed expression
    mutable const std::map<std::string, double>* constant_values_{nullptr};
  };


  /*======================================================================*/
  /* Lexer methods */

  /*----------------------------------------------------------------------*/
  /*!
  \brief method used to step through std::string funct_
         delivers its character at position pos_++
  */
  int Lexer::get_next()
  {
    if (pos_ < funct_.length())
    {
      return funct_[pos_++];
    }
    else
    {
      return EOF;
    }
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Identify current token
         type: tok_,
         value: integer_, real_,
         operator name: str_
  */
  void Lexer::lexan()
  {
    for (;;)
    {
      int t = get_next();
      if ((t == ' ') || (t == '\t'))
      {
        /* ignore whitespaces */
      }
      else if (t == '\n')
      {
        FOUR_C_THROW("newline in function definition");
      }
      else if (t == EOF)
      {
        tok_ = Lexer::Lexer::tok_done;
        return;
      }
      else
      {
        if (isdigit(t))
        {
          str_ = &(funct_[pos_ - 1]);
          while (isdigit(t))
          {
            t = get_next();
          }
          if ((t != '.') && (t != 'E') && (t != 'e'))
          {
            if (t != EOF)
            {
              pos_--;
            }
            integer_ = atoi(str_);
            tok_ = Lexer::tok_int;
            return;
          }
          if (t == '.')
          {
            t = get_next();
            if (isdigit(t))
            {
              while (isdigit(t))
              {
                t = get_next();
              }
            }
            else
            {
              FOUR_C_THROW("no digits after point at pos {}", pos_);
            }
          }
          if ((t == 'E') || (t == 'e'))
          {
            t = get_next();
            if ((t == '-') || (t == '+'))
            {
              t = get_next();
            }
            if (isdigit(t))
            {
              while (isdigit(t))
              {
                t = get_next();
              }
            }
            else
            {
              FOUR_C_THROW("no digits after exponent at pos {}", pos_);
            }
          }
          if (t != EOF)
          {
            pos_--;
          }
          real_ = strtod(str_, nullptr);
          tok_ = Lexer::tok_real;
          return;
        }
        else if (isalpha(t) || (t == '_'))
        {
          str_ = &(funct_[pos_ - 1]);
          while (isalnum(t) || (t == '_'))
          {
            t = get_next();
          }
          if (t != EOF)
          {
            pos_--;
          }
          tok_ = Lexer::tok_name;
          integer_ = &(funct_[pos_]) - str_;  // length of operator name, e.g. 'sin' has '3'
          return;
        }
        else if (t == '+')
        {
          tok_ = Lexer::tok_add;
          return;
        }
        else if (t == '-')
        {
          tok_ = Lexer::tok_sub;
          return;
        }
        else if (t == '*')
        {
          tok_ = Lexer::tok_mul;
          return;
        }
        else if (t == '/')
        {
          tok_ = Lexer::tok_div;
          return;
        }
        else if (t == '^')
        {
          tok_ = Lexer::tok_pow;
          return;
        }
        else if (t == '(')
        {
          tok_ = Lexer::tok_lpar;
          return;
        }
        else if (t == ')')
        {
          tok_ = Lexer::tok_rpar;
          return;
        }
        else if (t == ',')
        {
          tok_ = Lexer::tok_comma;
          return;
        }
        else
        {
          if (t >= 32)
            FOUR_C_THROW("unexpected char '{}' at pos {}", t, pos_);
          else
            FOUR_C_THROW("unexpected char '{}' at pos {}", t, pos_);
          tok_ = Lexer::tok_none;
          return;
        }
      }
    }
  }

  /*===================================== da=================================*/
  /* Parser methods */

  /*----------------------------------------------------------------------*/
  /*!
  \brief Constructor of parser object
  */
  template <class T>
  Parser<T>::Parser(std::string funct) : symbolicexpression_{std::move(funct)}
  {
    parse();
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Check if input is a variable name
  */
  template <class T>
  bool Parser<T>::is_variable(const std::string& varname) const
  {
    //! check if variable exists
    return (parsed_variable_constant_names_.count(varname) != 0);
  }

  template <class T>
  void Parser<T>::parse()
  {
    //! create Lexer which stores all token of funct
    Lexer lexer{symbolicexpression_};

    //! retrieve first token of funct
    lexer.lexan();

    //! create syntax tree equivalent to funct
    root_ = parse(lexer);
  }


  template <class T>
  T Parser<T>::evaluate_expression(const std::map<std::string, T>& variable_values) const
  {
    return evaluate(variable_values);
  }


  template <class T>
  template <typename T2>
    requires IsFAD<T2>::value
  T Parser<T>::evaluate_derivative(const std::map<std::string, T2>& variable_values,
      const std::map<std::string, double>& constants) const
  {
    return evaluate(variable_values, constants);
  }


  template <class T>
  T Parser<T>::evaluate(const std::map<std::string, T>& variable_values,
      const std::map<std::string, double>& constants) const
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    const bool all_required_variables_passed =
        std::all_of(parsed_variable_constant_names_.begin(), parsed_variable_constant_names_.end(),
            [&](const auto& var_name)
            {
              return (variable_values.count(std::string(var_name)) +
                         constants.count(std::string(var_name))) == 1;
            });

    if (!all_required_variables_passed)
    {
      std::string evaluate_variable_names = std::accumulate(variable_values.begin(),
          variable_values.end(), std::string(), [](const std::string& acc, const auto& v)
          { return acc.empty() ? v.first : acc + ", " + v.first; });

      std::string evaluate_constant_names = std::accumulate(constants.begin(), constants.end(),
          std::string(), [](const std::string& acc, const auto& v)
          { return acc.empty() ? v.first : acc + ", " + v.first; });

      FOUR_C_THROW(
          "Some variables that this parser encountered in the expression are not passed to "
          "the Evaluate function.\n\n"
          "Expression:  {} \n"
          "Variables passed to Evaluate: {} \n"
          "Constants passed to Evaluate: {}",
          symbolicexpression_.c_str(), evaluate_variable_names.c_str(),
          evaluate_constant_names.c_str());
    }
#endif

    //! safety check if variable_values_ and constant_values_ are nullptr
    FOUR_C_ASSERT(variable_values_ == nullptr, "Internal error");
    FOUR_C_ASSERT(constant_values_ == nullptr, "Internal error");

    //! set variable values
    variable_values_ = &variable_values;
    //! set constant values if map of constants is not empty
    if (constants.size() != 0) constant_values_ = &constants;

    //! check if function has been parsed
    FOUR_C_ASSERT(!node_arena_.empty(), "Internal error");

    //! evaluate syntax tree of function depending on set variables
    auto result = this->interpret(*root_);

    //! set variable_values_ and constant_values_ as nullptr again for safety reasons
    variable_values_ = nullptr;
    constant_values_ = nullptr;

    return result;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse primary entities, i.e. literals and unary operators,
         such as numbers, parentheses, independent variables, operator names
  */
  template <class T>
  auto Parser<T>::parse_primary(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs = nullptr;

    switch (lexer.tok_)
    {
      case Lexer::tok_lpar:
        lexer.lexan();
        lhs = parse_expr(lexer);
        if (lexer.tok_ != Lexer::tok_rpar) FOUR_C_THROW("')' expected");
        lexer.lexan();
        break;
      case Lexer::tok_int:
        lhs = create_node(NodeType::number);
        lhs->as.number = lexer.integer_;
        lexer.lexan();
        break;
      case Lexer::tok_real:
        lhs = create_node(NodeType::number);
        lhs->as.number = lexer.real_;
        lexer.lexan();
        break;
      case Lexer::tok_sub:
      {
        NodePtr rhs;
        lexer.lexan();
        rhs = parse_pow(lexer);
        // This is a unary minus operator.
        if (rhs->type == NodeType::number)
        {
          rhs->as.number *= -1;
          lhs = rhs;
        }
        else
        {
          lhs = create_node(NodeType::number);
          lhs->as.number = -1;
          lhs = create_node(NodeType::op, lhs, rhs);
          lhs->as.str = "*";
        }
        break;
      }
      case Lexer::tok_name:
      {
        // get substring starting from str_ with length of lexer.integer_
        std::string_view name(lexer.str_, lexer.integer_);
        if (name == "pi")
        {
          lhs = create_node(NodeType::number);
          lhs->as.number = M_PI;
          lexer.lexan();
          break;
        }
        else
        {
          if (name == "acos" or name == "asin" or name == "atan" or name == "cos" or
              name == "sin" or name == "tan" or name == "cosh" or name == "sinh" or
              name == "tanh" or name == "exp" or name == "log" or name == "log10" or
              name == "sqrt" or name == "ceil" or name == "heaviside" or name == "fabs" or
              name == "floor")
          {
            // Consume opening parenthesis
            lexer.lexan();
            if (lexer.tok_ != Lexer::tok_lpar)
              FOUR_C_THROW("'(' expected after function name '{}'", name);

            lexer.lexan();
            lhs = create_node(NodeType::function, parse_expr(lexer));
            lhs->as.str = name;

            // Consume closing parenthesis
            if (lexer.tok_ != Lexer::tok_rpar) FOUR_C_THROW("')' expected");
            lexer.lexan();
            break;
          }
          else if (name == "atan2")
          {
            // Consume opening parenthesis
            lexer.lexan();
            if (lexer.tok_ != Lexer::tok_lpar)
              FOUR_C_THROW("'(' expected after function name '{}'", name);

            lexer.lexan();
            auto first_arg = parse_expr(lexer);

            // Consume comma separating the two arguments
            if (lexer.tok_ != Lexer::tok_comma) FOUR_C_THROW("',' expected");
            lexer.lexan();
            auto second_arg = parse_expr(lexer);

            lhs = create_node(NodeType::function, first_arg, second_arg);
            lhs->as.str = name;

            // Consume closing parenthesis
            if (lexer.tok_ != Lexer::tok_rpar)
              FOUR_C_THROW("')' expected for function name '{}'", name);
            lexer.lexan();
            break;
          }
          else
          {
            lhs = create_node(NodeType::variable);
            lhs->as.str = name;
            parsed_variable_constant_names_.insert(name);
            lexer.lexan();
          }
        }
        break;
      }
      default:
        FOUR_C_THROW("unexpected token {}", lexer.tok_);
        break;
    }

    return lhs;
  }


  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by power: a^b
  */
  template <class T>
  auto Parser<T>::parse_pow(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;
    NodePtr rhs;

    lhs = parse_primary(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_pow)
      {
        lexer.lexan();
        rhs = parse_primary(lexer);
        lhs = create_node(NodeType::op, lhs, rhs);
        lhs->as.str = "^";
      }
      else
      {
        break;
      }
    }

    return lhs;
  }


  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by multiplication or division: a*b, a/b
  */
  template <class T>
  auto Parser<T>::parse_term(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;
    NodePtr rhs;

    lhs = parse_pow(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_mul)
      {
        lexer.lexan();
        rhs = parse_pow(lexer);
        lhs = create_node(NodeType::op, lhs, rhs);
        lhs->as.str = "*";
      }
      else if (lexer.tok_ == Lexer::tok_div)
      {
        lexer.lexan();
        rhs = parse_pow(lexer);
        lhs = create_node(NodeType::op, lhs, rhs);
        lhs->as.str = "/";
      }
      else
      {
        break;
      }
    }

    return lhs;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entity
  */
  template <class T>
  auto Parser<T>::parse(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;

    lhs = parse_expr(lexer);

    // check if parsing ended before processing the entire string
    if (lexer.tok_ != Lexer::tok_done)
    {
      FOUR_C_THROW("Invalid syntax: The remaining string '{}' is not parsed.",
          lexer.funct_.data() + lexer.pos_);
    }

    return lhs;
  }

  template <class T>
  typename Parser<T>::NodePtr Parser<T>::create_node(NodeType type, NodePtr lhs, NodePtr rhs)
  {
    SyntaxTreeNode& node = node_arena_.emplace_back();
    node.my_index = node_arena_.size() - 1;
    node.type = type;
    node.lhs_index = lhs ? lhs->my_index : invalid_index;
    node.rhs_index = rhs ? rhs->my_index : invalid_index;
    return &node;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by addition or subtraction: a+b, a-b
  */
  template <class T>
  auto Parser<T>::parse_expr(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;
    NodePtr rhs;

    lhs = parse_term(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_add)
      {
        lexer.lexan();
        rhs = parse_term(lexer);
        lhs = create_node(NodeType::op, lhs, rhs);
        lhs->as.str = "+";
      }
      else if (lexer.tok_ == Lexer::tok_sub)
      {
        lexer.lexan();
        rhs = parse_term(lexer);
        lhs = create_node(NodeType::op, lhs, rhs);
        lhs->as.str = "-";
      }
      else
      {
        break;
      }
    }

    return lhs;
  }


  /*----------------------------------------------------------------------*/
  /*!
  \brief Recursively extract corresponding number out of a syntax tree node
  */
  template <class T>
  T Parser<T>::interpret(const SyntaxTreeNode& node) const
  {
    T res = 0;  // the result

    switch (node.type)
    {
      // literal numbers: leaf of syntax tree node
      case NodeType::number:
        res = node.as.number;
        break;
      // binary operators: bifurcating branch of syntax tree node
      case NodeType::op:
      {
        T lhs;
        T rhs;

        // recursively visit branches and obtain sub-results
        lhs = interpret(at(node.lhs_index));
        rhs = interpret(at(node.rhs_index));

        // evaluate the node operator
        switch (node.as.str.front())
        {
          case '+':
            res = lhs + rhs;
            break;
          case '-':
            res = lhs - rhs;
            break;
          case '*':
            res = lhs * rhs;
            break;
          case '/':
            /* check for rhs==0.0? */
            res = lhs / rhs;
            break;
          case '^':
            res = std::pow(lhs, rhs);
            break;
          default:
            FOUR_C_THROW("unsupported operator '{}'", node.as.str);
        }
        break;
      }
      // independent variables: as set by user
      case NodeType::variable:
      {
        const auto& var = node.str();
        if (parsed_variable_constant_names_.count(var))
        {
          if (constant_values_ == nullptr)
          {
            if (!variable_values_->contains(std::string(var)))
            {
              FOUR_C_THROW("variable or constant '{}' not given as input in evaluate()", var);
            }
            else
            {
              res = variable_values_->at(std::string(var));
            }
          }
          else
          {
            if (variable_values_->contains(std::string(var)))
            {
              res = variable_values_->at(std::string(var));
            }
            else if (constant_values_->contains(std::string(var)))
            {
              res = constant_values_->at(std::string(var));
            }
            else
              FOUR_C_THROW("variable or constant '{}' not given as input in EvaluateDeriv()", var);
          }
        }
        else
          FOUR_C_THROW("unknown variable '{}'", var);
        break;
      }
      // unary operators
      case NodeType::function:
      {
        T arg;
        arg = interpret(at(node.lhs_index));
        const auto& function = node.str();
        if (function == "acos")
          res = acos(arg);
        else if (function == "asin")
          res = asin(arg);
        else if (function == "atan")
          res = atan(arg);
        else if (function == "cos")
          res = cos(arg);
        else if (function == "sin")
          res = sin(arg);
        else if (function == "tan")
          res = tan(arg);
        else if (function == "cosh")
          res = cosh(arg);
        else if (function == "sinh")
          res = sinh(arg);
        else if (function == "tanh")
          res = tanh(arg);
        else if (function == "exp")
          res = exp(arg);
        else if (function == "log")
          res = log(arg);
        else if (function == "log10")
          res = log10(arg);
        else if (function == "sqrt")
          res = sqrt(arg);
        else if (function == "atan2")
        {
          T arg2;
          // recursively visit branches and obtain sub-results
          arg2 = interpret(at(node.rhs_index));
          res = atan2(arg, arg2);
        }
        else if (function == "fabs")
          res = fabs(arg);
        else if (function == "heaviside")
        {
          if (arg > 0)
          {
            res = 1.0;
          }
          else
          {
            res = 0.0;
          }
        }
        else
          FOUR_C_THROW("unknown function_ '{}'", function);
        break;
      }
      default:
        FOUR_C_THROW("unknown syntax tree node type");
        break;
    }

    return res;
  }

}  // namespace Core::Utils::SymbolicExpressionDetails

template <typename T>
Core::Utils::SymbolicExpression<T>::SymbolicExpression(const std::string& expression)
    : parser_for_value_(
          std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<ValueType>>(expression)),
      parser_for_firstderivative_(
          std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<FirstDerivativeType>>(
              expression)),
      parser_for_secondderivative_(
          std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<SecondDerivativeType>>(
              expression))
{
}



template <typename T>
auto Core::Utils::SymbolicExpression<T>::value(
    const std::map<std::string, ValueType>& variable_values) const -> ValueType
{
  return parser_for_value_->evaluate_expression(variable_values);
}


template <typename T>
auto Core::Utils::SymbolicExpression<T>::first_derivative(
    std::map<std::string, FirstDerivativeType> variable_values,
    const std::map<std::string, ValueType>& constant_values) const -> FirstDerivativeType
{
  return parser_for_firstderivative_->evaluate_derivative(variable_values, constant_values);
}


template <typename T>
auto Core::Utils::SymbolicExpression<T>::second_derivative(
    const std::map<std::string, SecondDerivativeType>& variable_values,
    const std::map<std::string, ValueType>& constant_values) const -> SecondDerivativeType
{
  return parser_for_secondderivative_->evaluate_derivative(variable_values, constant_values);
}


template <typename Number>
Core::Utils::SymbolicExpression<Number>::SymbolicExpression(
    const Core::Utils::SymbolicExpression<Number>& other)
    : parser_for_value_{std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<ValueType>>(
          *other.parser_for_value_)},
      parser_for_firstderivative_{
          std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<FirstDerivativeType>>(
              *other.parser_for_firstderivative_)},
      parser_for_secondderivative_{
          std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<SecondDerivativeType>>(
              *other.parser_for_secondderivative_)}
{
}


template <typename Number>
Core::Utils::SymbolicExpression<Number>& Core::Utils::SymbolicExpression<Number>::operator=(
    const Core::Utils::SymbolicExpression<Number>& other)
{
  parser_for_value_ = std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<ValueType>>(
      *other.parser_for_value_);
  parser_for_firstderivative_ =
      std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<FirstDerivativeType>>(
          *other.parser_for_firstderivative_);
  parser_for_secondderivative_ =
      std::make_unique<Core::Utils::SymbolicExpressionDetails::Parser<SecondDerivativeType>>(
          *other.parser_for_secondderivative_);
  return *this;
}



template <typename Number>
Core::Utils::SymbolicExpression<Number>::~SymbolicExpression() = default;

// explicit instantiations
template class Core::Utils::SymbolicExpression<double>;

FOUR_C_NAMESPACE_CLOSE
