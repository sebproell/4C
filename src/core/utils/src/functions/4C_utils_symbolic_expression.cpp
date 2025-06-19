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


namespace Core::Utils::SymbolicExpressionDetails
{
  using IndexType = std::size_t;
  constexpr IndexType invalid_index = -1uL;

  enum class NodeType : unsigned char
  {
    // numeric value
    number,
    // unary function e.g. 'sin', 'cos', 'exp', ...
    unary_function,
    // binary function e.g. '+', '-', '*', '/', ...
    binary_function,
    // independent variable: 't', 'x', ...
    variable,
  };

  enum class UnaryFunctionType : unsigned char
  {
    acos,
    asin,
    atan,
    cos,
    sin,
    tan,
    cosh,
    sinh,
    tanh,
    exp,
    log,
    log10,
    sqrt,
    heaviside,
    fabs,
  };

  enum class BinaryFunctionType : unsigned char
  {
    plus,
    minus,
    mul,
    div,
    atan2,
    pow,
  };

  /**
   * \brief Syntax tree node holding binary, unary operator, variables or literals.
   */
  struct SyntaxTreeNode
  {
    /**
     * The different options that could be stored in the node.
     */
    union Value
    {
      double number{};
      UnaryFunctionType unary_function;
      BinaryFunctionType binary_function;
      IndexType var_index;  // refer to a variable stored in the parser
    } as;

    //! Lhs and rhs nodes
    IndexType lhs_index{-1u};
    IndexType rhs_index{-1u};

    //! type of the node, ie operator, literal, ...
    NodeType type;
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
      variable_constants_.clear();

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

   private:
    void parse();
    IndexType parse_primary(Lexer& lexer);
    IndexType parse_pow(Lexer& lexer);
    IndexType parse_term(Lexer& lexer);
    IndexType parse_expr(Lexer& lexer);
    IndexType parse(Lexer& lexer);

    //! Create a new node in the node arena and return a pointer to it.
    //! Optionally, a left-hand side and right-hand side node index can be provided.
    IndexType create_node(NodeType type, SyntaxTreeNode::Value value, IndexType lhs = invalid_index,
        IndexType rhs = invalid_index);

    //! Either return the index of @p var or give it a new one.
    IndexType create_variable(std::string_view var);

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

    SyntaxTreeNode& at(IndexType index)
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
    IndexType root_{invalid_index};

    //! set of all parsed variables
    std::vector<std::string_view> parsed_variable_constant_names_;

    //! During evaluation, this map hold the values of the variables and constants at
    //! a specified index.
    mutable std::vector<T> variable_constants_;
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
    //! check if function has been parsed
    FOUR_C_ASSERT(!node_arena_.empty(), "Internal error");

    variable_constants_.resize(parsed_variable_constant_names_.size());
    for (unsigned int i = 0; i < parsed_variable_constant_names_.size(); i++)
    {
      if (auto it = variable_values.find(std::string(parsed_variable_constant_names_[i]));
          it != variable_values.end())
      {
        variable_constants_[i] = it->second;
      }
      else if (auto it = constants.find(std::string(parsed_variable_constant_names_[i]));
          it != constants.end())
      {
        variable_constants_[i] = it->second;
      }
      else
      {
        std::string evaluate_variable_names = std::accumulate(variable_values.begin(),
            variable_values.end(), std::string(), [](const std::string& acc, const auto& v)
            { return acc.empty() ? v.first : acc + ", " + v.first; });

        std::string evaluate_constant_names = std::accumulate(constants.begin(), constants.end(),
            std::string(), [](const std::string& acc, const auto& v)
            { return acc.empty() ? v.first : acc + ", " + v.first; });

        FOUR_C_THROW(
            "Some variables that this parser encountered in the expression are not passed to "
            "the evaluate function.\n\n"
            "Expression:  {} \n"
            "Variables passed: {} \n"
            "Constants passed: {}",
            symbolicexpression_.c_str(), evaluate_variable_names.c_str(),
            evaluate_constant_names.c_str());
      }
    }


    //! evaluate syntax tree of function depending on set variables
    auto result = this->interpret(at(root_));

    return result;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse primary entities, i.e. literals and unary operators,
         such as numbers, parentheses, independent variables, operator names
  */
  template <class T>
  IndexType Parser<T>::parse_primary(Lexer& lexer)
  {
    IndexType lhs = invalid_index;
    switch (lexer.tok_)
    {
      case Lexer::tok_lpar:
        lexer.lexan();
        lhs = parse_expr(lexer);
        if (lexer.tok_ != Lexer::tok_rpar) FOUR_C_THROW("')' expected");
        lexer.lexan();
        break;
      case Lexer::tok_int:
        lhs = create_node(NodeType::number, {.number = static_cast<double>(lexer.integer_)});
        lexer.lexan();
        break;
      case Lexer::tok_real:
        lhs = create_node(NodeType::number, {.number = lexer.real_});
        lexer.lexan();
        break;
      case Lexer::tok_sub:
      {
        // This is a unary minus operator.
        lexer.lexan();
        IndexType rhs = parse_pow(lexer);
        auto& rhs_node = at(rhs);
        if (rhs_node.type == NodeType::number)
        {
          rhs_node.as.number *= -1;
          lhs = rhs;
        }
        else
        {
          lhs = create_node(NodeType::number, {.number = -1});
          lhs = create_node(
              NodeType::binary_function, {.binary_function = BinaryFunctionType::mul}, lhs, rhs);
        }
        break;
      }
      case Lexer::tok_name:
      {
        // get substring starting from str_ with length of lexer.integer_
        std::string_view name(lexer.str_, lexer.integer_);
        if (name == "pi")
        {
          lhs = create_node(NodeType::number, {.number = M_PI});
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
            auto function_name = *EnumTools::enum_cast<UnaryFunctionType>(name);
            lhs = create_node(
                NodeType::unary_function, {.unary_function = function_name}, parse_expr(lexer));

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

            lhs = create_node(NodeType::binary_function,
                {.binary_function = BinaryFunctionType::atan2}, first_arg, second_arg);

            // Consume closing parenthesis
            if (lexer.tok_ != Lexer::tok_rpar)
              FOUR_C_THROW("')' expected for function name '{}'", name);
            lexer.lexan();
            break;
          }
          else
          {
            lhs = create_node(NodeType::variable, {.var_index = create_variable(name)});
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
  IndexType Parser<T>::parse_pow(Lexer& lexer)
  {
    IndexType lhs;
    IndexType rhs;

    lhs = parse_primary(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_pow)
      {
        lexer.lexan();
        rhs = parse_primary(lexer);
        lhs = create_node(
            NodeType::binary_function, {.binary_function = BinaryFunctionType::pow}, lhs, rhs);
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
  IndexType Parser<T>::parse_term(Lexer& lexer)
  {
    IndexType lhs;
    IndexType rhs;

    lhs = parse_pow(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_mul)
      {
        lexer.lexan();
        rhs = parse_pow(lexer);
        lhs = create_node(
            NodeType::binary_function, {.binary_function = BinaryFunctionType::mul}, lhs, rhs);
      }
      else if (lexer.tok_ == Lexer::tok_div)
      {
        lexer.lexan();
        rhs = parse_pow(lexer);
        lhs = create_node(
            NodeType::binary_function, {.binary_function = BinaryFunctionType::div}, lhs, rhs);
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
  IndexType Parser<T>::parse(Lexer& lexer)
  {
    IndexType lhs = parse_expr(lexer);

    // check if parsing ended before processing the entire string
    if (lexer.tok_ != Lexer::tok_done)
    {
      FOUR_C_THROW("Invalid syntax: The remaining string '{}' is not parsed.",
          lexer.funct_.data() + lexer.pos_);
    }

    return lhs;
  }

  template <class T>
  IndexType Parser<T>::create_node(
      NodeType type, SyntaxTreeNode::Value value, IndexType lhs, IndexType rhs)
  {
    SyntaxTreeNode& node = node_arena_.emplace_back();
    node.as = value;
    node.type = type;
    node.lhs_index = lhs;
    node.rhs_index = rhs;
    // Index of the node that was just created.
    return node_arena_.size() - 1;
  }

  template <class T>
  IndexType Parser<T>::create_variable(std::string_view var)
  {
    auto it = std::ranges::find(parsed_variable_constant_names_, var);
    if (it != parsed_variable_constant_names_.end())
    {
      return std::distance(parsed_variable_constant_names_.begin(), it);
    }
    else
    {
      parsed_variable_constant_names_.emplace_back(var);
      return parsed_variable_constant_names_.size() - 1;
    }
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by addition or subtraction: a+b, a-b
  */
  template <class T>
  IndexType Parser<T>::parse_expr(Lexer& lexer)
  {
    IndexType lhs = parse_term(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_add)
      {
        lexer.lexan();
        IndexType rhs = parse_term(lexer);
        lhs = create_node(
            NodeType::binary_function, {.binary_function = BinaryFunctionType::plus}, lhs, rhs);
      }
      else if (lexer.tok_ == Lexer::tok_sub)
      {
        lexer.lexan();
        IndexType rhs = parse_term(lexer);
        lhs = create_node(
            NodeType::binary_function, {.binary_function = BinaryFunctionType::minus}, lhs, rhs);
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
      case NodeType::number:
      {
        res = node.as.number;
        break;
      }
      case NodeType::binary_function:
      {
        T lhs;
        T rhs;

        // recursively visit branches and obtain sub-results
        lhs = interpret(at(node.lhs_index));
        rhs = interpret(at(node.rhs_index));

        // evaluate the node operator
        switch (node.as.binary_function)
        {
          case BinaryFunctionType::plus:
            res = lhs + rhs;
            break;
          case BinaryFunctionType::minus:
            res = lhs - rhs;
            break;
          case BinaryFunctionType::mul:
            res = lhs * rhs;
            break;
          case BinaryFunctionType::div:
            res = lhs / rhs;
            break;
          case BinaryFunctionType::pow:
            res = std::pow(lhs, rhs);
            break;
          case BinaryFunctionType::atan2:
            res = std::atan2(lhs, rhs);
            break;
        }
        break;
      }
      case NodeType::variable:
      {
        res = variable_constants_[node.as.var_index];
        break;
      }
      case NodeType::unary_function:
      {
        T arg = interpret(at(node.lhs_index));
        switch (node.as.unary_function)
        {
          case UnaryFunctionType::acos:
            res = acos(arg);
            break;
          case UnaryFunctionType::asin:
            res = asin(arg);
            break;
          case UnaryFunctionType::atan:
            res = atan(arg);
            break;
          case UnaryFunctionType::cos:
            res = cos(arg);
            break;
          case UnaryFunctionType::sin:
            res = sin(arg);
            break;
          case UnaryFunctionType::tan:
            res = tan(arg);
            break;
          case UnaryFunctionType::cosh:
            res = cosh(arg);
            break;
          case UnaryFunctionType::sinh:
            res = sinh(arg);
            break;
          case UnaryFunctionType::tanh:
            res = tanh(arg);
            break;
          case UnaryFunctionType::exp:
            res = exp(arg);
            break;
          case UnaryFunctionType::log:
            res = log(arg);
            break;
          case UnaryFunctionType::log10:
            res = log10(arg);
            break;
          case UnaryFunctionType::sqrt:
            res = sqrt(arg);
            break;
          case UnaryFunctionType::fabs:
            res = fabs(arg);
            break;
          case UnaryFunctionType::heaviside:
          {
            if (arg > 0)
            {
              res = 1.0;
            }
            else
            {
              res = 0.0;
            }
            break;
          }
        }
        break;
      }
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

template <typename Number>
Core::Utils::SymbolicExpression<Number>::SymbolicExpression(
    SymbolicExpression&& other) noexcept = default;

template <typename Number>
Core::Utils::SymbolicExpression<Number>& Core::Utils::SymbolicExpression<Number>::operator=(
    SymbolicExpression&& other) noexcept = default;

// explicit instantiations
template class Core::Utils::SymbolicExpression<double>;

FOUR_C_NAMESPACE_CLOSE
