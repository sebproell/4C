// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_POSITION_HPP
#define FOUR_C_CUT_POSITION_HPP

#include "4C_config.hpp"

#include "4C_cut_element.hpp"
#include "4C_cut_kernel.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  class BoundingBox;
  class PositionFactory;

  /*----------------------------------------------------------------------------*/
  /** \brief Base class of all \e Position objects
   *
   *  */
  class Position
  {
   public:
    /** \brief Status of the position computation.
     *
     *  Please note that the ordering is important, since the assigned integers are
     *  used for comparison operations. As soon as all position cases are working
     *  flawlessly this enum list can and should be adapted.
     *
     *  */
    enum Status
    {
      position_outside_of_bbox =
          -3,  ///< The point is outside of the bounding-box. No computation has been performed.
      position_invalid = -2,     ///< The calculated position is invalid.
      position_zero_area = -1,   ///< A zero element was detected during the position calculation.
      position_unevaluated = 0,  ///< The Compute() routine has not been called, yet.
      position_distance_valid = 1,  ///< The distance but not the position is valid.
      position_valid = 2            ///< The distance and the position are valid.
    };

    //! Map Status enum to std::string
    static inline std::string status_to_string(enum Status pstatus)
    {
      switch (pstatus)
      {
        case position_outside_of_bbox:
          return "position_outside_of_bbox";
        case position_unevaluated:
          return "position_unevaluated";
        case position_zero_area:
          return "position_zero_area";
        case position_invalid:
          return "position_invalid";
        case position_valid:
          return "position_valid";
        case position_distance_valid:
          return "position_distance_valid";
        default:
          return "Unknown PositionStatus";
      }
      exit(EXIT_FAILURE);
    };

   public:
    /// @ name create methods
    /// @{

    /** \brief build variant #1
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param point   (in) : Given global point object.
     *  \param floattype (in) : Floattype to compute geometric operations. */
    static std::shared_ptr<Position> create(
        const Element& element, const Point& point, CutFloatType floattype = floattype_double);

    /** \brief build variant #2
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param floattype (in) : Floattype to compute geometric operations.
     *
     *  */
    template <unsigned rdim>
    static std::shared_ptr<Position> create(const Element& element,
        const Core::LinAlg::Matrix<rdim, 1>& xyz, CutFloatType floattype = floattype_double);
    /** \brief build variant #3
     *
     *  \param xyze    (in) : Global nodal positions of the element.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type.
     *  \param floattype (in) : Floattype to compute geometric operations.
     *
     *  */
    template <unsigned rdim, unsigned cdim, unsigned rdim_2>
    static std::shared_ptr<Position> create(const Core::LinAlg::Matrix<rdim, cdim>& xyze,
        const Core::LinAlg::Matrix<rdim_2, 1>& xyz, const Core::FE::CellType& distype,
        CutFloatType floattype = floattype_double);

    /// \brief build variant #3-1 (Core::LinAlg::SerialDenseMatrix)
    template <unsigned rdim>
    static std::shared_ptr<Position> create(const Core::LinAlg::SerialDenseMatrix& xyze,
        const Core::LinAlg::Matrix<rdim, 1>& xyz, const Core::FE::CellType& distype,
        CutFloatType floattype = floattype_double);

    /** \brief build variant #4
     *
     *  \param nodes   (in) : Nodes of the element we want to check.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type (optional).
     *  \param floattype (in) : Floattype to compute geometric operations.
     *
     *  */
    template <unsigned rdim>
    static std::shared_ptr<Position> create(const std::vector<Node*> nodes,
        const Core::LinAlg::Matrix<rdim, 1>& xyz,
        Core::FE::CellType distype = Core::FE::CellType::dis_none,
        CutFloatType floattype = floattype_double);

    /// @}
   public:
    /// constructor
    Position() { /* intentionally left blank */ };

    /// destructor
    virtual ~Position() = default;

    virtual unsigned n_dim() const = 0;

    virtual unsigned n_prob_dim() const = 0;

    virtual enum Status status() const = 0;

    /** \brief Default Compute method
     *
     *  \param tol (in)        : User defined tolerance for the within_limits check.
     *  \param allow_dist (in) : If TRUE, the method allows an offset in normal direction,
     *                           otherwise the higher-dimensional point \c px_ has to lie
     *                           on the lower dimensional element. */
    bool compute() { return compute(POSITIONTOL); }
    virtual bool compute(const double& tol) = 0;
    bool compute(const bool& allow_dist) { return compute(POSITIONTOL, allow_dist); }
    virtual bool compute(const double& tol, const bool& allow_dist) = 0;

    virtual bool is_given_point_within_element() = 0;

    virtual const double& newton_tolerance() const = 0;

    /*! \brief Return the local coordinates of the given point \c px_ */
    template <class T>
    void local_coordinates(T& rst)
    {
      if (rst.m() < n_dim())
        FOUR_C_THROW(
            "rst has the wrong row number! ( DIM = {} ( rst.m() = {} < DIM ) )", n_dim(), rst.m());
      local_coordinates(rst.data());

      std::fill(rst.data() + n_dim(), rst.data() + rst.m(), 0.0);
    }

    /*! \brief Return the scalar signed perpendicular distance between
     *         given point and embedded element */
    virtual double distance() const
    {
      FOUR_C_THROW("Unsupported for the standard case!");
      exit(EXIT_FAILURE);
    }

    /*! \brief Return the signed perpendicular distance vector between
     *         given point and embedded element.
     *
     *  If you are interested in the distance between a line and a point
     *  in 3-D, you will get 2 distance values by calling this function.
     *  The scalar variant won't work, because it is designed to give you
     *  a signed scalar distance value, what is not possible in this
     *  special case.
     *
     *  \param d (in) : access calculated distance value(s)
     *
     *  */
    template <class T>
    void distance(T& d) const
    {
      if (d.m() < (n_prob_dim() - n_dim()))
        FOUR_C_THROW(
            "The distance vector has the wrong row number! "
            "( DIM = {} ( rst.m() = {} < DIM ) )",
            n_dim(), d.m());

      distance(d.data());
    }

    bool within_limits() const { return within_limits_tol(POSITIONTOL); }
    virtual bool within_limits_tol(const double& Tol) const
    {
      return within_limits_tol(Tol, false);
    }
    bool within_limits(const bool& allow_dist) const
    {
      return within_limits_tol(POSITIONTOL, allow_dist);
    }
    virtual bool within_limits_tol(const double& Tol, const bool& allow_dist) const
    {
      FOUR_C_THROW("Unsupported for the standard case!");
      exit(EXIT_FAILURE);
    }

   protected:
    /*! \brief Return the local coordinates of the given point \c px_ */
    virtual void local_coordinates(double* rst) = 0;

    virtual void distance(double* distance) const = 0;

    /*! \brief Scaling and shifting of the input nodal positions and the location
     *         of the point \c px_
     *
     *  The corner points of the side and the location of the point
     *  are converted into the corresponding range of the natural coordinate
     *  system of the element. This also gives the proper initial choice for
     *  the Newton scheme to find the local coordinate system. */
    virtual void scale_and_shift() = 0;
  };  // class Position

  /*----------------------------------------------------------------------------*/
  /*! \brief Class to check whether a point lies inside an element, and to
   *  compute local coordinates of a point with respect to the element
   *
   *  This class is a template and capable of different kinds of embedded
   *  elements in a higher dimensional space (manifolds). Considered are 1-D elements
   *  embedded in a 2-D space, as well 2-D elements embedded in a 3-D space. Not
   *  considered is the case of 1-D elements embedded in a 3-D space. Nevertheless,
   *  the extension is possible.
   *
   *  The class is a template on the problem dimension \c probdim,
   *                             the element/side type \c eletype,
   *                             the number of nodes per element/side \c numNodesElement,
   *                             the element dimension \c dim.
   *
   *
   */
  template <unsigned probdim, Core::FE::CellType eletype,
      unsigned num_nodes_element = Core::FE::num_nodes<eletype>,
      unsigned dim = Core::FE::dim<eletype>, CutFloatType floattype = floattype_double>
  class PositionGeneric : public Position
  {
   public:
    /** \brief constructor
     *
     *  \param xyze (in) : Global nodal positions of the element.
     *  \param xyz  (in) : Global coordinates of the given point. */
    PositionGeneric(const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
        const Core::LinAlg::Matrix<probdim, 1>& xyz)
        : px_(xyz),
          xyze_(xyze),
          scale_(1.0),
          pos_status_(position_unevaluated),
          compute_tolerance_(-1.0),
          bbside_(nullptr)
    {
      scale_and_shift();
      construct_bounding_box();
    }

    /* Don't call any of these methods directly! Use the base class public
     * methods, instead. */
   protected:
    /// return the element dimension
    unsigned n_dim() const override { return dim; }

    /// return the problem dimension
    unsigned n_prob_dim() const override { return probdim; }

    enum Status status() const override { return pos_status_; }

    /*! \brief Return the local coordinates of the given point \c px_
     *
     *  The xsi_is_valid flag indicates, whether the Compute() function has
     *  been called successfully or not. In the embedded case, it is possible, that
     *  the fall-back method of the Compute function succeeded. Anyway, the
     *  local coordinates are not updated consistently in this case. If you
     *  run into this, you will have to take a closer look.
     *
     *  \param rst: local coordinates.
     *
     *  */
    void local_coordinates(double* rst) override
    {
      if (pos_status_ != position_valid)
      {
        std::ostringstream msg;
        msg << "The local coordinates are not valid. "
               "( Position::Status = "
            << status_to_string(pos_status_) << " )";
        FOUR_C_THROW("{}", msg.str());
      }

      std::copy(xsi_.data(), xsi_.data() + xsi_.m(), rst);
    }

    void distance(double* distance) const override = 0;

    bool compute(const double& tol) override = 0;

    bool is_given_point_within_element() override = 0;

    const double& newton_tolerance() const override
    {
      if (compute_tolerance_ < 0.0) FOUR_C_THROW("Call the Compute() routine first!");
      if (compute_tolerance_ > 1.0)
        FOUR_C_THROW(
            "The compute_tolerance seems not trustworthy, "
            "properly the Compute() routine failed in some way. "
            "( tol = {} )",
            compute_tolerance_);
      return compute_tolerance_;
    }

    /*! \brief Scaling and shifting of the input nodal positions and the location
     *         of the point \c px_
     *
     *  The corner points of the side and the location of the point
     *  are converted into the corresponding range of the natural coordinate
     *  system of the element. This leads to a  proper initial choice for
     *  the Newton scheme to find the local coordinate system. */
    void scale_and_shift() override
    {
      get_element_scale<probdim>(xyze_, scale_);

      px_.scale(1.0 / scale_);
      xyze_.scale(1.0 / scale_);

      get_element_shift<probdim>(xyze_, shift_);

      for (unsigned i = 0; i < num_nodes_element; ++i)
      {
        Core::LinAlg::Matrix<probdim, 1> x1(&xyze_(0, i), true);
        x1.update(-1, shift_, 1);
      }
      px_.update(-1, shift_, 1);
    }

    /** \brief Construct bounding box over the given element
     *  (after scaling) */
    void construct_bounding_box();

   protected:
    /// @name class internal variables
    /// @{

    /// given point location vector (scaled and shifted after constructor call)
    Core::LinAlg::Matrix<probdim, 1> px_;

    /// spatial nodal positions (scaled and shifted after constructor call)
    Core::LinAlg::Matrix<probdim, num_nodes_element> xyze_;

    /** initial scaling of the spatial nodal coordinates \c xyze_
     *  and the location vector \c px_ */
    double scale_;

    /** initial shifting of the spatial nodal coordinates \c xyze_
     *  and the location vector \c px_ */
    Core::LinAlg::Matrix<probdim, 1> shift_;

    /// local coordinates corresponding to \c px_
    Core::LinAlg::Matrix<dim, 1> xsi_;

    /// computation status of the position calculation
    enum Status pos_status_;

    /** contains the tolerance used for the internal Newton method
     *  ( see the Compute() routines ) */
    double compute_tolerance_;

    //! Bounding box over the given embedded element (after scaling is performed)
    std::shared_ptr<BoundingBox> bbside_;

    /// @}

  };  // class Position

  /*----------------------------------------------------------------------------*/
  /** \brief class for the position computation in the standard/non-embedded case
   *
   *  */
  template <unsigned probdim, Core::FE::CellType eletype,
      unsigned num_nodes_element = Core::FE::num_nodes<eletype>,
      unsigned dim = Core::FE::dim<eletype>, CutFloatType floattype = floattype_double>
  class ComputePosition
      : public PositionGeneric<probdim, eletype, num_nodes_element, dim, floattype>
  {
   public:
    /// constructor
    ComputePosition(const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
        const Core::LinAlg::Matrix<probdim, 1>& xyz)
        : PositionGeneric<probdim, eletype, num_nodes_element, dim, floattype>(xyze, xyz) {};

    /* Don't call any of these methods directly! Use the base class public
     * methods, instead. */
   protected:
    void distance(double* distance) const override
    {
      FOUR_C_THROW("Unsupported for the standard case!");
    }

    /** \brief Compute method for the standard case with user defined tolerance
     *
     *  \param Tol (in)        : User defined tolerance for the within_limits check.
     *  \param allow_dist (in) : If TRUE, the method allows an offset in normal direction,
     *                           otherwise the higher-dimensional point \c px_ has to lie
     *                           on the lower dimensional element. */
    bool compute(const double& Tol, const bool& allow_dist) override
    {
      if (allow_dist) FOUR_C_THROW("Compute: allow_dist for ComputePosition not possible!");
      return compute(Tol);
    }
    bool compute(const double& Tol) override;

    bool is_given_point_within_element() override { return Position::compute(); }

    bool within_limits_tol(const double& Tol) const override
    {
      return Kernel::within_limits<eletype>(this->xsi_, Tol);
    }
  };  // class ComputePosition

  /*----------------------------------------------------------------------------*/
  /** \brief class for the position computation in the embedded case
   *
   *  */
  template <unsigned probdim, Core::FE::CellType eletype,
      unsigned num_nodes_element = Core::FE::num_nodes<eletype>,
      unsigned dim = Core::FE::dim<eletype>, CutFloatType floattype = floattype_double>
  class ComputeEmbeddedPosition
      : public PositionGeneric<probdim, eletype, num_nodes_element, dim, floattype>
  {
   public:
    /// constructor
    ComputeEmbeddedPosition(const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
        const Core::LinAlg::Matrix<probdim, 1>& xyz)
        : PositionGeneric<probdim, eletype, num_nodes_element, dim, floattype>(xyze, xyz)
    {
      this->xsi_.set_view(xsi_aug_.data());
    };

    /* Don't call any of these methods directly! Use the base class public
     * methods, instead. */
   protected:
    /** \brief Default Compute method for the embedded case
     *
     *  \param Tol (in)        : User defined tolerance for the within_limits check.
     *  \param allow_dist (in) : If TRUE, the method allows an offset in normal direction,
     *                           otherwise the higher-dimensional point \c px_ has to lie
     *                           on the lower dimensional element. */
    bool compute(const double& tol) override { return compute(tol, false); }
    bool compute(const double& Tol, const bool& allow_dist) override;

    bool is_given_point_within_element() override;

    /*! \brief Return the perpendicular distance between given point to the side */
    double distance() const override
    {
      if (this->pos_status_ < Position::position_distance_valid)
      {
        std::ostringstream msg;
        msg << "Neither the position nor the distance value is valid! "
               "( Position::Status = "
            << Position::status_to_string(this->pos_status_) << " )";
        FOUR_C_THROW("{}", msg.str());
      }

      switch (probdim - dim)
      {
        case 1:
          return xsi_aug_(dim, 0);
        default:
          FOUR_C_THROW("A scalar signed distance value is not available!");
          exit(EXIT_FAILURE);
      }
    }

    void distance(double* distance) const override
    {
      if (this->pos_status_ < Position::position_distance_valid)
      {
        std::ostringstream msg;
        msg << "Neither the position nor the distance value is valid! "
               "( Position::Status = "
            << Position::status_to_string(this->pos_status_) << " )";
        FOUR_C_THROW("{}", msg.str());
      }

      std::copy(xsi_aug_.data() + dim, xsi_aug_.data() + probdim, distance);
    }

    bool within_limits_tol(const double& Tol, const bool& allow_dist) const override
    {
      double tol2 = Tol;
      double tol_xyze = this->xyze_.norm_inf();
      double tol_px = this->px_.norm_inf();
      // choose the weaker tolerance
      if (tol_xyze > tol_px)
        tol2 *= tol_xyze;
      else
        tol2 *= tol_px;

      return Kernel::within_limits_embedded_manifold<probdim, eletype>(
          xsi_aug_, Tol, allow_dist, tol2);
    }

   private:
    /** local coordinates corresponding to \c px_, the last n = ( probdim - dim )
     *  coordinates are the distance values (embedded case only) */
    Core::LinAlg::Matrix<probdim, 1> xsi_aug_;

  };  // class ComputeEmbeddedPosition

  /*----------------------------------------------------------------------------*/
  /** \class Position factory
   *
   *  */
  class PositionFactory
  {
   public:
    /// constructor
    PositionFactory();

    /// get the current problem dimension
    unsigned n_prob_dim() const { return probdim_; }

    /** \brief build variant #1
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param point   (in) : Given global point
     *  \param floattype (in) : Floattype to compute geometric operations. */
    std::shared_ptr<Position> create_position(const Element& element, const Point& point,
        CutFloatType floattype = floattype_double) const;

    /** \brief build variant #2
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param xyz     (in) : Global coordinates of the given point
     *  *  \param floattype (in) : Floattype to compute geometric operations. */
    std::shared_ptr<Position> create_position(
        const Element& element, const double* xyz, CutFloatType floattype = floattype_double) const;

    /** \brief build variant #3
     *
     *  \param xyze    (in) : Global nodal positions of the element.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type.
     *  \param floattype (in) : Floattype to compute geometric operations. */
    std::shared_ptr<Position> create_position(const double* xyze, const double* xyz,
        const Core::FE::CellType& distype, CutFloatType floattype = floattype_double) const;

    /** \brief build variant #4
     *
     *  \param nodes   (in) : Nodes of the element we want to check.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type. (optional)
     *  \param floattype (in) : Floattype to compute geometric operations. */
    std::shared_ptr<Position> create_position(const std::vector<Node*> nodes, const double* xyz,
        Core::FE::CellType distype = Core::FE::CellType::dis_none,
        CutFloatType floattype = floattype_double) const;

    /// \brief specify general floattype for all geometric operations in Cut::POSITION
    static void specify_general_pos_floattype(CutFloatType floattype)
    {
      general_pos_floattype_ = floattype;
    }
    static void specify_general_dist_floattype(CutFloatType floattype)
    {
      general_dist_floattype_ = floattype;
    }

   private:
    /// template class for the actual position creation
    template <bool isembedded, unsigned probdim, unsigned dim, Core::FE::CellType eletype,
        unsigned num_nodes_element = Core::FE::num_nodes<eletype>>
    class PositionCreator
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<probdim, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        return nullptr;
      }
    };

    /// @name Catch a wrong template combination (seems to be necessary for some compilers)
    ///       If more wrong combinations come up, add them here.            hiermeier 04/17
    /// @{

    template <unsigned dim, Core::FE::CellType eletype, unsigned num_nodes_element>
    class PositionCreator<true, dim, dim, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<dim, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<dim, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        FOUR_C_THROW(
            "Wrong template combination: ProbDim must be unequal"
            " element Dim for the embedded case.");
        exit(EXIT_FAILURE);
      }
    };

    template <Core::FE::CellType eletype, unsigned num_nodes_element>
    class PositionCreator<true, 1, 2, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<1, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<1, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        FOUR_C_THROW(
            "Wrong template combination: ProbDim must be larger"
            " than the element Dim for the embedded case.");
        exit(EXIT_FAILURE);
      }
    };

    template <Core::FE::CellType eletype, unsigned num_nodes_element>
    class PositionCreator<true, 1, 3, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<1, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<1, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        FOUR_C_THROW(
            "Wrong template combination: ProbDim must be larger"
            " than the element Dim for the embedded case.");
        exit(EXIT_FAILURE);
      }
    };

    template <Core::FE::CellType eletype, unsigned num_nodes_element>
    class PositionCreator<true, 2, 3, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<2, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<2, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        FOUR_C_THROW(
            "Wrong template combination: ProbDim must be larger"
            " than the element Dim for the embedded case.");
        exit(EXIT_FAILURE);
      }
    };

    template <unsigned probdim, unsigned dim, Core::FE::CellType eletype,
        unsigned num_nodes_element>
    class PositionCreator<false, probdim, dim, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<probdim, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        FOUR_C_THROW(
            "Wrong template combination: ProbDim must be equal"
            " element Dim for the standard case.");
        exit(EXIT_FAILURE);
      }
    };

    /// @}

    /// create an embedded position object ( dim < probdim )
    template <unsigned probdim, unsigned dim, Core::FE::CellType eletype,
        unsigned num_nodes_element>
    class PositionCreator<true, probdim, dim, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<probdim, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        switch (use_dist_floattype(floattype))
        {
          case floattype_double:
          {
            return new ComputeEmbeddedPosition<probdim, eletype, num_nodes_element,
                Core::FE::dim<eletype>, floattype_double>(xyze, xyz);
            break;
          }
          case floattype_cln:
          {
            return new ComputeEmbeddedPosition<probdim, eletype, num_nodes_element,
                Core::FE::dim<eletype>, floattype_cln>(xyze, xyz);
            break;
          }
          default:
            FOUR_C_THROW("Unsupported floattype!");
        }
        return nullptr;
      }
    };

    /// create a standard position object ( probdim == dim )
    template <unsigned dim, Core::FE::CellType eletype, unsigned num_nodes_element>
    class PositionCreator<false, dim, dim, eletype, num_nodes_element>
    {
     public:
      static Position* create(const Core::LinAlg::Matrix<dim, num_nodes_element>& xyze,
          const Core::LinAlg::Matrix<dim, 1>& xyz, CutFloatType floattype = floattype_double)
      {
        switch (use_pos_floattype(floattype))
        {
          case floattype_double:
          {
            return new ComputePosition<dim, eletype, num_nodes_element, Core::FE::dim<eletype>,
                floattype_double>(xyze, xyz);
            break;
          }
          case floattype_cln:
          {
            return new ComputePosition<dim, eletype, num_nodes_element, Core::FE::dim<eletype>,
                floattype_cln>(xyze, xyz);
            break;
          }
          default:
            FOUR_C_THROW("Unsupported floattype!");
        }
        return nullptr;
      }
    };

   private:
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #1                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{

    /** build variant #1
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param point   (in) : Given global point
     *
     *  */
    template <unsigned probdim, Core::FE::CellType eletype, unsigned dim = Core::FE::dim<eletype>,
        unsigned num_nodes_element = Core::FE::num_nodes<eletype>>
    std::shared_ptr<Position> build_position(
        const Element& element, const Point& point, CutFloatType floattype = floattype_double) const
    {
      Core::LinAlg::Matrix<probdim, num_nodes_element> xyze;
      element.coordinates(xyze);
      Core::LinAlg::Matrix<probdim, 1> px;
      point.coordinates(px.data());

      if (probdim > dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<true, probdim, dim, eletype>::create(xyze, px, floattype));
      }
      else if (probdim == dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<false, probdim, dim, eletype>::create(xyze, px, floattype));
      }
      else
        FOUR_C_THROW(
            "The element dimension is larger than the problem dimension! \n"
            "dim = {}, probdim = {}",
            dim, probdim);

      exit(EXIT_FAILURE);
    }

    /// concrete create variant #1
    template <Core::FE::CellType eletype>
    std::shared_ptr<Position> create_concrete_position(
        const Element& element, const Point& point, CutFloatType floattype = floattype_double) const
    {
      const int dim = Core::FE::dim<eletype>;
      if (dim > probdim_)
        FOUR_C_THROW(
            "The element dimension is larger than the problem dimension! \n"
            "dim = {}, probdim = {}",
            dim, probdim_);

      switch (probdim_)
      {
        case 2:
          return build_position<2, eletype>(element, point, floattype);
        case 3:
          return build_position<3, eletype>(element, point, floattype);
        default:
          FOUR_C_THROW("Unsupported problem dimension! (probdim = {})", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #2                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{
   public:
    /** \brief build variant #2
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param xyz     (in) : Global coordinates of the given point.
     *
     *  */
    template <unsigned probdim, Core::FE::CellType eletype, unsigned dim = Core::FE::dim<eletype>,
        unsigned num_nodes_element = Core::FE::num_nodes<eletype>>
    static std::shared_ptr<Position> build_position(const Element& element,
        const Core::LinAlg::Matrix<probdim, 1>& xyz, CutFloatType floattype = floattype_double)
    {
      Core::LinAlg::Matrix<probdim, num_nodes_element> xyze;
      element.coordinates(xyze);

      if (probdim > dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<true, probdim, dim, eletype>::create(xyze, xyz, floattype));
      }
      else if (probdim == dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<false, probdim, dim, eletype>::create(xyze, xyz, floattype));
      }
      else
        FOUR_C_THROW(
            "The element dimension is larger than the problem dimension! \n"
            "dim = {}, probdim = {}",
            dim, probdim);

      exit(EXIT_FAILURE);
    }

   private:
    /// concrete create variant #2
    template <Core::FE::CellType eletype>
    std::shared_ptr<Position> create_concrete_position(
        const Element& element, const double* xyz, CutFloatType floattype = floattype_double) const
    {
      const int dim = Core::FE::dim<eletype>;
      if (dim > probdim_)
        FOUR_C_THROW(
            "The element dimension is larger than the problem dimension! \n"
            "dim = {}, probdim = {}",
            dim, probdim_);

      switch (probdim_)
      {
        case 2:
        {
          Core::LinAlg::Matrix<2, 1> xyz_mat(xyz, true);
          return build_position<2, eletype>(element, xyz_mat, floattype);
        }
        case 3:
        {
          Core::LinAlg::Matrix<3, 1> xyz_mat(xyz, true);
          return build_position<3, eletype>(element, xyz_mat, floattype);
        }
        default:
          FOUR_C_THROW("Unsupported problem dimension! (probdim = {})", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #3                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{
   public:
    /** \brief build variant #3
     *
     *  \param xyze    (in) : Global nodal positions of the element.
     *  \param xyz     (in) : Global coordinates of the given point.
     *
     *  */
    template <unsigned probdim, Core::FE::CellType eletype, unsigned dim = Core::FE::dim<eletype>,
        unsigned num_nodes_element = Core::FE::num_nodes<eletype>>
    static std::shared_ptr<Position> build_position(
        const Core::LinAlg::Matrix<probdim, num_nodes_element>& xyze,
        const Core::LinAlg::Matrix<probdim, 1>& xyz, CutFloatType floattype = floattype_double)
    {
      if (probdim > dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<true, probdim, dim, eletype>::create(xyze, xyz, floattype));
      }
      else if (probdim == dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<false, probdim, dim, eletype>::create(xyze, xyz, floattype));
      }
      else
        FOUR_C_THROW(
            "The element dimension is larger than the problem dimension! \n"
            "dim = {}, probdim = {}",
            dim, probdim);

      exit(EXIT_FAILURE);
    }

   private:
    /// concrete create variant #3
    template <Core::FE::CellType eletype>
    std::shared_ptr<Position> create_concrete_position(
        const double* xyze, const double* xyz, CutFloatType floattype = floattype_double) const
    {
      const unsigned num_nodes_ele = Core::FE::num_nodes<eletype>;
      switch (probdim_)
      {
        case 2:
        {
          Core::LinAlg::Matrix<2, num_nodes_ele> xyze_mat(xyze, true);
          Core::LinAlg::Matrix<2, 1> xyz_mat(xyz, true);
          return build_position<2, eletype>(xyze_mat, xyz_mat, floattype);
        }
        case 3:
        {
          Core::LinAlg::Matrix<3, num_nodes_ele> xyze_mat(xyze, true);
          Core::LinAlg::Matrix<3, 1> xyz_mat(xyz, true);
          return build_position<3, eletype>(xyze_mat, xyz_mat, floattype);
        }
        default:
          FOUR_C_THROW("Unsupported problem dimension! (probdim = {})", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #4                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{
   public:
    /** \brief build variant #4
     *
     *  \param nodes (in) : Nodes of the element we want to check.
     *  \param xyz   (in) : Global coordinates of the given point.
     *
     *  */
    template <unsigned probdim, Core::FE::CellType eletype, unsigned dim = Core::FE::dim<eletype>,
        unsigned num_nodes_element = Core::FE::num_nodes<eletype>>
    static std::shared_ptr<Position> build_position(const std::vector<Node*> nodes,
        const Core::LinAlg::Matrix<probdim, 1>& xyz, CutFloatType floattype = floattype_double)
    {
      if (not(nodes.size() == num_nodes_element))
        FOUR_C_THROW(
            "node number for this element is not correct\n"
            "numNodesElement = {}, nodes.size() = {}",
            num_nodes_element, nodes.size());

      Core::LinAlg::Matrix<probdim, num_nodes_element> xyze;
      for (unsigned i = 0; i < num_nodes_element; ++i)
      {
        Node* n = nodes[i];
        n->coordinates(&xyze(0, i));
      }

      if (probdim > dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<true, probdim, dim, eletype>::create(xyze, xyz, floattype));
      }
      else if (probdim == dim)
      {
        return std::shared_ptr<Position>(
            PositionCreator<false, probdim, dim, eletype>::create(xyze, xyz, floattype));
      }
      else
        FOUR_C_THROW(
            "The element dimension is larger than the problem dimension! \n"
            "dim = {}, probdim = {}",
            dim, probdim);

      exit(EXIT_FAILURE);
    }

   private:
    /// concrete create variant #4
    template <Core::FE::CellType eletype>
    std::shared_ptr<Position> create_concrete_position(const std::vector<Node*> nodes,
        const double* xyz, CutFloatType floattype = floattype_double) const
    {
      switch (probdim_)
      {
        case 2:
        {
          Core::LinAlg::Matrix<2, 1> xyz_mat(xyz, true);
          return build_position<2, eletype>(nodes, xyz_mat, floattype);
        }
        case 3:
        {
          Core::LinAlg::Matrix<3, 1> xyz_mat(xyz, true);
          return build_position<3, eletype>(nodes, xyz_mat, floattype);
        }
        default:
          FOUR_C_THROW("Unsupported problem dimension! (probdim = {})", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}

    /*! \brief get general floattype for all geometric operations in Cut::POSITION*/
    static CutFloatType use_pos_floattype(CutFloatType floattype);
    static CutFloatType use_dist_floattype(CutFloatType floattype);

    /*! \brief general floattype for all geometric operations in Cut::POSITION*/
    static CutFloatType general_pos_floattype_;
    static CutFloatType general_dist_floattype_;

    /// problem dimension
    unsigned probdim_;
  };  // class PositionFactory


}  // namespace Cut



FOUR_C_NAMESPACE_CLOSE

#endif
