/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of Mortar matrices on the overlap
       of two MORTAR::Elements in 1D and 2D (derived version for contact)

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_INTEGRATOR_HPP
#define FOUR_C_CONTACT_INTEGRATOR_HPP

#include "baci_config.hpp"

#include "baci_inpar_wear.hpp"
#include "baci_mortar_integrator.hpp"
#include "baci_utils_pairedvector.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
}
namespace MORTAR
{
  class ParamsInterface;
}

namespace CONTACT
{
  // forward declaration
  class ParamsInterface;
  /*!
   \brief A class to perform Gaussian integration and assembly of Mortar
   matrices on the overlap of two MORTAR::Elements (1 Slave, 1 Master)
   in 1D (which is equivalent to a 2D coupling problem) and in 2D
   (which is equivalent to a 3D coupling problem).

   This is a derived class from MORTAR::Integrator which does
   the contact-specific stuff for 3d mortar coupling.

   */

  class Integrator
  {
   public:
    /*!
     \brief Constructor  with shape function specification

     Constructs an instance of this class using a specific type of shape functions.<br>
     Note that this is \b not a collective call as overlaps are
     integrated in parallel by individual processes.<br>
     Note also that this constructor relies heavily on the
     CORE::FE::IntegrationPoints structs to get Gauss points
     and corresponding weights.

     */
    Integrator(Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm);

    /*!
     \brief Destructor

     */
    virtual ~Integrator() = default;

    //! don't want = operator
    Integrator operator=(const Integrator& old) = delete;
    //! don't want copy constructor
    Integrator(const Integrator& old) = delete;

    //! get specified integration type
    inline enum INPAR::MORTAR::IntType IntegrationType() const { return integrationtype_; }

    const Epetra_Comm& Comm() const { return Comm_; }

    //! @name 2D and 3D integration methods

    /*!
     \brief check for boundary segmentation in 2D

     */
    bool BoundarySegmCheck2D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles);

    /*!
     \brief check for boundary segmentation in 2D

     */
    bool BoundarySegmCheck3D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles);


    /*!
     \brief Build all integrals and linearizations without segmentation -- 2D
     (i.e. M, g, LinM, Ling and possibly D, LinD)

     */
    virtual void IntegrateDerivEle2D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    virtual void IntegrateDerivEle2D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    /*!
     \brief integrate D matrix without lin...

     */
    void IntegrateD(MORTAR::Element& sele, const Epetra_Comm& comm, bool lin = false);

    /*!
     \brief Build all integrals and linearizations on a 1D slave /
     master overlap (i.e. M, g, LinM, Ling and possibly D, LinD and
     wear)

     */
    virtual void IntegrateDerivSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
        MORTAR::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    virtual void IntegrateDerivSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
        MORTAR::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    /*!
     \brief Build all integrals and linearizations without segmentation -- 3D
     (i.e. M, g, LinM, Ling and possibly D, LinD)

     */
    virtual void IntegrateDerivEle3D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    virtual void IntegrateDerivEle3D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    /*!
     \brief Build all integrals and linearizations on a 2D slave /
     master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)
     for the auxiliary plane coupling case

     */
    virtual void IntegrateDerivCell3DAuxPlane(MORTAR::Element& sele, MORTAR::Element& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    virtual void IntegrateDerivCell3DAuxPlane(MORTAR::Element& sele, MORTAR::Element& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    /*!
     \brief Build all integrals and linearizations on a 2D slave /
     master integration cell (i.e. M, g, LinM, Ling) for
     the auxiliary plane coupling case with quadratic interpolation

     */
    void IntegrateDerivCell3DAuxPlaneQuad(MORTAR::Element& sele, MORTAR::Element& mele,
        MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn);

    /*!
     \brief ....

     */
    void IntegrateDerivCell3DAuxPlaneLTS(MORTAR::Element& sele, MORTAR::Element& lsele,
        MORTAR::Element& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
        const Epetra_Comm& comm);

    /*!
     \brief ....

     */
    void IntegrateDerivCell3DAuxPlaneSTL(MORTAR::Element& mele, MORTAR::Element& lele,
        MORTAR::Element& sele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
        const Epetra_Comm& comm);

    /*!
     \brief Compute penalty scaling factor kappa on slave element

     */
    void IntegrateKappaPenalty(MORTAR::Element& sele, double* sxia, double* sxib,
        Teuchos::RCP<CORE::LINALG::SerialDenseVector> gseg);


    /*!
     \brief Compute penalty scaling factor kappa on slave element for LTS algorithm
            for last converged configuration

     */
    void IntegrateKappaPenaltyLTS(MORTAR::Element& sele);

    /*!
     \brief Compute penalty scaling factor kappa on slave integration element
     (special version for the 3D quadratic case)

     */
    void IntegrateKappaPenalty(MORTAR::Element& sele, MORTAR::IntElement& sintele, double* sxia,
        double* sxib, Teuchos::RCP<CORE::LINALG::SerialDenseVector> gseg);

    //@}

    //! @name 2D and 3D linearization methods

    /*!
     \brief Compute directional derivative of segment end coordinates
     Xi on a 1D slave / master overlap

     */
    void DerivXiAB2D(MORTAR::Element& sele, double& sxia, double& sxib, MORTAR::Element& mele,
        double& mxia, double& mxib, std::vector<CORE::GEN::pairedvector<int, double>>& derivxi,
        bool& startslave, bool& endslave, int& linsize);

    /*!
     \brief Compute directional derivative of master Gauss point
     coordinates XiGP on a 1D slave / master overlap

     */
    void DerivXiGP2D(MORTAR::Element& sele, MORTAR::Element& mele, double sxigp, double mxigp,
        const CORE::GEN::pairedvector<int, double>& derivsxi,
        CORE::GEN::pairedvector<int, double>& derivmxi, int& linsize);

    /*!
     \brief Compute directional derivative of master Gauss point
     coordinates XiGP on a 2D slave / master integration cell

     */
    void DerivXiGP3D(MORTAR::Element& sele, MORTAR::Element& mele, const double* sxigp,
        const double* mxigp, const std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi, double alpha);

    /*!
     \brief Compute directional derivative of slave / master Gauss point
     coordinates XiGP on a 2D slave / master integration cell
     (This is the AuxPlane version, thus master and slave are projected)

     */
    void DerivXiGP3DAuxPlane(MORTAR::Element& ele, double* xigp, double* auxn,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivxi, double& alpha,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivauxn,
        CORE::GEN::pairedvector<int, CORE::LINALG::Matrix<3, 1>>& derivgp);

    /*!
     \brief Assemble g~ contribution of current overlap into slave nodes

     */
    bool AssembleG(
        const Epetra_Comm& comm, MORTAR::Element& sele, CORE::LINALG::SerialDenseVector& gseg);

    /*!
     \brief Assemble g~ contribution of current overlap into slave nodes
     (special version for 3D quadratic mortar with piecewise linear LM interpolation)

     */
    bool AssembleG(const Epetra_Comm& comm, MORTAR::IntElement& sintele,
        CORE::LINALG::SerialDenseVector& gseg);

    // GP calls
    /*!
     \brief Return number of Gauss points for this instance

     */
    int& nGP() { return ngp_; }

    /*!
     \brief Return coordinates of a specific GP in 1D/2D CElement

     */
    double& Coordinate(int gp, int dir) { return coords_(gp, dir); }

    /*!
     \brief Return weight of a specific GP in 1D/2D CElement

     */
    double& Weight(int gp) { return weights_[gp]; }

    /*!
     \brief Get problem dimension

     Note that only 2D and 3D are possible here as this refers to the global
     problem dimension. On integration level this corresponds to 1D integration
     (dim_==2) and 2D integration (dim_==3) on the interface!

     */
    int Dim() const { return dim_; };

   protected:
    /*!
     \brief Initialize Gauss rule (points, weights) for this MORTAR::Integrator

     */
    void InitializeGP(CORE::FE::CellType eletype);

    /*!
     * @brief Perform integration at Gauss point for 3D problems.
     * This is where the distinction between methods should be, i.e. mortar, augmented, gpts,...
     *
     * @param[in] sele     current mortar slave element
     * @param[in] mele     current mortar master element
     * @param[in] sval     slave side shape function evaluated at current Gauss point
     * @param[in] lmval    Lagrangian multiplier shape function evaluated at current Gauss point
     * @param[in] mval     master side shape function evaluated at current Gauss point
     * @param[in] sderiv   slave side shape function derivative at current Gauss point
     * @param[in] mderiv   master side shape function derivative at current Gauss point
     * @param[in] lmderiv  Lagrangian multiplier shape function derivative evaluated at current
     *                     Gauss point
     * @param[in] dualmap  directional derivative of dual shape functions
     * @param[in] wgt      Gauss point weight
     * @param[in] jac           Jacobian determinant of integration cell
     * @param[in] derivjac      directional derivative of cell Jacobian
     * @param[in] normal        integration cell normal
     * @param[in] dnmap_unit    directional derivative of integration cell normal
     * @param[in] gap           gap
     * @param[in] deriv_gap     directional derivative of gap
     * @param[in] sxi       slave side Gauss point coordinates
     * @param[in] mxi       master side Gauss point coordinates
     * @param[in] derivsxi  directional derivative of slave side Gauss point coordinates
     * @param[in] derivmxi  directional derivative of master side Gauss point coordinates
     */
    virtual void IntegrateGP_3D(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi);

    /*!
     * @brief Perform integration at Gauss point for 2D problems.
     * This is where the distinction between methods should be, i.e. mortar, augmented, gpts,...
     *
     * @param[in] sele     mortar slave element
     * @param[in] mele     mortar master element
     * @param[in] sval     slave side shape function evaluated at current Gauss point
     * @param[in] lmval    Lagrangian multiplier shape function evaluated at current Gauss point
     * @param[in] mval     master side shape function evaluated at current Gauss point
     * @param[in] sderiv   slave side shape function derivative at current Gauss point
     * @param[in] mderiv   master side shape function derivative at current Gauss point
     * @param[in] lmderiv  Lagrangian multiplier shape function derivative evaluated at current
     *                     Gauss point
     * @param[in] dualmap  directional derivative of dual shape functions
     * @param[in] wgt      Gauss point weight
     * @param[in] jac           Jacobian determinant of integration cell
     * @param[in] derivjac      directional derivative of cell Jacobian
     * @param[in] normal        integration cell normal
     * @param[in] dnmap_unit    directional derivative of integration cell normal
     * @param[in] gap           gap
     * @param[in] deriv_gap     directional derivative of gap
     * @param[in] sxi       slave side Gauss point coordinates
     * @param[in] mxi       master side Gauss point coordinates
     * @param[in] derivsxi  directional derivative of slave side Gauss point coordinates
     * @param[in] derivmxi  directional derivative of master side Gauss point coordinates
     */
    virtual void IntegrateGP_2D(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi);

    /*!
     \brief evaluate D2-matrix entries at GP

     */
    void inline GP_D2(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& lm2val, CORE::LINALG::SerialDenseVector& m2val,
        double& jac, double& wgt, const Epetra_Comm& comm);

    /*!
     \brief evaluate D/M-matrix entries at GP

     */
    void GP_DM(MORTAR::Element& sele, MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval, double& jac,
        double& wgt, bool& bound);

    /*!
     \brief evaluate D/M-matrix entries at GP (3D quadratic)

     */
    void inline GP_3D_DM_Quad(MORTAR::Element& sele, MORTAR::Element& mele,
        MORTAR::IntElement& sintele, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& lmintval, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, const double& jac, double& wgt, const int& nrow,
        const int& nintrow, const int& ncol, const int& ndof, bool& bound);

    /*!
     \brief lin D/M-matrix entries at GP for bound case

     */
    void inline GP_2D_DM_Lin_bound(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double& wgt, const CORE::GEN::pairedvector<int, double>& derivjac,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin D/M-matrix entries at GP for bound case

     */
    void inline GP_2D_DM_Lin(int& iter, bool& bound, bool& linlm, MORTAR::Element& sele,
        MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        const CORE::GEN::pairedvector<int, double>& derivjac,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin D/M-matrix entries at GP for elebased integration

     */
    void inline GP_2D_DM_Ele_Lin(int& iter, bool& bound, MORTAR::Element& sele,
        MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& mderiv, double& dxdsxi, double& wgt,
        const CORE::GEN::pairedvector<int, double>& dmxigp,
        const CORE::GEN::pairedvector<int, double>& derivjac,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin D/M-matrix entries at GP

     */
    void GP_3D_DM_Lin(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& wgt, double& jac, std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin D/M-matrix entries at GP for bound case

     */
    void inline GP_3D_DM_Lin_bound(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
        double& jac, double& wgt, const CORE::GEN::pairedvector<int, double>& derivjac,
        std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin D/M-matrix entries at GP for bound case (3D quad)

     */
    void inline GP_3D_DM_Quad_Lin(bool& duallin, MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& svalmod,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, double& wgt, double& jac,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dpsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dpmxigp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        bool dualquad3d);

    void inline GP_3D_DM_Quad_pwlin_Lin(int& iter, MORTAR::Element& sele, MORTAR::Element& sintele,
        MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseVector& lmintval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseMatrix& lmintderiv, double& wgt, double& jac,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dpsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dpmxigp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap);

    /*!
     \brief evaluate weighted Gap entries at GP

     */
    void GP_3D_wGap(MORTAR::Element& sele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& lmval, double* gap, double& jac, double& wgt,
        bool quadratic, int nintrow = 0);

    /*!
     \brief evaluate weighted Gap entries at GP

     */
    void inline GP_2D_wGap(MORTAR::Element& sele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& lmval, double* gap, double& jac, double& wgt);

    /*!
     \brief evaluate geometrical gap at GP
     */
    void Gap_3D(MORTAR::Element& sele, MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
        std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        CORE::GEN::pairedvector<int, double>& dgapgp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit);


    /*!
     \brief evaluate geometrical gap at GP
     */
    void Gap_2D(MORTAR::Element& sele, MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
        std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        CORE::GEN::pairedvector<int, double>& dgapgp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit);

    void inline GP_2D_G_Lin(int& iter, MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, double& gap, double* gpn, double& jac,
        double& wgt, CORE::GEN::pairedvector<int, double>& dgapgp,
        CORE::GEN::pairedvector<int, double>& jacintcellmap,
        std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief evaluate weighted Gap entries at GP (quad-pwlin)

     */
    void inline GP_3D_G_Quad_pwlin(MORTAR::Element& sele, MORTAR::IntElement& sintele,
        MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseVector& lmintval,
        CORE::LINALG::SerialDenseMatrix& scoord, CORE::LINALG::SerialDenseMatrix& mcoord,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
        double* gap, double* gpn, double* lengthn, double& jac, double& wgt,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        CORE::GEN::pairedvector<int, double>& dgapgp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit);

    /*!
     \brief evaluate weighted Gap entries at GP

     */
    void GP_G_Lin(int& iter, MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, double& gap, double* gpn, double& jac,
        double& wgt, CORE::GEN::pairedvector<int, double>& dgapgp,
        CORE::GEN::pairedvector<int, double>& jacintcellmap,
        std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief evaluate weighted Gap entries at GP (quad)

     */
    void inline GP_3D_G_Quad_Lin(int& iter, MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& svalmod,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, double& gap, double* gpn, double& jac,
        double& wgt, bool& duallin, const CORE::GEN::pairedvector<int, double>& dgapgp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dpsxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        bool dualquad3d);

    /*!
     \brief evaluate weighted Gap entries at GP (quad)

     */
    void inline GP_3D_G_Quad_pwlin_Lin(int& iter, MORTAR::IntElement& sintele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmintval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmintderiv,
        double& gap, double* gpn, double& jac, double& wgt,
        const CORE::GEN::pairedvector<int, double>& dgapgp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp);

    /*!
     \brief evaluate and lin slipincr at GP

     */
    void inline GP_2D_SlipIncr(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, double& jac, double& wgt, double* jumpvalv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        CORE::GEN::pairedvector<int, double>& dslipgp, int& linsize);

    /*!
     \brief evaluate and lin slipincr at GP

     */
    void inline GP_3D_SlipIncr(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, double& jac, double& wgt, double* jumpvalv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dslipgp);

    /*!
     \brief evaluate and lin slipincr at GP at node

     */
    void inline GP_2D_SlipIncr_Lin(int& iter, MORTAR::Element& sele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double& wgt, double* jumpvalv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, double>& dslipgp,
        const CORE::GEN::pairedvector<int, double>& derivjac,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    void inline GP_3D_SlipIncr_Lin(int& iter, MORTAR::Element& sele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double& wgt, double* jumpvalv,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dslipgp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);
    /*!
     \brief evaluate  T and E matrix

     */
    void inline GP_TE(MORTAR::Element& sele, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& sval, double& jac, double& wgt, double* jumpval);

    /*!
     \brief evaluate  T and E matrix

     */
    void inline GP_TE_Master(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseVector& lm2val,
        CORE::LINALG::SerialDenseVector& mval, double& jac, double& wgt, double* jumpval,
        const Epetra_Comm& comm);

    /*!
     \brief evaluate Lin T and E matrix

     */
    void inline GP_2D_TE_Lin(int& iter, MORTAR::Element& sele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double& wgt, double* jumpval,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, double>& derivjac,
        const CORE::GEN::pairedvector<int, double>& dsliptmatrixgp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief evaluate Lin T and E matrix

     */
    void inline GP_2D_TE_Master_Lin(int& iter,  // like k
        MORTAR::Element& sele, MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& dsxideta, double& dxdsxi, double& dxdsxidsxi, double& wgt, double* jumpval,
        const CORE::GEN::pairedvector<int, double>& dsxigp,
        const CORE::GEN::pairedvector<int, double>& dmxigp,
        const CORE::GEN::pairedvector<int, double>& derivjac,
        const CORE::GEN::pairedvector<int, double>& dsliptmatrixgp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& ximaps,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        const Epetra_Comm& comm);

    /*!
     \brief evaluate Lin T and E matrix

     */
    void inline GP_3D_TE_Lin(int& iter, MORTAR::Element& sele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double& wgt, double* jumpval,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const CORE::GEN::pairedvector<int, double>& dsliptmatrixgp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief evaluate Lin T and E matrix (Master)

     */
    void inline GP_3D_TE_Master_Lin(int& iter, MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseVector& lm2val,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, CORE::LINALG::SerialDenseMatrix& lm2deriv,
        double& jac, double& wgt, double* jumpval,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const CORE::GEN::pairedvector<int, double>& dsliptmatrixgp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dual2map,
        const Epetra_Comm& comm);

    /*!
     \brief evaluate wear + lin at GP

     */
    void inline GP_2D_Wear(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& lmderiv,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> lagmult, double* gpn, double& jac,
        double& wgt, double* jumpval, double* wearval,
        CORE::GEN::pairedvector<int, double>& dsliptmatrixgp,
        CORE::GEN::pairedvector<int, double>& dweargp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief evaluate wear + lin at GP

     */
    void inline GP_3D_Wear(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& lmderiv,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> lagmult, double* gpn, double& jac,
        double& wgt, double* jumpval, double* wearval,
        CORE::GEN::pairedvector<int, double>& dsliptmatrixgp,
        CORE::GEN::pairedvector<int, double>& dweargp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin weighted wear at GP

     */
    void inline GP_2D_Wear_Lin(int& iter, MORTAR::Element& sele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double* gpn, double& wgt, double& wearval, double* jumpval,
        const CORE::GEN::pairedvector<int, double>& dweargp,
        const CORE::GEN::pairedvector<int, double>& derivjac,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
     \brief lin weighted wear at GP

     */
    void inline GP_3D_Wear_Lin(int& iter, MORTAR::Element& sele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        double& jac, double* gpn, double& wgt, double& wearval, double* jumpval,
        const CORE::GEN::pairedvector<int, double>& dweargp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);


    /*!
    \brief evaluate scalar normal coupling condition for poro no penetration entries at GP
    (poro-contact)

    */
    void inline GP_NCOUP_DERIV(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, double* ncoup, double* gpn, double& jac,
        double& wgt, double* gpcoord,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        std::map<int, double>& dncoupgp, std::map<int, double>& dvelncoupgp,
        std::map<int, double>& dpresncoupgp,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, bool quadratic,
        int nintrow = 0);

    /*!
    \brief evaluate weighted normal coupling entries at GP

    */
    void inline GP_NCOUP_LIN(int& iter, MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& mval,
        CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& lmderiv, double& ncoup, double* gpn, double& jac,
        double& wgt, const std::map<int, double>& dncoupgp,
        const std::map<int, double>& dvelncoupgp, const std::map<int, double>& dpresncoupgp,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxigp,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap);

    /*!
    \brief Calculate Determinate of the Deformation Gradient at GP

    */
    double DetDeformationGradient(
        MORTAR::Element& sele, double& wgt, double* gpcoord, std::map<int, double>& JLin);

    /*!
    \brief Templated Calculate Determinate of the Deformation Gradient at GP

    */
    template <CORE::FE::CellType parentdistype, int dim>
    double TDetDeformationGradient(
        MORTAR::Element& sele, double& wgt, double* gpcoord, std::map<int, double>& JLin);

    /*!
     \brief Return the Wear shape fcn type (wear weighting...)

     */
    INPAR::WEAR::WearShape WearShapeFcn() { return wearshapefcn_; }

    /*!
     \brief Return type of wear surface definition

     */
    INPAR::WEAR::WearSide WearSide() { return wearside_; }

    /*!
     \brief Return type of wear algorithm

     */
    INPAR::WEAR::WearType WearType() { return weartype_; }

    /*!
     \brief Return the LM shape fcn type

     */
    INPAR::MORTAR::ShapeFcn ShapeFcn() { return shapefcn_; }

    /*!
     \brief Return the LM interpolation / testing type for quadratic FE

     */
    INPAR::MORTAR::LagMultQuad LagMultQuad() { return lagmultquad_; }
    //@}

    //! containing contact input parameters
    Teuchos::ParameterList& imortar_;
    //! communicator
    const Epetra_Comm& Comm_;

    //! number of Gauss points
    int ngp_;
    //! Gauss point coordinates
    CORE::LINALG::SerialDenseMatrix coords_;
    //! Gauss point weights
    std::vector<double> weights_;
    //! dimension of problem (2D or 3D)
    int dim_;

    // inputs from parameter list
    //! lm shape function type
    INPAR::MORTAR::ShapeFcn shapefcn_;
    //! type of lm interpolation for quadr. FE
    INPAR::MORTAR::LagMultQuad lagmultquad_;
    //! gp-wise evaluated slip increment
    bool gpslip_;
    //! contact algorithm
    INPAR::MORTAR::AlgorithmType algo_;
    //! solution stratety
    INPAR::CONTACT::SolvingStrategy stype_;
    //! flag for closest point normal -> change in linsize
    bool cppnormal_;

    // wear inputs from parameter list
    //! type of wear law
    INPAR::WEAR::WearLaw wearlaw_;
    //! flag for implicit wear algorithm
    bool wearimpl_;
    //! definition of wear surface
    INPAR::WEAR::WearSide wearside_;
    //! definition of contact wear algorithm
    INPAR::WEAR::WearType weartype_;
    //! type of wear shape function
    INPAR::WEAR::WearShape wearshapefcn_;
    //! flag for steady state wear
    bool sswear_;
    //! wear coefficient
    double wearcoeff_;
    //! wear coefficient master
    double wearcoeffm_;
    //! fixed slip for steady state wear
    double ssslip_;

    //! flag for non-smooth contact
    bool nonsmooth_;
    //! flag is true if (self) contact surface is non-smooth
    const bool nonsmoothselfcontactsurface_;

   private:
    //! integration type from the parameter-list
    INPAR::MORTAR::IntType integrationtype_;
  };  // class Integrator
}  // namespace CONTACT


BACI_NAMESPACE_CLOSE

#endif