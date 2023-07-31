/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beam3_kirchhoff.H"

#include "baci_beam3_triad_interpolation_local_rotation_vectors.H"

#include "baci_beaminteraction_periodic_boundingbox.H"

#include "baci_structure_new_elements_paramsinterface.H"

// Todo @grill: check for obsolete header inclusions
#include "baci_utils_fad.H"
#include "baci_lib_discret.H"
#include "baci_lib_utils.H"
#include "baci_utils_exceptions.H"
#include "baci_lib_globalproblem.H"
#include "baci_inpar_validparameters.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_discretization_fem_general_largerotations.H"
#include "baci_lib_linedefinition.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_beaminteraction_periodic_boundingbox.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3kType DRT::ELEMENTS::Beam3kType::instance_;

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3kType& DRT::ELEMENTS::Beam3kType::Instance() { return instance_; }

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::Beam3kType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Beam3k* object = new DRT::ELEMENTS::Beam3k(-1, -1);
  object->Unpack(data);
  return object;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3kType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "BEAM3K")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3k(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3kType::Create(const int id, const int owner)

{
  return Teuchos::rcp(new Beam3k(id, owner));
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3kType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  dserror("method 'NodalBlockInformation' not implemented for element type beam3k!");
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::SerialDenseMatrix<int, double> DRT::ELEMENTS::Beam3kType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Teuchos::SerialDenseMatrix<int, double> nullspace;
  dserror("method ComputeNullSpace not implemented for element type beam3k!");
  return nullspace;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3kType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["BEAM3K"];

  defs["LINE2"]
      .AddIntVector("LINE2", 2)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 6)
      .AddOptionalTag("FAD");

  defs["LINE3"]
      .AddIntVector("LINE3", 3)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 9)
      .AddOptionalTag("FAD");

  defs["LINE4"]
      .AddIntVector("LINE4", 4)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 12)
      .AddOptionalTag("FAD");
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 01/16|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3kType::Initialize(DRT::Discretization& dis)
{
  // setting up geometric variables for Beam3k elements
  for (int num = 0; num < dis.NumMyColElements(); ++num)
  {
    // in case that current element is not a Beam3k element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    // if we get so far current element is a Beam3k element and  we get a pointer at it
    DRT::ELEMENTS::Beam3k* currele = dynamic_cast<DRT::ELEMENTS::Beam3k*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3k* failed");

    // reference node position
    std::vector<CORE::LINALG::Matrix<3, 1>> xrefe;

    const int nnode = currele->NumNode();

    // resize xrefe for the number of nodes to store
    xrefe.resize(nnode);

    // the next section is needed in case of periodic boundary conditions and a shifted
    // configuration (i.e. elements cut by the periodic boundary) in the input file
    Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> periodic_boundingbox =
        Teuchos::rcp(new CORE::GEO::MESHFREE::BoundingBox());
    periodic_boundingbox->Init();  // no Setup() call needed here

    std::vector<double> disp_shift;
    int numdof = currele->NumDofPerNode(*(currele->Nodes()[0]));
    disp_shift.resize(numdof * nnode);
    for (unsigned int i = 0; i < disp_shift.size(); ++i) disp_shift[i] = 0.0;
    if (periodic_boundingbox->HavePBC())
      currele->UnShiftNodePosition(disp_shift, *periodic_boundingbox);

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int node = 0; node < nnode; ++node)  // element has k nodes
        for (int dof = 0; dof < 3; ++dof)       // element node has three coordinates x1, x2 and x3
        {
          xrefe[node](dof) = currele->Nodes()[node]->X()[dof] + disp_shift[node * numdof + dof];
        }
    }

    // Set up all geometrical (triads, curvatures, jacobians etc.) quantities describing the
    // (initial) reference geometry
    currele->SetUpReferenceGeometry(xrefe);
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::Beam3k(int id, int owner)
    : DRT::ELEMENTS::Beam3Base(id, owner),
      useFAD_(false),
      isinit_(false),
      T_(0),
      theta0_(0),
      Qrefconv_(0),
      Qrefnew_(0),
      K0_(0),
      length_(0.0),
      jacobi_(0.0),
      jacobi2_(0.0),
      jacobi_cp_(0.0),
      jacobi2_cp_(0.0),
      rotvec_(false),
      weakkirchhoff_(false),
      Eint_(0.0),
      Ekin_(0.0),
      Qconvmass_(0),
      Qnewmass_(0),
      wconvmass_(0),
      wnewmass_(0),
      aconvmass_(0),
      anewmass_(0),
      amodconvmass_(0),
      amodnewmass_(0),
      rttconvmass_(0),
      rttnewmass_(0),
      rttmodconvmass_(0),
      rttmodnewmass_(0),
      rtconvmass_(0),
      rtnewmass_(0),
      rconvmass_(0),
      rnewmass_(0),
      axial_strain_GP_(0),
      twist_GP_(0),
      curvature_2_GP_(0),
      curvature_3_GP_(0),
      axial_force_GP_(0),
      torque_GP_(0),
      bending_moment_2_GP_(0),
      bending_moment_3_GP_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::Beam3k(const DRT::ELEMENTS::Beam3k& old)
    : DRT::ELEMENTS::Beam3Base(old),
      useFAD_(old.useFAD_),
      isinit_(old.isinit_),
      T_(old.T_),
      theta0_(old.theta0_),
      Qrefconv_(old.Qrefconv_),
      Qrefnew_(old.Qrefnew_),
      K0_(old.K0_),
      length_(old.length_),
      jacobi_(old.jacobi_),
      jacobi2_(old.jacobi2_),
      jacobi_cp_(old.jacobi_cp_),
      jacobi2_cp_(old.jacobi_cp_),
      rotvec_(old.rotvec_),
      weakkirchhoff_(old.weakkirchhoff_),
      Eint_(old.Eint_),
      Ekin_(old.Ekin_),
      Qconvmass_(old.Qconvmass_),
      Qnewmass_(old.Qnewmass_),
      wconvmass_(old.wconvmass_),
      wnewmass_(old.wnewmass_),
      aconvmass_(old.aconvmass_),
      anewmass_(old.anewmass_),
      amodconvmass_(old.amodconvmass_),
      amodnewmass_(old.amodnewmass_),
      rttconvmass_(old.rttconvmass_),
      rttnewmass_(old.rttnewmass_),
      rttmodconvmass_(old.rttmodconvmass_),
      rttmodnewmass_(old.rttmodnewmass_),
      rtconvmass_(old.rtconvmass_),
      rtnewmass_(old.rtnewmass_),
      rconvmass_(old.rconvmass_),
      rnewmass_(old.rnewmass_),
      axial_strain_GP_(old.axial_strain_GP_),
      twist_GP_(old.twist_GP_),
      curvature_2_GP_(old.curvature_2_GP_),
      curvature_3_GP_(old.curvature_3_GP_),
      axial_force_GP_(old.axial_force_GP_),
      torque_GP_(old.torque_GP_),
      bending_moment_2_GP_(old.bending_moment_2_GP_),
      bending_moment_3_GP_(old.bending_moment_3_GP_)
{
  return;
}
/*--------------------------------------------------------------------------------*
 |  Deep copy this instance of Beam3k and return pointer to it (public) |
 |                                                                    meier 05/12 |
 *--------------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3k::Clone() const
{
  DRT::ELEMENTS::Beam3k* newelement = new DRT::ELEMENTS::Beam3k(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::~Beam3k() { return; }

/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Print(std::ostream& os) const { return; }

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3k::Shape() const
{
  int numnodes = NumNode();
  switch (numnodes)
  {
    case 2:
      return line2;
      break;
    case 3:
      return line3;
      break;
    case 4:
      return line4;
      break;
    default:
      dserror("Only Line2, Line3 and Line4 elements are implemented.");
      break;
  }

  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           meier 05/12/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Beam3Base::Pack(data);

  // add all class variables
  AddtoPack(data, useFAD_);
  AddtoPack(data, isinit_);
  AddtoPack<3, 1>(data, Tref_);
  AddtoPack<3, 1>(data, T_);
  AddtoPack<3, 1>(data, theta0_);
  AddtoPack<4, 1>(data, Qrefconv_);
  AddtoPack<4, 1>(data, Qrefnew_);
  AddtoPack<3, 1>(data, K0_);
  AddtoPack(data, length_);
  AddtoPack(data, jacobi_);
  AddtoPack(data, jacobi2_);
  AddtoPack(data, jacobi_cp_);
  AddtoPack(data, jacobi2_cp_);
  AddtoPack(data, rotvec_);
  AddtoPack(data, weakkirchhoff_);
  AddtoPack(data, Eint_);
  AddtoPack(data, Ekin_);
  AddtoPack<4, 1>(data, Qconvmass_);
  AddtoPack<4, 1>(data, Qnewmass_);
  AddtoPack<3, 1>(data, wconvmass_);
  AddtoPack<3, 1>(data, wnewmass_);
  AddtoPack<3, 1>(data, aconvmass_);
  AddtoPack<3, 1>(data, anewmass_);
  AddtoPack<3, 1>(data, amodconvmass_);
  AddtoPack<3, 1>(data, amodnewmass_);
  AddtoPack<3, 1>(data, rttconvmass_);
  AddtoPack<3, 1>(data, rttnewmass_);
  AddtoPack<3, 1>(data, rttmodconvmass_);
  AddtoPack<3, 1>(data, rttmodnewmass_);
  AddtoPack<3, 1>(data, rtconvmass_);
  AddtoPack<3, 1>(data, rtnewmass_);
  AddtoPack<3, 1>(data, rconvmass_);
  AddtoPack<3, 1>(data, rnewmass_);
  AddtoPack(data, axial_strain_GP_);
  AddtoPack(data, twist_GP_);
  AddtoPack(data, curvature_2_GP_);
  AddtoPack(data, curvature_3_GP_);
  AddtoPack(data, axial_force_GP_);
  AddtoPack(data, torque_GP_);
  AddtoPack(data, bending_moment_2_GP_);
  AddtoPack(data, bending_moment_3_GP_);

  return;
}
/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Beam3Base::Unpack(basedata);

  // extract all class variables of beam3k element
  useFAD_ = ExtractInt(position, data);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack<3, 1>(position, data, Tref_);
  ExtractfromPack<3, 1>(position, data, T_);
  ExtractfromPack<3, 1>(position, data, theta0_);
  ExtractfromPack<4, 1>(position, data, Qrefconv_);
  ExtractfromPack<4, 1>(position, data, Qrefnew_);
  ExtractfromPack<3, 1>(position, data, K0_);
  ExtractfromPack(position, data, length_);
  ExtractfromPack(position, data, jacobi_);
  ExtractfromPack(position, data, jacobi2_);
  ExtractfromPack(position, data, jacobi_cp_);
  ExtractfromPack(position, data, jacobi2_cp_);
  rotvec_ = ExtractInt(position, data);
  weakkirchhoff_ = ExtractInt(position, data);
  ExtractfromPack(position, data, Eint_);
  ExtractfromPack(position, data, Ekin_);
  ExtractfromPack<4, 1>(position, data, Qconvmass_);
  ExtractfromPack<4, 1>(position, data, Qnewmass_);
  ExtractfromPack<3, 1>(position, data, wconvmass_);
  ExtractfromPack<3, 1>(position, data, wnewmass_);
  ExtractfromPack<3, 1>(position, data, aconvmass_);
  ExtractfromPack<3, 1>(position, data, anewmass_);
  ExtractfromPack<3, 1>(position, data, amodconvmass_);
  ExtractfromPack<3, 1>(position, data, amodnewmass_);
  ExtractfromPack<3, 1>(position, data, rttconvmass_);
  ExtractfromPack<3, 1>(position, data, rttnewmass_);
  ExtractfromPack<3, 1>(position, data, rttmodconvmass_);
  ExtractfromPack<3, 1>(position, data, rttmodnewmass_);
  ExtractfromPack<3, 1>(position, data, rtconvmass_);
  ExtractfromPack<3, 1>(position, data, rtnewmass_);
  ExtractfromPack<3, 1>(position, data, rconvmass_);
  ExtractfromPack<3, 1>(position, data, rnewmass_);
  ExtractfromPack(position, data, axial_strain_GP_);
  ExtractfromPack(position, data, twist_GP_);
  ExtractfromPack(position, data, curvature_2_GP_);
  ExtractfromPack(position, data, curvature_3_GP_);
  ExtractfromPack(position, data, axial_force_GP_);
  ExtractfromPack(position, data, torque_GP_);
  ExtractfromPack(position, data, bending_moment_2_GP_);
  ExtractfromPack(position, data, bending_moment_3_GP_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          meier 05/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Beam3k::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpInitialRotations(const std::vector<double>& nodal_thetas)
{
  theta0_.resize(BEAM3K_COLLOCATION_POINTS);
  for (int i = 0; i < BEAM3K_COLLOCATION_POINTS; i++)
  {
    for (int j = 0; j < 3; j++) (theta0_[i])(j) = nodal_thetas[3 * i + j];

    // Shift angles by 2PI in case these angles are not in the interval [-PI,PI].
    if (theta0_[i].Norm2() > M_PI)
    {
      CORE::LINALG::Matrix<4, 1> Q(true);
      CORE::LARGEROTATIONS::angletoquaternion(theta0_[i], Q);
      CORE::LARGEROTATIONS::quaterniontoangle(Q, theta0_[i]);
    }
  }
}

/*----------------------------------------------------------------------*
 | sets up geometric data from current nodal position as reference
 | position; this method can be used by the register class or when ever
 | a new beam element is generated for which some reference configuration
 | has to be stored; prerequesite for applying this method is that the
 | element nodes are already known (public)                   meier 01/16|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometry(
    const std::vector<CORE::LINALG::Matrix<3, 1>>& xrefe, const bool secondinit)
{
  if (weakkirchhoff_)
    SetUpReferenceGeometryWK(xrefe, secondinit);
  else
    SetUpReferenceGeometrySK(xrefe, secondinit);
}

/*--------------------------------------------------------------------------------------------*
 |  Set up the reference geometry of the case of a weak Kirchhoff constraint     meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometryWK(
    const std::vector<CORE::LINALG::Matrix<3, 1>>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe and rotrefe; note: the
   *isinit_ flag is important for avoiding reinitialization upon restart. However, it should be
   *possible to conduct a second initilization in principle (e.g. for periodic boundary conditions*/

  // TODO: It turned out that the approximation of initially curved centerlines via third order
  // Hermite polynomials by simply setting the nodal positions and tangents at the boundary nodes to
  // the corresponding values of the analytical geometry at these nodes yields a worse approximation
  // than that following from a third order Lagrange polynomial interpolation based on nodal
  // positions that coincide with the analytical values at the four element nodes. This worse
  // approximation appears in terms of a larger error in the jacobian, the derivative of the
  // relative angle theta_s and finally of the initial curvature K0_. For strongly curved initial
  // configurations where the error in the initial geometry representation might dominate the
  // discretization error it can be useful to consider strategies which yield in a better Hermite
  // approximation of the initial geometry e.g. by determining the initial nodal values (or the
  // shape function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if (!isinit_ || secondinit)
  {
    const int nnode = 2;  // number of nodes

    // Calculate the (initial reference triads) = (initial material triads) at the CPs out of the
    // angles theta0_. So far the initial value for the relative angle is set to zero, i.e. material
    // coordinate system and reference system in the reference configuration coincidence.
    std::vector<CORE::LINALG::Matrix<3, 3>> Gref(BEAM3K_COLLOCATION_POINTS);
    for (int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      Gref[node].Clear();
      CORE::LARGEROTATIONS::angletotriad(theta0_[node], Gref[node]);
    }

    Tref_.resize(2);
    T_.resize(2);
    // write initial nodal tangents in extra vector
    for (int i = 0; i < 3; i++)
    {
      (Tref_[0])(i) = (Gref[0])(i, 0);
      (Tref_[1])(i) = (Gref[1])(i, 0);
      (T_[0])(i) = (Gref[0])(i, 0);
      (T_[1])(i) = (Gref[1])(i, 0);
    }

    // Get integration points for exact integration
    CORE::DRT::UTILS::IntegrationPoints1D gausspoints =
        CORE::DRT::UTILS::IntegrationPoints1D(MYGAUSSRULEBEAM3K);

    // Vector holding angle theta of triads
    std::vector<CORE::LINALG::Matrix<3, 1>> theta_cp;
    theta_cp.resize(BEAM3K_COLLOCATION_POINTS);
    CORE::LINALG::Matrix<3, 1> theta(true);
    CORE::LINALG::Matrix<3, 1> theta_s(true);
    CORE::LINALG::Matrix<3, 3> triad_mat(true);  // material triad at gp

    // Resize vectors for storage of time integration quantities
    ResizeClassVariables(gausspoints.nquad);

    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    axial_strain_GP_.resize(gausspoints.nquad);
    std::fill(axial_strain_GP_.begin(), axial_strain_GP_.end(), 0.0);
    twist_GP_.resize(gausspoints.nquad);
    std::fill(twist_GP_.begin(), twist_GP_.end(), 0.0);
    curvature_2_GP_.resize(gausspoints.nquad);
    std::fill(curvature_2_GP_.begin(), curvature_2_GP_.end(), 0.0);
    curvature_3_GP_.resize(gausspoints.nquad);
    std::fill(curvature_3_GP_.begin(), curvature_3_GP_.end(), 0.0);

    axial_force_GP_.resize(gausspoints.nquad);
    std::fill(axial_force_GP_.begin(), axial_force_GP_.end(), 0.0);
    torque_GP_.resize(gausspoints.nquad);
    std::fill(torque_GP_.begin(), torque_GP_.end(), 0.0);
    bending_moment_2_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_2_GP_.begin(), bending_moment_2_GP_.end(), 0.0);
    bending_moment_3_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_3_GP_.begin(), bending_moment_3_GP_.end(), 0.0);

    // calculate the length of the element via Newton iteration
    CORE::LINALG::Matrix<12, 1> pos_ref_centerline;
    AddRefValuesDispCenterline<2, 2>(pos_ref_centerline);
    length_ = CalcReflength<2, 2>(pos_ref_centerline);

    // Matrices to store the function values of the Lagrange shape functions used to interpolate
    // theta
    CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i;
    CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i_xi;

    // Matrices to store the (derivative of) the Hermite shape functions
    CORE::LINALG::Matrix<1, 2 * nnode> N_i;
    CORE::LINALG::Matrix<1, 2 * nnode> N_i_xi;

    // Matrices to store r and dr/dxi
    CORE::LINALG::Matrix<3, 1> r;
    CORE::LINALG::Matrix<3, 1> r_xi;

    // storage index for collocation point
    unsigned int ind = 0;

    // Calculate initial material triads at the collocation points
    for (unsigned int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      // colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi = (double)node / (BEAM3K_COLLOCATION_POINTS - 1) * 2 - 1.0;

      // Get values of shape functions
      N_i_xi.Clear();
      CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);

      // Determine storage position for the node colpt
      ind = CORE::LARGEROTATIONS::NumberingTrafo(node + 1, BEAM3K_COLLOCATION_POINTS);

      // current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();

      for (int i = 0; i < 3; i++)
      {
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + Tref_[0](i) * N_i_xi(1) +
                   Tref_[1](i) * N_i_xi(3);
      }

      jacobi_cp_[ind] = r_xi.Norm2();

      // rotate (initial reference triad) = (initial material triad) at the interior CPs on
      // tangential line resulting from the Hermite interpolation. This is necessary for initially
      // curved geometries for which the Hermite interpolation does not deliver the exact tangent
      // values of the analytical representation of the initial curve at the CPs.
      if (ind > 1)  // only for internal CPs
      {
        CORE::LINALG::Matrix<3, 3> G_aux(true);
        CORE::LARGEROTATIONS::CalculateSRTriads<double>(r_xi, Gref[ind], G_aux);
        // rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind] = G_aux;
        CORE::LARGEROTATIONS::triadtoquaternion(G_aux, Qrefconv_[ind]);
        CORE::LARGEROTATIONS::quaterniontoangle(Qrefconv_[ind], theta0_[ind]);
      }
      else
      {
        CORE::LARGEROTATIONS::triadtoquaternion(Gref[ind], Qrefconv_[ind]);
      }
      Qrefnew_[ind] = Qrefconv_[ind];
    }

    // SETUP INTERPOLATION via calculation of difference angle
    for (int colpt = 0; colpt < BEAM3K_COLLOCATION_POINTS; colpt++)
    {
      theta_cp[colpt].Clear();
      CORE::LARGEROTATIONS::triadtoangleright(theta_cp[colpt], Gref[REFERENCE_NODE], Gref[colpt]);
    }

    // Loop through all GPs and computation of all relevant values at each gp
    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      // Get values of shape functions
      L_i.Clear();
      L_i_xi.Clear();
      N_i_xi.Clear();
      N_i.Clear();
      CORE::DRT::UTILS::shape_function_1D(L_i, xi, Shape());
      CORE::DRT::UTILS::shape_function_1D_deriv1(L_i_xi, xi, Shape());
      CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
      CORE::DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);

      // current value of derivatives at GP (derivatives in xi!)
      r.Clear();
      r_xi.Clear();
      theta.Clear();
      theta_s.Clear();

      for (int i = 0; i < 3; i++)
      {
        r(i) += xrefe[0](i) * N_i(0) + xrefe[1](i) * N_i(2) + Tref_[0](i) * N_i(1) +
                Tref_[1](i) * N_i(3);
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + Tref_[0](i) * N_i_xi(1) +
                   Tref_[1](i) * N_i_xi(3);
      }

      // calculate jacobi jacobi_=|r'_0|
      jacobi_[numgp] = r_xi.Norm2();

      // calculate interpolated angle
      for (int i = 0; i < BEAM3K_COLLOCATION_POINTS; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          theta(j) += L_i(i) * theta_cp[i](j);
          theta_s(j) += L_i_xi(i) * theta_cp[i](j) / jacobi_[numgp];
        }
      }
      computestrain(theta, theta_s, K0_[numgp]);

      CORE::LINALG::Matrix<3, 3> temp_triad(true);
      CORE::LARGEROTATIONS::angletotriad(theta, temp_triad);
      triad_mat.Multiply(Gref[REFERENCE_NODE], temp_triad);
      SetInitialDynamicClassVariables(numgp, triad_mat, r);
    }

    isinit_ = true;
  }
}

/*--------------------------------------------------------------------------------------------*
 |  Set up the reference geometry of the case of a strong Kirchhoff constraint     meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometrySK(
    const std::vector<CORE::LINALG::Matrix<3, 1>>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe and rotrefe; note: the
   *isinit_ flag is important for avoiding reinitialization upon restart. However, it should be
   *possible to conduct a second initilization in principle (e.g. for periodic boundary conditions*/

  // TODO: It turned out that the approximation of initially curved centerlines via third order
  // Hermite polynomials by simply setting the nodal positions and tangents at the boundary nodes to
  // the corresponding values of the analytical geometry at these nodes yields a worse approximation
  // than that following from a third order Lagrange polynomial interpolation based on nodal
  // positions that coincide with the analytical values at the four element nodes. This worse
  // approximation appears in terms of a larger error in the jacobian, the derivative of the
  // relative angle theta_s and finally of the initial curvature K0_. For strongly curved initial
  // configurations where the error in the initial geometry representation might dominate the
  // discretization error it can be useful to consider strategies which yield in a better Hermite
  // approximation of the initial geometry e.g. by determining the initial nodal values (or the
  // shape function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if (!isinit_ || secondinit)
  {
    const int nnode = 2;  // number of nodes

    // Calculate the (initial reference triads) = (initial material triads) at the CPs out of the
    // angles theta0_. So far the initial value for the relative angle is set to zero, i.e. material
    // coordinate system and reference system in the reference configuration coincidence.
    std::vector<CORE::LINALG::Matrix<3, 3>> Gref(BEAM3K_COLLOCATION_POINTS);
    for (int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      Gref[node].Clear();
      CORE::LARGEROTATIONS::angletotriad(theta0_[node], Gref[node]);
    }

    Tref_.resize(2);
    T_.resize(2);
    // write initial nodal tangents in extra vector
    for (int i = 0; i < 3; i++)
    {
      (Tref_[0])(i) = (Gref[0])(i, 0);
      (Tref_[1])(i) = (Gref[1])(i, 0);
      (T_[0])(i) = (Gref[0])(i, 0);
      (T_[1])(i) = (Gref[1])(i, 0);
    }

    // Get integration points for exact integration
    CORE::DRT::UTILS::IntegrationPoints1D gausspoints =
        CORE::DRT::UTILS::IntegrationPoints1D(MYGAUSSRULEBEAM3K);

    // Vector holding angle theta of triads
    std::vector<double> phi_cp;
    phi_cp.resize(BEAM3K_COLLOCATION_POINTS);
    double phi(true);
    double phi_s(true);
    CORE::LINALG::Matrix<3, 3> triad_mat(true);  // material triad at gp

    // Resize vectors for storage of time integration quantities
    ResizeClassVariables(gausspoints.nquad);

    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    axial_strain_GP_.resize(gausspoints.nquad);
    std::fill(axial_strain_GP_.begin(), axial_strain_GP_.end(), 0.0);
    twist_GP_.resize(gausspoints.nquad);
    std::fill(twist_GP_.begin(), twist_GP_.end(), 0.0);
    curvature_2_GP_.resize(gausspoints.nquad);
    std::fill(curvature_2_GP_.begin(), curvature_2_GP_.end(), 0.0);
    curvature_3_GP_.resize(gausspoints.nquad);
    std::fill(curvature_3_GP_.begin(), curvature_3_GP_.end(), 0.0);

    axial_force_GP_.resize(gausspoints.nquad);
    std::fill(axial_force_GP_.begin(), axial_force_GP_.end(), 0.0);
    torque_GP_.resize(gausspoints.nquad);
    std::fill(torque_GP_.begin(), torque_GP_.end(), 0.0);
    bending_moment_2_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_2_GP_.begin(), bending_moment_2_GP_.end(), 0.0);
    bending_moment_3_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_3_GP_.begin(), bending_moment_3_GP_.end(), 0.0);

    // calculate the length of the element via Newton iteration
    CORE::LINALG::Matrix<12, 1> disp_refe_centerline;
    AddRefValuesDispCenterline<2, 2>(disp_refe_centerline);
    length_ = CalcReflength<2, 2>(disp_refe_centerline);

    // Matrices to store the function values of the Lagrange shape functions used to interpolate
    // theta
    CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i;
    CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i_xi;

    // Matrices to store the (derivative of) the Hermite shape functions
    CORE::LINALG::Matrix<1, 2 * nnode> N_i;
    CORE::LINALG::Matrix<1, 2 * nnode> N_i_xi;
    CORE::LINALG::Matrix<1, 2 * nnode> N_i_xixi;

    // Matrices to store r, dr/dxi and d^2r/dxi^2
    CORE::LINALG::Matrix<3, 1> r;
    CORE::LINALG::Matrix<3, 1> r_xi;
    CORE::LINALG::Matrix<3, 1> r_xixi;
    CORE::LINALG::Matrix<3, 1> r_s;
    CORE::LINALG::Matrix<3, 1> r_ss;
    // centerline curvature
    CORE::LINALG::Matrix<3, 1> kappacl;

    // storage index for collocation point
    unsigned int ind = 0;

    // Calculate initial material triads at the collocation points
    for (unsigned int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      // colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi = (double)node / (BEAM3K_COLLOCATION_POINTS - 1) * 2 - 1.0;

      // Get values of shape functions
      L_i.Clear();
      N_i_xi.Clear();
      N_i_xixi.Clear();
      CORE::DRT::UTILS::shape_function_1D(L_i, xi, Shape());
      CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
      CORE::DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi, xi, length_, line2);

      // Determine storage position for the node colpt
      ind = CORE::LARGEROTATIONS::NumberingTrafo(node + 1, BEAM3K_COLLOCATION_POINTS);

      // current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();
      r_xixi.Clear();

      for (int i = 0; i < 3; i++)
      {
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + Tref_[0](i) * N_i_xi(1) +
                   Tref_[1](i) * N_i_xi(3);
        r_xixi(i) += xrefe[0](i) * N_i_xixi(0) + xrefe[1](i) * N_i_xixi(2) +
                     Tref_[0](i) * N_i_xixi(1) + Tref_[1](i) * N_i_xixi(3);
      }

      // calculate jacobi_=||r'_0|| and jacobi2_=r'_0^T r''_0
      jacobi_cp_[ind] = r_xi.Norm2();
      jacobi2_cp_[ind] = r_xi.Dot(r_xixi);

      // rotate (initial reference triad) = (initial material triad) at the interior CPs on
      // tangential line resulting from the Hermite interpolation. This is necessary for initially
      // curved geometries for which the Hermite interpolation does not deliver the exact tangent
      // values of the analytical representation of the initial curve at the CPs.
      if (ind > 1)  // only for internal CPs
      {
        CORE::LINALG::Matrix<3, 3> G_aux(true);
        CORE::LARGEROTATIONS::CalculateSRTriads<double>(r_xi, Gref[ind], G_aux);
        // rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind] = G_aux;
        CORE::LARGEROTATIONS::triadtoquaternion(G_aux, Qrefconv_[ind]);
        CORE::LARGEROTATIONS::quaterniontoangle(Qrefconv_[ind], theta0_[ind]);
      }
      else
      {
        CORE::LARGEROTATIONS::triadtoquaternion(Gref[ind], Qrefconv_[ind]);
      }
      Qrefnew_[ind] = Qrefconv_[ind];
    }  //(int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)

    // SETUP INTERPOLATION via calculation of difference angle
    for (int colpt = 0; colpt < BEAM3K_COLLOCATION_POINTS; colpt++)
    {
      CORE::LINALG::Matrix<3, 3> Lambdabarref(true);
      CORE::LINALG::Matrix<3, 1> tangentref(true);
      CORE::LINALG::Matrix<3, 1> phivec(true);
      for (int i = 0; i < 3; i++)
      {
        tangentref(i) = Gref[colpt](i, 0);
      }
      CORE::LARGEROTATIONS::CalculateSRTriads<double>(
          tangentref, Gref[REFERENCE_NODE], Lambdabarref);
      CORE::LARGEROTATIONS::triadtoangleleft(phivec, Lambdabarref, Gref[colpt]);
      phi_cp[colpt] = 0.0;
      for (int i = 0; i < 3; i++)
      {
        phi_cp[colpt] += tangentref(i) * phivec(i);
      }
    }

    // Loop through all GPs and computation of all relevant values at each gp
    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      // Get values of shape functions
      L_i.Clear();
      L_i_xi.Clear();
      N_i.Clear();
      N_i_xi.Clear();
      N_i_xixi.Clear();

      CORE::DRT::UTILS::shape_function_1D(L_i, xi, Shape());
      CORE::DRT::UTILS::shape_function_1D_deriv1(L_i_xi, xi, Shape());
      CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
      CORE::DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi, xi, length_, line2);
      CORE::DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);

      // current value of derivatives at GP (derivatives in xi!)
      r.Clear();
      r_xi.Clear();
      r_xixi.Clear();
      r_s.Clear();
      r_ss.Clear();
      kappacl.Clear();
      phi = 0.0;
      phi_s = 0.0;

      for (int i = 0; i < 3; i++)
      {
        r(i) += xrefe[0](i) * N_i(0) + xrefe[1](i) * N_i(2) + Tref_[0](i) * N_i(1) +
                Tref_[1](i) * N_i(3);
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + Tref_[0](i) * N_i_xi(1) +
                   Tref_[1](i) * N_i_xi(3);
        r_xixi(i) += xrefe[0](i) * N_i_xixi(0) + xrefe[1](i) * N_i_xixi(2) +
                     Tref_[0](i) * N_i_xixi(1) + Tref_[1](i) * N_i_xixi(3);
      }

      // calculate jacobi_=||r'_0|| and jacobi2_=r'_0^T r''_0
      jacobi_[numgp] = r_xi.Norm2();
      jacobi2_[ind] = r_xi.Dot(r_xixi);

      // calculate interpolated angle
      for (int i = 0; i < BEAM3K_COLLOCATION_POINTS; i++)
      {
        phi += L_i(i) * phi_cp[i];
        phi_s += L_i_xi(i) * phi_cp[i] / jacobi_[numgp];
      }

      // calculate derivatives in s
      r_s = r_xi;
      r_s.Scale(1 / jacobi_[numgp]);
      for (int i = 0; i < 3; i++)
      {
        r_ss(i) = r_xixi(i) / pow(jacobi_[numgp], 2.0) -
                  r_xi(i) * jacobi2_[numgp] / pow(jacobi_[numgp], 4.0);
      }

      triad_mat.Clear();
      ComputeTriadSK(phi, r_s, Gref[REFERENCE_NODE], triad_mat);
      Calculate_clcurvature(r_s, r_ss, kappacl);

      computestrainSK(phi_s, kappacl, Gref[REFERENCE_NODE], triad_mat, K0_[numgp]);

      SetInitialDynamicClassVariables(numgp, triad_mat, r);
    }

    isinit_ = true;
  }  // if(!isinit_)

}  // DRT::ELEMENTS::Beam3k::SetUpReferenceGeometrySK()

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3k::GetJacobiFacAtXi(const double& xi) const
{
  const int nnode = 2;

  // Matrices to store the the Hermite shape function derivative values
  CORE::LINALG::Matrix<1, 2 * nnode> N_i_xi;
  CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);

  // jacobi = ds/dxi = ||r'_0||
  CORE::LINALG::Matrix<3, 1> r_xi;

  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    r_xi(dim) += Nodes()[0]->X()[dim] * N_i_xi(0) + Nodes()[1]->X()[dim] * N_i_xi(2) +
                 Tref_[0](dim) * N_i_xi(1) + Tref_[1](dim) * N_i_xi(3);
  }

  return r_xi.Norm2();
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetPosAtXi(
    CORE::LINALG::Matrix<3, 1>& pos, const double& xi, const std::vector<double>& disp) const
{
  /* we expect the (centerline) displacement state vector (positions and tangents at boundary nodes)
   * here; for flexibility, we also accept complete DoF vector of this element (2*7+1 DoFs) and do
   * the rest automatically here NOTE: the latter option is favorable in case of rotvec_==true: in
   * case of ROTVEC=true, we need 3 positions, 3 absolute rotvec DOFs AND the length of tangent
   * vector (7th DoF) */
  CORE::LINALG::Matrix<12, 1> disp_totlag_centerline(true);

  if (disp.size() == 15)
  {
    // in this case, we need to "add" reference values first, because if rotvec_==true,
    // we can extract tangent vectors only from total rotation vectors
    CORE::LINALG::Matrix<15, 1, double> disp_totlag(disp.data());
    AddRefValuesDisp<2, double>(disp_totlag);
    this->ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
        disp_totlag, disp_totlag_centerline);
  }
  else if (disp.size() == 12)
  {
    /* in this case, we expect the position and tangent DOF values for both boundary nodes;
     * for rotvec_==true, the tangents are NOT nodal DOFs, so they need to be pre-calculated
     * before calling this method */
    disp_totlag_centerline = CORE::LINALG::Matrix<12, 1>(disp.data());
    AddRefValuesDispCenterline<2, 2, double>(disp_totlag_centerline);
  }
  else
  {
    dserror(
        "size mismatch: expected either 12 values for disp_centerline or 15 for "
        "full element disp vector and got %d",
        disp.size());
  }

  Beam3Base::GetPosAtXi<2, 2>(pos, xi, disp_totlag_centerline);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetTriadAtXi(
    CORE::LINALG::Matrix<3, 3>& triad, const double& xi, const std::vector<double>& disp) const
{
  if (not weakkirchhoff_)
    dserror("method GetTriadAtXi is limited to WK so far! extend to SK if needed");

  if (disp.size() != 15)
    dserror("size mismatch: expected 15 values for element disp vector and got %d", disp.size());


  // Dof vector in total Lagrangian style, i.e. "displacement + reference values"
  CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double> disp_totlag(true);

  UpdateDispTotlag<2, double>(disp, disp_totlag);

  // material triads at collocation points
  std::vector<CORE::LINALG::Matrix<3, 3, double>> triad_mat_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<3, 3, double>(true));
  CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double> dummy(true);
  std::vector<CORE::LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<4, 1>(true));

  // Todo @grill:
  //    this method uses Qrefconv_[node] as reference triads; so we can only call this method if
  //    those were calculated correctly in a preceding time step. (be careful in post-processing!)
  UpdateNodalVariables<2, double>(disp_totlag, dummy, triad_mat_cp,
      Qref_dummy);  // Todo @grill split/adapt method and avoid dummy variables

  // create object of triad interpolation scheme
  Teuchos::RCP<
      LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS, double>>
      triad_interpolation_scheme_ptr = Teuchos::rcp(
          new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              double>());

  // reset scheme with nodal triads
  triad_interpolation_scheme_ptr->Reset(triad_mat_cp);

  triad_interpolation_scheme_ptr->GetInterpolatedTriadAtXi(triad, xi);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetScaledSecondAndThirdBaseVectorAtXi(const double& xi,
    const std::vector<double>& disp, CORE::LINALG::Matrix<3, 2>& scaledbasevectors) const
{
  // Todo @grill: delete or update this method. kept it for now as we might need it soon
  dserror(
      "Beam3k: method GetScaledSecondAndThirdBaseVectorAtXi is deprecated for now because it "
      "was only valid for the implicit assumption of a rectangular cross-section! If need be, "
      "generalize the definition of cross-section shape and dimensions and adapt this method. "
      "For now, the cross-section shape and dimensions are only required explicitly if beam "
      "interactions (contact, potentials, drag in background fluid) are to be evaluated. In this "
      "case, we only support and hence assume a circular cross-section with a radius which is "
      "either "
      "explicitly specified in the material definition line of the input file or per default "
      "computed "
      "from the area moment of inertia Iyy.");

  //  CORE::LINALG::Matrix<3,3> triad(true);
  //
  //  GetTriadAtXi(triad,xi,disp);
  //
  //  // ToDo careful, this is a hack for rectangular cross-sections ?!?
  //  double Ly=std::pow(pow(12.0*Izz_,3)/(12.0*Iyy_),1.0/8.0);
  //  double Lz=12.0*Izz_/pow(Ly,3);
  //
  //  for(int i=0;i<3;i++)
  //  {
  //    scaledbasevectors(i,0) = 0.5 * Ly * triad(i,1);
  //    scaledbasevectors(i,1) = 0.5 * Lz * triad(i,2);
  //  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetGeneralizedInterpolationMatrixVariationsAtXi(
    CORE::LINALG::SerialDenseMatrix& Ivar, const double& xi, const std::vector<double>& disp) const
{
  if (not weakkirchhoff_) dserror("method is limited to WK so far! extend to SK if needed");

  if (rotvec_)
    dserror("method is limited to tangent-based parametrization! extend to rotvec if needed");

  const unsigned int ndim = 3;
  const unsigned int nnodecl = 2;
  const unsigned int vpernode = 2;
  const unsigned int numdof = ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS;

  // safety check
  if ((unsigned int)Ivar.numRows() != 2 * ndim or (unsigned int) Ivar.numCols() != numdof)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", 6, numdof, Ivar.numRows(),
        Ivar.numCols());

  Ivar.putScalar(0.0);

  // *******************************************************************************
  // concerning interpolation of variation of CENTERLINE POSITION \vardelta r:
  // *******************************************************************************

  CORE::LINALG::Matrix<ndim, numdof, double> N(true);
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> N_i(true);


  CORE::DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);
  AssembleShapefunctionsN(N_i, N);

  // this part is associated with the variation of the centerline position
  // (first three rows of Ivar)
  for (unsigned int irow = 0; irow < N.numRows(); ++irow)
    for (unsigned int icol = 0; icol < N.numCols(); ++icol) Ivar(irow, icol) += N(irow, icol);

  // *******************************************************************************
  // concerning interpolation of variation of CENTERLINE ORIENTATION \vardelta \theta:
  // *******************************************************************************

  // define and initialize variables
  // *******************************************************************************
  // position index where CP quantities have to be stored (according to numbering convention)
  unsigned int ind = 0;
  double xi_cp = 0.0;


  CORE::LINALG::Matrix<ndim, numdof, double> N_s(true);
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> N_i_xi(true);

  CORE::LINALG::Matrix<1, numdof, double> L(true);
  CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, double> L_i(true);


  CORE::LINALG::Matrix<numdof, 1, double> disp_totlag(true);
  CORE::LINALG::Matrix<numdof, 1, double> disp_totlag_centerline(true);

  std::vector<CORE::LINALG::Matrix<3, 3, double>> triad_dummy(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<3, 3, double>(true));
  std::vector<CORE::LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<4, 1>(true));


  CORE::LINALG::Matrix<3, 1, double> r_s(true);  // r' vector
  double abs_r_s = 0.0;                          // ||r'||


  std::vector<CORE::LINALG::Matrix<numdof, ndim, double>> v_thetaperp_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<numdof, ndim, double>(true));
  std::vector<CORE::LINALG::Matrix<numdof, ndim, double>> v_thetapar_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<numdof, ndim, double>(true));


  // re-interpolated spin vector variation: v_theta_bar
  CORE::LINALG::Matrix<numdof, ndim, double> v_theta_bar(true);


  // set nodal / cp quantities: positions, tangents, triads
  // *******************************************************************************

  // Set current positions and orientations at all nodes:
  UpdateDispTotlag<nnodecl, double>(disp, disp_totlag);


  UpdateNodalVariables<nnodecl, double>(disp_totlag, disp_totlag_centerline, triad_dummy,
      Qref_dummy);  // Todo @grill: make this nicer and avoid dummies !


  // compute quantities at collocation points
  // *******************************************************************************
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    // Determine storage position for this cp
    ind = CORE::LARGEROTATIONS::NumberingTrafo(icp + 1, BEAM3K_COLLOCATION_POINTS);

    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi_cp = (double)icp / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;

    // get value of interpolating function for theta (Lagrange polynomials) at xi_cp
    L_i.Clear();
    CORE::DRT::UTILS::shape_function_1D(L_i, xi_cp, Shape());

    L.Clear();
    AssembleShapefunctionsL(L_i, L);


    N_i_xi.Clear();
    CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi_cp, length_, line2);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);


    // Calculation of r' at xi
    r_s.Clear();
    r_s.Multiply(N_s, disp_totlag_centerline);

    abs_r_s = CORE::FADUTILS::Norm(r_s);


    Calc_v_thetaperp<nnodecl>(v_thetaperp_cp[ind], N_s, r_s, abs_r_s);

    Calc_v_thetapartheta<nnodecl>(v_thetapar_cp[ind], L, r_s, abs_r_s);
  }


  // re-interpolation of quantities at xi based on CP values
  // *******************************************************************************

  L_i.Clear();
  CORE::DRT::UTILS::shape_function_1D(L_i, xi, Shape());

  v_theta_bar.Clear();
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    v_theta_bar.Update(L_i(icp), v_thetaperp_cp[icp], 1.0);
    v_theta_bar.Update(L_i(icp), v_thetapar_cp[icp], 1.0);
  }

  // finally assemble the generalized interpolation matrix for the variations Ivar
  // *******************************************************************************

  // this part is associated with the increment of the centerline orientation
  // (expressed as rotation vector theta)
  // (rows 4-6 of Iinc)
  // note: we need the transposed of v_theta_bar (rows <-> cols)
  for (unsigned int irow = 0; irow < v_theta_bar.numCols(); ++irow)
    for (unsigned int icol = 0; icol < v_theta_bar.numRows(); ++icol)
      Ivar(ndim + irow, icol) += v_theta_bar(icol, irow);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
    CORE::LINALG::SerialDenseMatrix& stiffmat, const double& xi, const std::vector<double>& disp,
    const CORE::LINALG::SerialDenseVector& force) const
{
  if (not weakkirchhoff_) dserror("method is limited to WK so far! extend to SK if needed");

  if (rotvec_)
    dserror("method is limited to tangent-based parametrization! extend to rotvec if needed");


  const unsigned int ndim = 3;
  const unsigned int nnodecl = 2;
  const unsigned int vpernode = 2;
  const unsigned int numdof = ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS;

  // safety check
  if ((unsigned int)stiffmat.numRows() != numdof or (unsigned int) stiffmat.numCols() != numdof)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", numdof, numdof,
        stiffmat.numRows(), stiffmat.numCols());

  stiffmat.putScalar(0.0);

  // create an auxiliary fixed size matrix and set as a view to original data in stiffmat
  CORE::LINALG::Matrix<numdof, numdof, double> stiffmat_fixedsize(stiffmat, true);


  CORE::LINALG::Matrix<ndim, 1, double> moment(&(force(3)));

  CORE::LINALG::Matrix<ndim, ndim, double> S_of_moment(true);
  CORE::LARGEROTATIONS::computespin<double>(S_of_moment, moment);


  // define and initialize variables
  // *******************************************************************************
  // position index where CP quantities have to be stored (according to numbering convention)
  unsigned int ind = 0;
  double xi_cp = 0.0;


  CORE::LINALG::Matrix<ndim, numdof, double> N_s(true);
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> N_i_xi(true);


  CORE::LINALG::Matrix<1, numdof, double> L(true);
  CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, double> L_i(true);


  CORE::LINALG::Matrix<numdof, 1, double> disp_totlag(true);
  CORE::LINALG::Matrix<numdof, 1, double> disp_totlag_centerline(true);

  std::vector<CORE::LINALG::Matrix<3, 3, double>> triad_mat_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<3, 3, double>(true));
  std::vector<CORE::LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<4, 1>(true));


  CORE::LINALG::Matrix<ndim, 1, double> r_s(true);  // r' vector
  double abs_r_s = 0.0;                             // ||r'||

  // first base vector at CP
  CORE::LINALG::Matrix<ndim, 1, double> g_1_cp(true);


  std::vector<CORE::LINALG::Matrix<numdof, numdof, double>> lin_v_thetaperp_moment_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<numdof, numdof, double>(true));

  std::vector<CORE::LINALG::Matrix<numdof, numdof, double>> lin_v_thetapar_moment_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<numdof, numdof, double>(true));


  // set nodal / cp quantities: positions, tangents, triads
  // *******************************************************************************

  // Set current positions and orientations at all nodes:
  UpdateDispTotlag<nnodecl, double>(disp, disp_totlag);


  UpdateNodalVariables<nnodecl, double>(disp_totlag, disp_totlag_centerline, triad_mat_cp,
      Qref_dummy);  // Todo make this nicer and avoid dummies!


  // compute quantities at collocation points
  // *******************************************************************************
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    // Determine storage position for this cp
    ind = CORE::LARGEROTATIONS::NumberingTrafo(icp + 1, BEAM3K_COLLOCATION_POINTS);

    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi_cp = (double)icp / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;

    // get value of interpolating function for theta (Lagrange polynomials) at xi_cp
    L_i.Clear();
    CORE::DRT::UTILS::shape_function_1D(L_i, xi_cp, Shape());

    L.Clear();
    AssembleShapefunctionsL(L_i, L);


    N_i_xi.Clear();
    CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi_cp, length_, line2);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);


    // Calculation of first base vector at xi
    r_s.Clear();
    r_s.Multiply(N_s, disp_totlag_centerline);

    abs_r_s = CORE::FADUTILS::Norm(r_s);

    g_1_cp.Clear();
    g_1_cp.Update(std::pow(abs_r_s, -1.0), r_s);


    Calc_lin_v_thetaperp_moment<nnodecl>(
        lin_v_thetaperp_moment_cp[ind], N_s, g_1_cp, abs_r_s, S_of_moment);

    Calc_lin_v_thetapar_moment<nnodecl>(
        lin_v_thetapar_moment_cp[ind], L, N_s, g_1_cp, abs_r_s, moment);
  }


  // re-interpolation of quantities at xi based on CP values
  // *******************************************************************************
  L_i.Clear();
  CORE::DRT::UTILS::shape_function_1D(L_i, xi, Shape());

  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    stiffmat_fixedsize.Update(L_i(icp), lin_v_thetaperp_moment_cp[ind], 1.0);
    stiffmat_fixedsize.Update(L_i(icp), lin_v_thetapar_moment_cp[ind], 1.0);
  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetGeneralizedInterpolationMatrixIncrementsAtXi(
    CORE::LINALG::SerialDenseMatrix& Iinc, const double& xi, const std::vector<double>& disp) const
{
  if (not weakkirchhoff_) dserror("method is limited to WK so far! extend to SK if needed");

  if (rotvec_)
    dserror("method is limited to tangent-based parametrization! extend to rotvec if needed");


  const unsigned int ndim = 3;
  const unsigned int nnodecl = 2;
  const unsigned int vpernode = 2;
  const unsigned int numdof = ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS;

  // safety check
  if ((unsigned int)Iinc.numRows() != 2 * ndim or (unsigned int) Iinc.numCols() != numdof)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", 6, numdof, Iinc.numRows(),
        Iinc.numCols());

  Iinc.putScalar(0.0);


  // *******************************************************************************
  // concerning interpolation of increment of CENTERLINE POSITION \Delta r:
  // *******************************************************************************

  CORE::LINALG::Matrix<ndim, numdof, double> N(true);
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> N_i(true);


  CORE::DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);
  AssembleShapefunctionsN(N_i, N);

  // this part is associated with the increment of the centerline position
  // (first three rows of Iinc)
  for (unsigned int irow = 0; irow < N.numRows(); ++irow)
    for (unsigned int icol = 0; icol < N.numCols(); ++icol) Iinc(irow, icol) += N(irow, icol);


  // *******************************************************************************
  // concerning interpolation of increment of CENTERLINE ORIENTATION \Delta \theta:
  // *******************************************************************************

  // define and initialize variables
  // *******************************************************************************
  // position index where CP quantities have to be stored (according to numbering convention)
  unsigned int ind = 0;
  double xi_cp = 0.0;


  CORE::LINALG::Matrix<ndim, numdof, double> N_s(true);
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> N_i_xi(true);

  CORE::LINALG::Matrix<1, numdof, double> L(true);
  CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, double> L_i(true);


  CORE::LINALG::Matrix<numdof, 1, double> disp_totlag(true);
  CORE::LINALG::Matrix<numdof, 1, double> disp_totlag_centerline(true);

  std::vector<CORE::LINALG::Matrix<3, 3, double>> triad_mat_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<3, 3, double>(true));
  std::vector<CORE::LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<4, 1>(true));


  CORE::LINALG::Matrix<3, 1, double> r_s(true);  // r' vector
  double abs_r_s = 0.0;                          // ||r'||

  // first base vector at CP
  CORE::LINALG::Matrix<3, 1, double> g_1_cp(true);
  // first base vector of the triad, from which the new intermediate triad is obtained via
  // smallest rotation (SR); this triad is arbitrary, but we choose the intermediate triad
  // of the last time step; see Dissertation Meier, p.25
  CORE::LINALG::Matrix<3, 1, double> g_1_cp_bar(true);


  CORE::LINALG::Matrix<ndim, numdof, double> lin_theta_perp_cp(true), lin_theta_par_cp(true);

  // lin_theta_cp = lin_theta_perp_cp + lin_theta_par_cp
  std::vector<CORE::LINALG::Matrix<ndim, numdof, double>> lin_theta_cp(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<ndim, numdof, double>(true));

  // re-interpolated lin_theta:
  CORE::LINALG::Matrix<ndim, numdof, double> lin_theta_bar(true);


  // set nodal / cp quantities: positions, tangents, triads
  // *******************************************************************************

  // Set current positions and orientations at all nodes:
  UpdateDispTotlag<nnodecl, double>(disp, disp_totlag);


  UpdateNodalVariables<nnodecl, double>(disp_totlag, disp_totlag_centerline, triad_mat_cp,
      Qref_dummy);  // Todo @grill: make this nicer and avoid dummies!


  // compute quantities at collocation points
  // *******************************************************************************
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    // Determine storage position for this cp
    ind = CORE::LARGEROTATIONS::NumberingTrafo(icp + 1, BEAM3K_COLLOCATION_POINTS);

    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi_cp = (double)icp / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;

    // get value of interpolating function for theta (Lagrange polynomials) at xi_cp
    L_i.Clear();
    CORE::DRT::UTILS::shape_function_1D(L_i, xi_cp, Shape());

    L.Clear();
    AssembleShapefunctionsL(L_i, L);

    N_i_xi.Clear();
    CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi_cp, length_, line2);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);

    // Calculation of r' at xi
    r_s.Clear();
    r_s.Multiply(N_s, disp_totlag_centerline);

    abs_r_s = CORE::FADUTILS::Norm(r_s);


    // lin_thetaperp
    Calc_lin_thetaperp<nnodecl>(lin_theta_perp_cp, N_s, r_s, abs_r_s);


    // v_lin_thetapar
    g_1_cp.Clear();
    g_1_cp.Update(std::pow(abs_r_s, -1.0), r_s);

    CORE::LINALG::Matrix<3, 3, double> triad_ref_conv_cp(true);
    CORE::LARGEROTATIONS::quaterniontotriad(Qrefconv_[ind], triad_ref_conv_cp);

    g_1_cp_bar.Clear();
    for (unsigned int idim = 0; idim < ndim; ++idim) g_1_cp_bar(idim) = triad_ref_conv_cp(idim, 0);

    Calc_lin_thetapar<nnodecl>(lin_theta_par_cp, L, N_s, g_1_cp, g_1_cp_bar, abs_r_s);

    // lin_theta
    lin_theta_cp[ind].Clear();
    lin_theta_cp[ind].Update(1.0, lin_theta_par_cp, 1.0, lin_theta_perp_cp);
  }


  // re-interpolation of quantities at xi based on CP values
  // *******************************************************************************
  std::vector<CORE::LINALG::Matrix<3, 3, double>> Itilde(
      BEAM3K_COLLOCATION_POINTS, CORE::LINALG::Matrix<3, 3, double>(true));

  // create object of triad interpolation scheme
  Teuchos::RCP<
      LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS, double>>
      triad_interpolation_scheme_ptr = Teuchos::rcp(
          new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              double>());

  // reset scheme with nodal triads
  triad_interpolation_scheme_ptr->Reset(triad_mat_cp);

  // compute Itilde matrices required for re-interpolation of CP values of lin_theta
  triad_interpolation_scheme_ptr->GetNodalGeneralizedRotationInterpolationMatricesAtXi(Itilde, xi);


  CORE::LINALG::Matrix<3, numdof, double> auxmatrix(true);

  lin_theta_bar.Clear();
  for (unsigned int inode = 0; inode < BEAM3K_COLLOCATION_POINTS; ++inode)
  {
    auxmatrix.Clear();

    auxmatrix.Multiply(Itilde[inode], lin_theta_cp[inode]);

    lin_theta_bar.Update(1.0, auxmatrix, 1.0);
  }


  // finally assemble the generalized interpolation matrix for the increments Iinc
  // *******************************************************************************

  // this part is associated with the increment of the centerline orientation
  // (expressed as rotation vector theta)
  // (rows 4-6 of Iinc)
  for (unsigned int irow = 0; irow < lin_theta_bar.numRows(); ++irow)
    for (unsigned int icol = 0; icol < lin_theta_bar.numCols(); ++icol)
      Iinc(ndim + irow, icol) += lin_theta_bar(irow, icol);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector(
    const std::vector<double>& dofvec, std::vector<double>& dofvec_centerline,
    bool add_reference_values) const
{
  if (dofvec.size() != 15)
    dserror("size mismatch: expected 15 values for element state vector and got %d", dofvec.size());

  dofvec_centerline.resize(12, 0.0);

  // we use the method for CORE::LINALG fixed size matrix and create it as a view on the STL vector
  CORE::LINALG::Matrix<15, 1, double> dofvec_fixedsize(dofvec.data());
  CORE::LINALG::Matrix<12, 1, double> dofvec_centerline_fixedsize(dofvec_centerline.data(), true);

  this->ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
      dofvec_fixedsize, dofvec_centerline_fixedsize, add_reference_values);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector(
    const CORE::LINALG::Matrix<3 * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& dofvec,
    CORE::LINALG::Matrix<3 * vpernode * nnodecl, 1, T>& dofvec_centerline,
    bool add_reference_values) const
{
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  const int dofperboundarynode = 3 * vpernode + 1;

  // get current values for position DOFs directly
  for (unsigned int dim = 0; dim < 3; ++dim)
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      dofvec_centerline(3 * vpernode * node + dim) = dofvec(dofperboundarynode * node + dim);
    }

  // Hermite interpolation: get values of tangent DOFs as well
  if (rotvec_ == false)
  {
    for (unsigned int dim = 0; dim < 3; ++dim)
      for (unsigned int node = 0; node < nnodecl; ++node)
      {
        // in case of rotvec_==false, tangents are equivalent to nodal DOFs
        dofvec_centerline(3 * vpernode * node + 3 + dim) =
            dofvec(dofperboundarynode * node + 3 + dim);
      }
  }
  else
  {
    // values for tangent DOFs must be transformed in case of rotvec_==true
    CORE::LINALG::Matrix<3, 1, T> theta(true);
    CORE::LINALG::Matrix<3, 3, T> triad(true);

    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      for (unsigned int dim = 0; dim < 3; ++dim)
      {
        // get values for rotation vector DOFs
        theta(dim) = dofvec(dofperboundarynode * node + 3 + dim);
      }
      // transform to triad
      CORE::LARGEROTATIONS::angletotriad(theta, triad);

      // direction of tangent is equivalent to first base vector of triad; length of tangent is 7th
      // nodal DOF
      for (unsigned int dim = 0; dim < 3; ++dim)
      {
        dofvec_centerline(3 * vpernode * node + 3 + dim) =
            triad(dim, 0) * dofvec(dofperboundarynode * node + 6);
      }
    }
  }

  if (add_reference_values) AddRefValuesDispCenterline<nnodecl, vpernode, T>(dofvec_centerline);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::AddRefValuesDisp(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& dofvec) const
{
  // For the boundary nodes we have to add the initial values. For the interior nodes we are already
  // done since the initial relative angles are zero: alpha_i(t=0)=0;
  if (rotvec_ == false)
  {
    // Calculate total displacements = positions vectors, tangents and relative angles
    for (unsigned int node = 0; node < 2; ++node)  // loop over boundary nodes
    {
      for (unsigned int ndof = 0; ndof < 7; ++ndof)  // loop over dofs per node
      {
        if (ndof < 3)
        {
          dofvec(7 * node + ndof) += (Nodes()[node])->X()[ndof];
        }
        else if (ndof < 6)
        {
          dofvec(7 * node + ndof) += (Tref()[node])(ndof - 3);
        }
        else
        {
          // nothing to do here: alpha_i(t=0)=0;
        }
      }
    }
  }
  else
  {
    // Calculate total displacements = positions vectors, tangents and relative angles
    for (unsigned int node = 0; node < 2; ++node)  // loop over boundary nodes
    {
      for (unsigned int ndof = 0; ndof < 7; ++ndof)  // loop over dofs per node
      {
        if (ndof < 3)
        {
          dofvec(7 * node + ndof) += (Nodes()[node])->X()[ndof];
        }
        else if (ndof < 6)
        {
          // Nothing to do here, rotations are treated below
        }
        else
        {
          // here we have to add the initial length of the tangents at the boundary nodes,
          // i.e. ||r'_i(t=0)||=1:
          dofvec(7 * node + ndof) += 1.0;
        }
      }  // for (int ndof=0;ndof<7;ndof++)//loop over dofs per node

      CORE::LINALG::Matrix<3, 1> disptheta(true);
      CORE::LINALG::Matrix<3, 1> thetanew(true);
      CORE::LINALG::Matrix<4, 1> deltaQ(true);
      CORE::LINALG::Matrix<4, 1> Qnew(true);
      for (unsigned int i = 0; i < 3; ++i)
      {
        disptheta(i) = CORE::FADUTILS::CastToDouble(dofvec(7 * node + 3 + i));
      }

      CORE::LARGEROTATIONS::angletoquaternion(disptheta, deltaQ);

      CORE::LINALG::Matrix<4, 1> Q0;
      CORE::LARGEROTATIONS::angletoquaternion(Theta0()[node], Q0);
      CORE::LARGEROTATIONS::quaternionproduct(Q0, deltaQ, Qnew);

      // renormalize quaternion to keep its absolute value one even in case of long simulations
      // and intricate calculations
      Qnew.Scale(1.0 / Qnew.Norm2());

      // Calculate the new nodal angle thetanew \in ]-PI,PI] -> Here, thetanew \in ]-PI,PI] by
      // quaterniontoangle()
      CORE::LARGEROTATIONS::quaterniontoangle(Qnew, thetanew);

      // Finally set rotation values in disp_totlag
      for (unsigned int i = 0; i < 3; ++i)
      {
        dofvec(7 * node + 3 + i) = thetanew(i);
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::UpdateDispTotlag(const std::vector<double>& disp,
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag) const
{
  disp_totlag.Clear();
  for (unsigned int dof = 0; dof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++dof)
    disp_totlag(dof) = disp[dof];

  AddRefValuesDisp<nnodecl, T>(disp_totlag);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::UpdateNodalVariables(
    const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline,
    std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
    std::vector<CORE::LINALG::Matrix<4, 1>>& Qref_new) const
{
  // Set positions vectors and tangents and triads at boundary nodes
  SetPositionsAtBoundaryNodes<nnodecl, T>(disp_totlag, disp_totlag_centerline);

  // next, set triads and tangents at boundary nodes
  SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes<nnodecl, T>(
      disp_totlag, disp_totlag_centerline, triad_mat_cp, Qref_new);

  // finally, set triads at remaining CPs (all except boundary nodes)
  SetTriadsAndReferenceTriadsAtRemainingCollocationPoints<nnodecl, T>(
      disp_totlag, disp_totlag_centerline, triad_mat_cp, Qref_new);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::SetPositionsAtBoundaryNodes(
    const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline)
    const
{
  for (unsigned int i = 0; i < 3; ++i)
  {
    disp_totlag_centerline(i) = disp_totlag(i);
    disp_totlag_centerline(7 + i) = disp_totlag(7 + i);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes(
    const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline,
    std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
    std::vector<CORE::LINALG::Matrix<4, 1>>& Qref_new) const
{
  if (rotvec_ == false)
  {
    CORE::LINALG::Matrix<3, 1, T> tangent(true);
    CORE::LINALG::Matrix<3, 3, T> triad_ref(true);
    CORE::LINALG::Matrix<3, 3> triad_aux(true);
    // Todo @grill: get rid of auxiliary matrix of different type
    CORE::LINALG::Matrix<3, 3, T> triad_aux2(true);
    T alpha = 0.0;

    for (unsigned int node = 0; node < 2; ++node)
    {
      // tangent
      for (unsigned int i = 0; i < 3; ++i)
      {
        disp_totlag_centerline(7 * node + 3 + i) = disp_totlag(7 * node + 3 + i);
        tangent(i) = disp_totlag(7 * node + 3 + i);
      }

      alpha = disp_totlag(7 * node + 6);

      // calculate new sr triad
      triad_ref.Clear();
      triad_aux.Clear();
      triad_aux2.Clear();
      CORE::LARGEROTATIONS::quaterniontotriad(Qrefconv_[node], triad_aux);
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j) triad_aux2(i, j) = triad_aux(i, j);

      CORE::LARGEROTATIONS::CalculateSRTriads<T>(tangent, triad_aux2, triad_ref);

      // Store nodal reference triad
      CORE::LINALG::Matrix<4, 1, T> Qref(true);
      CORE::LARGEROTATIONS::triadtoquaternion(triad_ref, Qref);
      for (unsigned int i = 0; i < 4; ++i)
        Qref_new[node](i) = CORE::FADUTILS::CastToDouble(Qref(i));

      // calculate material triad
      triad_mat_cp[node].Clear();
      CORE::LARGEROTATIONS::RotateTriad(triad_ref, alpha, triad_mat_cp[node]);
    }
  }
  else
  {
    CORE::LINALG::Matrix<3, 1, T> theta(true);
    for (unsigned int node = 0; node < 2; ++node)
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        theta(i) = disp_totlag(7 * node + 3 + i);
      }
      triad_mat_cp[node].Clear();
      CORE::LARGEROTATIONS::angletotriad(theta, triad_mat_cp[node]);

      // tangent
      for (unsigned int i = 0; i < 3; ++i)
      {
        disp_totlag_centerline(7 * node + 3 + i) =
            (triad_mat_cp[node])(i, 0) * disp_totlag(7 * node + 6);
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::SetTriadsAndReferenceTriadsAtRemainingCollocationPoints(
    const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
        disp_totlag_centerline,
    std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
    std::vector<CORE::LINALG::Matrix<4, 1>>& Qref_new) const
{
  CORE::LINALG::Matrix<3, 1, T> tangent(true);
  CORE::LINALG::Matrix<3, 3, T> triad_ref(true);
  CORE::LINALG::Matrix<3, 3> triad_aux(true);
  // Todo @grill: get rid of auxiliary matrix of different type
  CORE::LINALG::Matrix<3, 3, T> triad_aux2(true);
  T alpha = 0.0;

  CORE::LINALG::Matrix<1, 4, T> N_i_xi(true);
  CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, T> N_s(true);
  CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, T> L_i;
  double xi = 0.0;
  unsigned int ind = 0;

  //********begin: evaluate quantities at collocation points********************************
  // intermediate nodes: all but first (node=0) and last node (node=BEAM3K_COLLOCATION_POINTS)
  for (unsigned int node = 1; node < BEAM3K_COLLOCATION_POINTS - 1; ++node)
  {
    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0 node=2->xi=1
    xi = (double)node / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;
    N_i_xi.Clear();
    CORE::DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
    L_i.Clear();
    CORE::DRT::UTILS::shape_function_1D(L_i, xi, Shape());

    // Determine storage position for the node node
    ind = CORE::LARGEROTATIONS::NumberingTrafo(node + 1, BEAM3K_COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);

    tangent.Clear();
    // Calculation of r' at xi
    tangent.Multiply(N_s, disp_totlag_centerline);

    alpha = disp_totlag(7 * 2 + ind - 2);

    // calculate new sr triads
    triad_ref.Clear();
    triad_aux.Clear();
    triad_aux2.Clear();
    CORE::LARGEROTATIONS::quaterniontotriad(Qrefconv_[ind], triad_aux);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j) triad_aux2(i, j) = triad_aux(i, j);

    CORE::LARGEROTATIONS::CalculateSRTriads<T>(tangent, triad_aux2, triad_ref);

    // Store nodal reference triad
    CORE::LINALG::Matrix<4, 1, T> Qref(true);
    CORE::LARGEROTATIONS::triadtoquaternion(triad_ref, Qref);
    for (unsigned int i = 0; i < 4; ++i) Qref_new[ind](i) = CORE::FADUTILS::CastToDouble(Qref(i));

    // calculate material triad
    triad_mat_cp[ind].Clear();
    CORE::LARGEROTATIONS::RotateTriad(triad_ref, alpha, triad_mat_cp[ind]);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::SetAutomaticDifferentiationVariables(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>& disp_totlag) const
{
  for (unsigned int dof = 0; dof < nnodecl * 6 + BEAM3K_COLLOCATION_POINTS; dof++)
  {
    disp_totlag(dof).diff(dof, nnodecl * 6 + BEAM3K_COLLOCATION_POINTS);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3k::Calc_velocity(
    const CORE::LINALG::Matrix<ndim * vpernode * nnodecl, 1, double>& velocity_dofvec,
    const CORE::LINALG::Matrix<1, vpernode * nnodecl, double>& N_i,
    CORE::LINALG::Matrix<ndim, 1, double>& velocity,
    const CORE::LINALG::Matrix<ndim, 1, double>& position, int gausspoint_index) const
{
  DRT::UTILS::BEAM::CalcInterpolation<nnodecl, vpernode, ndim, double>(
      velocity_dofvec, N_i, velocity);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3k::Calc_velocity(
    const CORE::LINALG::Matrix<ndim * vpernode * nnodecl, 1, double>& velocity_dofvec,
    const CORE::LINALG::Matrix<1, vpernode * nnodecl, double>& N_i,
    CORE::LINALG::Matrix<ndim, 1, FAD>& velocity,
    const CORE::LINALG::Matrix<ndim, 1, FAD>& position, int gausspoint_index)
{
  /* if we use FAD, we need to track the dependency of velocity vector on primary variables, i.e.
   * we must calculate it like this using the class variable rconvmass_ for position vector at GP
   * of last converged state this again is tricky in case of periodic boundary conditions, because
   * the position of this element (nodes) might have been shifted outside; we therefore manually
   * adapt the calculated velocity and compare it with the velocity calculated in time integrator
   * and handed in from outside for safety reasons */
  CORE::LINALG::Matrix<ndim, 1, double> velocity_test;
  DRT::UTILS::BEAM::CalcInterpolation<nnodecl, vpernode, ndim, double>(
      velocity_dofvec, N_i, velocity_test);

  // get time step size
  const double dt = ParamsInterface().GetDeltaTime();

  CORE::LINALG::Matrix<3, 1> diff(true);

  CORE::LINALG::Matrix<ndim, 1, FAD> delta_r_ost(true);
  Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> pbb =
      BrownianDynParamsInterface().GetPeriodicBoundingBox();

  CORE::LINALG::Matrix<3, 1> unshiftedrconvmass_i(true), position_i_double(true);

  for (unsigned int idim = 0; idim < ndim; ++idim)
  {
    unshiftedrconvmass_i(idim) = rconvmass_[gausspoint_index](idim);
    position_i_double(idim) = CORE::FADUTILS::CastToDouble(position(idim));
  }

  // difference in position of this GP as compared to last time step
  pbb->UnShift3D(unshiftedrconvmass_i, position_i_double);

  for (unsigned int idim = 0; idim < ndim; ++idim)
  {
    // difference in position of this GP as compared to last time step
    delta_r_ost(idim) = position(idim) - unshiftedrconvmass_i(idim);

    // velocity according to Backward Euler scheme
    velocity(idim) = delta_r_ost(idim) / dt;

    diff(idim) = velocity(idim).val() - velocity_test(idim);

    // set class variable, such that rconvmass_ is available in next time step
    rnewmass_[gausspoint_index](idim) = position(idim).val();
  }

  // safety check
  if (diff.NormInf() > 1e-11)
  {
    std::cout << "\nrnewmass = " << position;
    std::cout << "\nrconvmass_ = " << rconvmass_[gausspoint_index];
    std::cout << "\nvel = " << velocity;
    std::cout << "\nvel_test = " << velocity_test;
    std::cout << "\nabs(diff) = " << diff.Norm2();
    std::cout << "\n*** SERIOUS WARNING: ***";
    std::cout << "\nvelocity vector at GP computed locally in beam3k differs from OneStepTheta "
                 "velocity vector (see values above)!\n\n";
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::Calc_v_thetaperp(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>& v_thetaperp,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>& N_s,
    const CORE::LINALG::Matrix<3, 1, T>& r_s, T abs_r_s) const
{
  v_thetaperp.Clear();

  CORE::LINALG::Matrix<3, 3, T> S_of_r_s(true);
  CORE::LARGEROTATIONS::computespin<T>(S_of_r_s, r_s);

  v_thetaperp.MultiplyTN(N_s, S_of_r_s);
  v_thetaperp.Scale(-1.0 * std::pow(abs_r_s, -2.0));
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::Calc_v_thetapartheta(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>& v_thetapartheta,
    const CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>& L,
    const CORE::LINALG::Matrix<3, 1, T>& r_s, T abs_r_s) const
{
  v_thetapartheta.Clear();

  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      v_thetapartheta(idof, jdim) = L(idof) * r_s(jdim) / abs_r_s;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_thetaperp(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_thetaperp,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& r_s, double abs_r_s) const
{
  // Todo @grill: maybe re-use method for v_thetaperp here, which is simply the transpose of this
  // term

  lin_thetaperp.Clear();

  CORE::LINALG::Matrix<3, 3, double> S_of_r_s(true);
  CORE::LARGEROTATIONS::computespin<double>(S_of_r_s, r_s);

  lin_thetaperp.Multiply(S_of_r_s, N_s);
  lin_thetaperp.Scale(std::pow(abs_r_s, -2.0));
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_thetapar(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_thetapar,
    const CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1,
    const CORE::LINALG::Matrix<3, 1, double>& g_1_bar, double abs_r_s) const
{
  const unsigned int numdof = 6 * nnodecl + BEAM3K_COLLOCATION_POINTS;

  lin_thetapar.Clear();

  CORE::LINALG::Matrix<3, 3, double> S_of_g_1(true);
  CORE::LARGEROTATIONS::computespin<double>(S_of_g_1, g_1);


  // Todo @grill: decide about alternatives, apparently no change in results

  // *********** alternative 1 ************
  CORE::LINALG::Matrix<3, 3, double> g_1_dyadicproduct_g_1_bar_T(true);

  for (unsigned int irow = 0; irow < 3; ++irow)
    for (unsigned int icol = 0; icol < 3; ++icol)
      g_1_dyadicproduct_g_1_bar_T(irow, icol) = g_1(irow) * g_1_bar(icol);

  CORE::LINALG::Matrix<3, 3, double> g_1_dyadicproduct_g_1_bar_T_S_of_g_1(true);
  g_1_dyadicproduct_g_1_bar_T_S_of_g_1.Multiply(g_1_dyadicproduct_g_1_bar_T, S_of_g_1);

  lin_thetapar.Multiply(g_1_dyadicproduct_g_1_bar_T_S_of_g_1, N_s);

  // *********** alternative 2 ************

  //  CORE::LINALG::Matrix<1,3,double> g_1_bar_T_S_of_g_1(true);
  //  g_1_bar_T_S_of_g_1.MultiplyTN(g_1_bar, S_of_g_1);
  //
  //
  //  CORE::LINALG::Matrix<1,numdof,double> g_1_bar_T_S_of_g_1_N_s(true);
  //  g_1_bar_T_S_of_g_1_N_s.Multiply(g_1_bar_T_S_of_g_1, N_s);
  //
  //
  //  // dyadic product
  //  for (unsigned int idim=0; idim<3; ++idim)
  //    for (unsigned int jdof=0; jdof<numdof; ++jdof)
  //      lin_thetapar(idim,jdof) += g_1(idim) * g_1_bar_T_S_of_g_1_N_s(jdof) ;

  // *************************************

  lin_thetapar.Scale(-1.0 / (1.0 + g_1.Dot(g_1_bar)));
  lin_thetapar.Scale(1.0 / abs_r_s);


  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdof = 0; jdof < numdof; ++jdof)
      lin_thetapar(idim, jdof) += g_1(idim) * L(jdof);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_tangent_tilde,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const
{
  lin_tangent_tilde.Clear();

  CORE::LINALG::Matrix<3, 3, double> auxmatrix;

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - 2.0 * g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -2.0));

  lin_tangent_tilde.Multiply(auxmatrix, N_s);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde_s(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_tangent_tilde_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
    const CORE::LINALG::Matrix<3, 1, double>& r_s, const CORE::LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s) const
{
  lin_tangent_tilde_s.Clear();

  CORE::LINALG::Matrix<3, 3, double> auxmatrix;
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> summand(true);

  // first summand
  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - 2.0 * g_1(idim) * g_1(jdim);

  auxmatrix.Scale(2.0 * r_s.Dot(r_ss) * std::pow(abs_r_s, -4.0));

  lin_tangent_tilde_s.Multiply(auxmatrix, N_s);

  // second summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = g_1_s(idim) * g_1(jdim) + g_1(idim) * g_1_s(jdim);

  auxmatrix.Scale(-2.0 * std::pow(abs_r_s, -2.0));

  summand.Multiply(auxmatrix, N_s);

  lin_tangent_tilde_s.Update(1.0, summand, 1.0);

  // third summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - 2.0 * g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -2.0));

  summand.Clear();

  summand.Multiply(auxmatrix, N_ss);

  lin_tangent_tilde_s.Update(1.0, summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_g_1(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_g_1,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const
{
  lin_g_1.Clear();

  CORE::LINALG::Matrix<3, 3, double> auxmatrix;

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -1.0));

  lin_g_1.Multiply(auxmatrix, N_s);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_g_1_s(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_g_1_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
    const CORE::LINALG::Matrix<3, 1, double>& r_s, const CORE::LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s) const
{
  lin_g_1_s.Clear();

  CORE::LINALG::Matrix<3, 3, double> auxmatrix;
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> summand(true);

  // first summand
  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - g_1(idim) * g_1(jdim);

  auxmatrix.Scale(-1.0 * r_s.Dot(r_ss) * std::pow(abs_r_s, -3.0));

  lin_g_1_s.Multiply(auxmatrix, N_s);

  // second summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = g_1_s(idim) * g_1(jdim) + g_1(idim) * g_1_s(jdim);

  auxmatrix.Scale(-1.0 * std::pow(abs_r_s, -1.0));

  summand.Multiply(auxmatrix, N_s);

  lin_g_1_s.Update(1.0, summand, 1.0);

  // third summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -1.0));

  summand.Clear();

  summand.Multiply(auxmatrix, N_ss);

  lin_g_1_s.Update(1.0, summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_epsilon(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_epsilon,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const
{
  lin_v_epsilon.Clear();

  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1(true);

  Calc_lin_g_1<nnodecl>(lin_g_1, N_s, g_1, abs_r_s);

  lin_v_epsilon.MultiplyTN(N_s, lin_g_1);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_moment_resultant(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_moment_resultant,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta_s,
    const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment,
    const CORE::LINALG::Matrix<3, 3, double>& cm) const
{
  lin_moment_resultant.Clear();

  // first summand
  lin_moment_resultant.Multiply(spinmatrix_of_moment, lin_theta);
  lin_moment_resultant.Scale(-1.0);

  // second summand
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> auxmatrix(true);

  auxmatrix.Multiply(cm, lin_theta_s);

  lin_moment_resultant.Update(1.0, auxmatrix, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_moment_inertia(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_moment_inertia,
    const CORE::LINALG::Matrix<3, 3, double>& triad_mat,
    const CORE::LINALG::Matrix<3, 3, double>& triad_mat_conv,
    const CORE::LINALG::Matrix<3, 1, double>& deltatheta,
    const CORE::LINALG::Matrix<3, 1, double>& angular_velocity_material,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
    const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment,
    const CORE::LINALG::Matrix<3, 3, double>& C_rho, double lin_prefactor_acc,
    double lin_prefactor_vel) const
{
  lin_moment_inertia.Clear();

  // first summand
  lin_moment_inertia.Multiply(spinmatrix_of_moment, lin_theta);

  // second summand
  CORE::LINALG::Matrix<3, 3, double> auxmatrix(true);

  CORE::LINALG::Matrix<3, 3, double> spinmatrix(true);
  CORE::LARGEROTATIONS::computespin(spinmatrix, angular_velocity_material);

  auxmatrix.Multiply(spinmatrix, C_rho);

  CORE::LINALG::Matrix<3, 1, double> C_rho_W(true);
  C_rho_W.Multiply(C_rho, angular_velocity_material);

  spinmatrix.Clear();
  CORE::LARGEROTATIONS::computespin(spinmatrix, C_rho_W);
  auxmatrix.Update(-1.0, spinmatrix, 1.0);
  auxmatrix.Scale(lin_prefactor_vel);

  auxmatrix.Update(lin_prefactor_acc, C_rho, 1.0);


  CORE::LINALG::Matrix<3, 3, double> Tmat_of_deltatheta = CORE::LARGEROTATIONS::Tmatrix(deltatheta);
  CORE::LINALG::Matrix<3, 3, double> Lambda_conv_Tmat_of_deltatheta(true);
  Lambda_conv_Tmat_of_deltatheta.MultiplyTN(triad_mat_conv, Tmat_of_deltatheta);

  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_THETA_tilde(true);
  lin_THETA_tilde.Multiply(Lambda_conv_Tmat_of_deltatheta, lin_theta);

  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> auxmatrix2(true);
  auxmatrix2.Multiply(auxmatrix, lin_THETA_tilde);


  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> auxmatrix3(true);
  auxmatrix3.Multiply(triad_mat, auxmatrix2);

  lin_moment_inertia.Update(1.0, auxmatrix3, 1.0);

  lin_moment_inertia.Scale(-1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_moment_viscous(
    CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_moment_viscous,
    const CORE::LINALG::Matrix<3, 3, double>& triad_mat,
    const CORE::LINALG::Matrix<3, 3, double>& triad_mat_conv,
    const CORE::LINALG::Matrix<3, 1, double>& deltatheta,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
    const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment, double gamma_polar,
    double dt) const
{
  lin_moment_viscous.Clear();

  // first summand
  lin_moment_viscous.Multiply(spinmatrix_of_moment, lin_theta);
  lin_moment_viscous.Scale(-1.0);

  // second summand
  CORE::LINALG::Matrix<3, 3, double> auxmatrix(true);

  CORE::LINALG::Matrix<3, 3, double> gamma_g1_g1_conv;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      gamma_g1_g1_conv(i, j) = triad_mat(i, 0) * triad_mat_conv(j, 0) * gamma_polar / dt;

  CORE::LINALG::Matrix<3, 3, double> Tmat_of_deltatheta = CORE::LARGEROTATIONS::Tmatrix(deltatheta);

  auxmatrix.Multiply(gamma_g1_g1_conv, Tmat_of_deltatheta);

  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> second_summand(true);
  second_summand.Multiply(auxmatrix, lin_theta);

  lin_moment_viscous.Update(1.0, second_summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_moment(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetaperp_moment,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
    const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment) const
{
  lin_v_thetaperp_moment.Clear();

  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_tangent_tilde(true);

  Calc_lin_tangent_tilde<nnodecl>(lin_tangent_tilde, N_s, g_1, abs_r_s);


  CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  // Todo @grill: check: is the order of matrix products relevant?
  auxmatrix.MultiplyTN(N_s, spinmatrix_of_moment);

  lin_v_thetaperp_moment.Multiply(auxmatrix, lin_tangent_tilde);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_s_moment(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetaperp_s_moment,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
    const CORE::LINALG::Matrix<3, 1, double>& r_s, const CORE::LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s, const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment) const
{
  lin_v_thetaperp_s_moment.Clear();

  // first summand
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_tangent_tilde(true);

  Calc_lin_tangent_tilde<nnodecl>(lin_tangent_tilde, N_s, g_1, abs_r_s);

  CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  auxmatrix.MultiplyTN(N_ss, spinmatrix_of_moment);

  lin_v_thetaperp_s_moment.Multiply(auxmatrix, lin_tangent_tilde);

  // second summand
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_tangent_tilde_s(
      true);

  Calc_lin_tangent_tilde_s<nnodecl>(lin_tangent_tilde_s, N_s, N_ss, g_1, g_1_s, r_s, r_ss, abs_r_s);

  auxmatrix.Clear();
  auxmatrix.MultiplyTN(N_s, spinmatrix_of_moment);

  CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
      6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>
      second_summand(true);

  second_summand.Multiply(auxmatrix, lin_tangent_tilde_s);

  lin_v_thetaperp_s_moment.Update(1.0, second_summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_moment(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetapar_moment,
    CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
    const CORE::LINALG::Matrix<3, 1, double>& moment) const
{
  lin_v_thetapar_moment.Clear();

  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1(true);

  Calc_lin_g_1<nnodecl>(lin_g_1, N_s, g_1, abs_r_s);

  CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  // Todo @grill: check: is the order of matrix products relevant?
  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim) auxmatrix(idof, jdim) = L(idof) * moment(jdim);

  lin_v_thetapar_moment.Multiply(auxmatrix, lin_g_1);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_s_moment(
    CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetapar_s_moment,
    CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
    CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
    const CORE::LINALG::Matrix<3, 1, double>& r_s, const CORE::LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s, const CORE::LINALG::Matrix<3, 1, double>& moment) const
{
  lin_v_thetapar_s_moment.Clear();

  // first summand
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1(true);

  Calc_lin_g_1<nnodecl>(lin_g_1, N_s, g_1, abs_r_s);

  CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim) auxmatrix(idof, jdim) = L_s(idof) * moment(jdim);

  lin_v_thetapar_s_moment.Multiply(auxmatrix, lin_g_1);

  // second summand
  CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1_s(true);

  Calc_lin_g_1_s<nnodecl>(lin_g_1_s, N_s, N_ss, g_1, g_1_s, r_s, r_ss, abs_r_s);

  auxmatrix.Clear();

  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim) auxmatrix(idof, jdim) = L(idof) * moment(jdim);

  CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
      6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>
      second_summand(true);

  second_summand.Multiply(auxmatrix, lin_g_1_s);

  lin_v_thetapar_s_moment.Update(1.0, second_summand, 1.0);
}

// explicit template instantiations
template void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector<2, 2,
    Sacado::Fad::DFad<double>>(
    const CORE::LINALG::Matrix<12 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    CORE::LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&, bool) const;
template void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
    const CORE::LINALG::Matrix<12 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    CORE::LINALG::Matrix<12, 1, double>&, bool) const;

template void DRT::ELEMENTS::Beam3k::AddRefValuesDisp<2, double>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&) const;
template void DRT::ELEMENTS::Beam3k::AddRefValuesDisp<2, Sacado::Fad::DFad<double>>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3k::UpdateDispTotlag<2, double>(const std::vector<double>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&) const;
template void DRT::ELEMENTS::Beam3k::UpdateDispTotlag<2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3k::UpdateNodalVariables<2, double>(
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1>>&) const;
template void DRT::ELEMENTS::Beam3k::UpdateNodalVariables<2, Sacado::Fad::DFad<double>>(
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>>&,
    std::vector<CORE::LINALG::Matrix<4, 1>>&) const;

template void DRT::ELEMENTS::Beam3k::SetPositionsAtBoundaryNodes<2, double>(
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&) const;
template void DRT::ELEMENTS::Beam3k::SetPositionsAtBoundaryNodes<2, Sacado::Fad::DFad<double>>(
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3k::SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes<2,
    double>(const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1>>&) const;
template void DRT::ELEMENTS::Beam3k::SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes<2,
    Sacado::Fad::DFad<double>>(
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>>&,
    std::vector<CORE::LINALG::Matrix<4, 1>>&) const;

template void DRT::ELEMENTS::Beam3k::SetTriadsAndReferenceTriadsAtRemainingCollocationPoints<2,
    double>(const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1>>&) const;
template void DRT::ELEMENTS::Beam3k::SetTriadsAndReferenceTriadsAtRemainingCollocationPoints<2,
    Sacado::Fad::DFad<double>>(
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    const CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>>&,
    std::vector<CORE::LINALG::Matrix<4, 1>>&) const;

template void DRT::ELEMENTS::Beam3k::SetAutomaticDifferentiationVariables<2>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, FAD>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_velocity<2, 2, 3>(
    const CORE::LINALG::Matrix<3 * 2 * 2, 1, double>&,
    const CORE::LINALG::Matrix<1, 2 * 2, double>&, CORE::LINALG::Matrix<3, 1, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, int) const;

template void DRT::ELEMENTS::Beam3k::Calc_velocity<2, 2, 3>(
    const CORE::LINALG::Matrix<3 * 2 * 2, 1, double>&,
    const CORE::LINALG::Matrix<1, 2 * 2, double>&, CORE::LINALG::Matrix<3, 1, FAD>&,
    const CORE::LINALG::Matrix<3, 1, FAD>&, int);

template void DRT::ELEMENTS::Beam3k::Calc_v_thetaperp<2, double>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double) const;
template void DRT::ELEMENTS::Beam3k::Calc_v_thetaperp<2, Sacado::Fad::DFad<double>>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, Sacado::Fad::DFad<double>>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, Sacado::Fad::DFad<double>>&,
    const CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>>&, Sacado::Fad::DFad<double>) const;

template void DRT::ELEMENTS::Beam3k::Calc_v_thetapartheta<2, double>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, double>&,
    const CORE::LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double) const;
template void DRT::ELEMENTS::Beam3k::Calc_v_thetapartheta<2, Sacado::Fad::DFad<double>>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, Sacado::Fad::DFad<double>>&,
    const CORE::LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, Sacado::Fad::DFad<double>>&,
    const CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>>&, Sacado::Fad::DFad<double>) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_thetaperp<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_thetapar<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde_s<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_g_1<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_g_1_s<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_epsilon<2>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS,
        double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_moment_resultant<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 3, double>&, const CORE::LINALG::Matrix<3, 3, double>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_moment_inertia<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 3, double>&, const CORE::LINALG::Matrix<3, 3, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 3, double>&, const CORE::LINALG::Matrix<3, 3, double>&, double,
    double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_moment_viscous<2>(
    CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 3, double>&, const CORE::LINALG::Matrix<3, 3, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 3, double>&, double, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_moment<2>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS,
        double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, double,
    const CORE::LINALG::Matrix<3, 3, double>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_moment<2>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_thetapar_moment,
    CORE::LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& L,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
    const CORE::LINALG::Matrix<3, 1, double>& moment) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_s_moment<2>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS,
        double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&,
    const CORE::LINALG::Matrix<3, 1, double>&, const CORE::LINALG::Matrix<3, 1, double>&, double,
    const CORE::LINALG::Matrix<3, 3, double>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_s_moment<2>(
    CORE::LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_thetapar_s_moment,
    CORE::LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& L,
    CORE::LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& L_s,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const CORE::LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const CORE::LINALG::Matrix<3, 1, double>& g_1, const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
    const CORE::LINALG::Matrix<3, 1, double>& r_s, const CORE::LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s, const CORE::LINALG::Matrix<3, 1, double>& moment) const;