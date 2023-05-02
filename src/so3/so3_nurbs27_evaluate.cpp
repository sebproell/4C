/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation methods of the quadratic NURBS 27 element

\level 2

*----------------------------------------------------------------------*/
#include "so3_nurbs27.H"
#include "lib_discret.H"
#include "nurbs_discret.H"
#include "lib_utils.H"
#include "lib_dserror.H"
#include "linalg_utils_sparse_algebra_math.H"
#include "linalg_serialdensevector.H"
#include <Epetra_SerialDenseSolver.h>
#include "discretization_fem_general_utils_integration.H"
#include "discretization_fem_general_utils_fem_shapefunctions.H"
#include "discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "mat_so3_material.H"
#include "lib_globalproblem.H"
#include "lib_elements_paramsinterface.H"
#include "so3_utils.H"

#include "lib_element_vtk_cell_type_register.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "nurbs_discret_nurbs_utils.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  LINALG::Matrix<81, 81> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<81, 81> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<81, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<81, 1> elevec2(elevec2_epetra.A(), true);

  // start with "none"
  DRT::ELEMENTS::NURBS::So_nurbs27::ActionType act = So_nurbs27::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = So_nurbs27::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = So_nurbs27::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_nurbs27::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = So_nurbs27::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_nurbs27::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_eleload")
    act = So_nurbs27::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = So_nurbs27::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = So_nurbs27::calc_struct_update_istep;
  else if (action == "calc_stc_matrix")
    act = So_nurbs27::calc_stc_matrix;
  else if (action == "calc_stc_matrix_inverse")
    act = So_nurbs27::calc_stc_matrix_inverse;
  else if (action == "calc_struct_reset_istep")
    act = So_nurbs27::calc_struct_reset_istep;
  else if (action == "calc_struct_energy")
    act = So_nurbs27::calc_struct_energy;
  else if (action == "calc_struct_nlnstifflmass")
    act = So_nurbs27::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action '%s' for So_nurbs27", action.c_str());
  // what should the element do
  switch (act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (double& i : mydisp) i = 0.0;
      std::vector<double> myres(lm.size());
      for (double& myre : myres) myre = 0.0;
      sonurbs27_nlnstiffmass(
          lm, discretization, mydisp, myres, &elemat1, nullptr, &elevec1, params);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      LINALG::Matrix<81, 81>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      sonurbs27_nlnstiffmass(lm, discretization, mydisp, myres, matptr, nullptr, &elevec1, params);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<81, 81> myemat(true);
      sonurbs27_nlnstiffmass(lm, discretization, mydisp, myres, &myemat, nullptr, &elevec1, params);
    }
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      sonurbs27_nlnstiffmass(
          lm, discretization, mydisp, myres, &elemat1, &elemat2, &elevec1, params);

      if (act == calc_struct_nlnstifflmass) lumpmass(&elemat2);
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    case calc_struct_update_istep:
    {
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    case calc_stc_matrix_inverse:
    {
      const auto stc_scaling = DRT::INPUT::get<INPAR::STR::STC_Scale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        dserror("To scale or not to scale, that's the query!");
      else
      {
        CalcSTCMatrix(elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, true);
      }
    }
    break;

    case calc_stc_matrix:
    {
      const auto stc_scaling = DRT::INPUT::get<INPAR::STR::STC_Scale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        dserror("To scale or not to scale, that's the query!");
      else
      {
        CalcSTCMatrix(
            elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, false);
      }
    }
    break;

    case calc_struct_energy:
    {
      if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");

      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      elevec1_epetra(0) = CalcIntEnergy(discretization, mydisp, params);
      break;
    }

    default:
      dserror("Unknown type of action for So_nurbs27");
  }
  return 0;
}  // DRT::ELEMENTS::So_nurbs27::Evaluate


/*----------------------------------------------------------------------*
 | calc. scaled thickness matrix for thin shell-like structs   (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::CalcSTCMatrix(LINALG::Matrix<81, 81>& elemat1,
    const INPAR::STR::STC_Scale stc_scaling, const int stc_layer, std::vector<int>& lm,
    DRT::Discretization& discretization, bool do_inverse)
{
  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if (nurbsdis == nullptr)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele)
  {
    return;
  }

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    auto* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }


  // --------------------------------------------------
  // determine the lengths in r-, s- and t-direction

  // compute coordinates  of corners 0,2,6,18

  LINALG::Matrix<27, 1> funct;

  LINALG::Matrix<3, 1> x0;
  LINALG::Matrix<3, 1> x2;
  LINALG::Matrix<3, 1> x6;
  LINALG::Matrix<3, 1> x18;

  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = -1.0;
    gpa(1) = -1.0;
    gpa(2) = -1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x0(isd) = val;
    }
  }

  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = 1.0;
    gpa(1) = -1.0;
    gpa(2) = -1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x2(isd) = val;
    }
  }
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = 1.0;
    gpa(1) = 1.0;
    gpa(2) = -1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x6(isd) = val;
    }
  }
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = -1.0;
    gpa(1) = -1.0;
    gpa(2) = 1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x18(isd) = val;
    }
  }

  LINALG::Matrix<3, 1> deltaX;

  deltaX.Update(1.0, x2, -1.0, x0);
  const double length_r = deltaX.Norm2();
  deltaX.Update(1.0, x6, -1.0, x0);
  const double length_s = deltaX.Norm2();
  deltaX.Update(1.0, x18, -1.0, x0);
  const double length_t = deltaX.Norm2();

  double ratio = 1.0;

  std::vector<int> topnodeids;
  std::vector<int> midnodeids;
  std::vector<int> botnodeids;

  if (length_t <= length_r && length_t <= length_s)
  {
    for (int i = 0; i < 9; ++i) botnodeids.push_back(i);
    for (int i = 9; i < 18; ++i) midnodeids.push_back(i);
    for (int i = 18; i < 27; ++i) topnodeids.push_back(i);

    ratio = (length_r + length_s) / (2.0 * length_t);
  }
  else if (length_s <= length_r && length_s <= length_t)
  {
    botnodeids.push_back(0);
    botnodeids.push_back(1);
    botnodeids.push_back(2);
    botnodeids.push_back(9);
    botnodeids.push_back(10);
    botnodeids.push_back(11);
    botnodeids.push_back(18);
    botnodeids.push_back(19);
    botnodeids.push_back(20);

    midnodeids.push_back(3);
    midnodeids.push_back(4);
    midnodeids.push_back(5);
    midnodeids.push_back(12);
    midnodeids.push_back(13);
    midnodeids.push_back(14);
    midnodeids.push_back(21);
    midnodeids.push_back(22);
    midnodeids.push_back(23);

    topnodeids.push_back(6);
    topnodeids.push_back(7);
    topnodeids.push_back(8);
    topnodeids.push_back(15);
    topnodeids.push_back(16);
    topnodeids.push_back(17);
    topnodeids.push_back(24);
    topnodeids.push_back(25);
    topnodeids.push_back(26);

    ratio = (length_r + length_t) / (2.0 * length_s);
  }
  else if (length_r <= length_s && length_r <= length_t)
  {
    for (int i = 0; i < 27; i += 3) botnodeids.push_back(i);

    for (int i = 1; i < 27; i += 3) midnodeids.push_back(i);

    for (int i = 2; i < 27; i += 3) topnodeids.push_back(i);

    ratio = (length_t + length_s) / (2.0 * length_r);
  }


  double C = 1.0;
  if (stc_scaling == INPAR::STR::stc_currsym)
  {
    C = ratio;
  }
  else
  {
    C = ratio * ratio;
  }


  double fac1 = 0.0;
  double fac2 = 0.0;

  if (do_inverse)
  {
    fac1 = (1.0 - C);
    fac2 = C;
  }
  else
  {
    fac1 = (C - 1.0) / (C);
    fac2 = 1.0 / C;
  }

  LINALG::Matrix<27, 1> adjele(true);

  for (int i = 0; i < 27; i++)
  {
    adjele(i, 0) = nodes[i]->NumElement();
  }
  /*
    // loop row midnode
    for(int i=0; i<9; i++)
      {
        int dvi=3*midnodeids[i];
        int dui=3*topnodeids[i];
        int dwi=3*botnodeids[i];

        for(int j=0; j<3; j++)
        {
          elemat1(dvi+j,dvi+j)+=fac2/adjele(midnodeids[i],0);
          elemat1(dvi+j,dui+j)+=fac1/adjele(midnodeids[i],0);
          elemat1(dvi+j,dwi+j)+=fac1/adjele(midnodeids[i],0);
        }
      }

    // loop row botnode
    for(int i=0; i<9; i++)
      {
        int dvi=3*botnodeids[i];

        for(int j=0; j<3; j++)
          {
            elemat1(dvi+j,dvi+j)+=1.0/adjele(botnodeids[i],0);
          }
      }

    // loop row topnode
    for(int i=0; i<9; i++)
      {
        int dvi=3*topnodeids[i];

        for(int j=0; j<3; j++)
          {
            elemat1(dvi+j,dvi+j)+=1.0/adjele(topnodeids[i],0);
          }
      }

  */

  // loop row midnode
  for (int i = 0; i < 9; i++)
  {
    int dvi = 3 * midnodeids[i];

    for (int j = 0; j < 3; j++) elemat1(dvi + j, dvi + j) += 1.0 / adjele(midnodeids[i], 0);
  }

  // loop row botnode
  for (int i = 0; i < 9; i++)
  {
    int dvi = 3 * botnodeids[i];
    int dui = 3 * midnodeids[i];

    for (int j = 0; j < 3; j++)
    {
      elemat1(dvi + j, dvi + j) += fac2 * 1.0 / adjele(botnodeids[i], 0);
      elemat1(dvi + j, dui + j) += fac1 * 1.0 / adjele(botnodeids[i], 0);
    }
  }

  // loop row topnode
  for (int i = 0; i < 9; i++)
  {
    int dvi = 3 * topnodeids[i];
    int dui = 3 * midnodeids[i];

    for (int j = 0; j < 3; j++)
    {
      elemat1(dvi + j, dvi + j) += fac2 * 1.0 / adjele(topnodeids[i], 0);
      elemat1(dvi + j, dui + j) += fac1 * 1.0 / adjele(topnodeids[i], 0);
    }
  }

  return;
}  // CalcSTCMatrix



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);
  // get values and switches from the condition
  const auto* onoff = condition.Get<std::vector<int>>("onoff");
  const auto* val = condition.Get<std::vector<double>>("val");

  /*
   **    TIME CURVE BUSINESS
   */
  // find out whether we will use a time curve
  double time = -1.0;
  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SONURBS27)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SONURBS27; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const auto* funct = condition.Get<std::vector<int>>("funct");
  LINALG::Matrix<NUMDIM_SONURBS27, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SONURBS27; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if (nurbsdis == nullptr) dserror("So_nurbs27 appeared in non-nurbs discretisation\n");

  // there is nothing to be done for zero sized elements in knotspan
  if ((*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id())) return (0);

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
    weights(inode) = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode])->W();

  /*------------------------------------------------------------------*/
  /*                   update element geometry                        */
  /*------------------------------------------------------------------*/

  // material coord. of element
  LINALG::Matrix<27, 3> xrefe;
  for (int i = 0; i < 27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  const int numgp = 27;
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::GaussRule3D::hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  LINALG::Matrix<3, 1> gpa;

  LINALG::Matrix<27, 1> shape;
  LINALG::Matrix<3, 27> deriv;

  for (int gp = 0; gp < numgp; ++gp)
  {
    gpa(0) = intpoints.qxg[gp][0];
    gpa(1) = intpoints.qxg[gp][1];
    gpa(2) = intpoints.qxg[gp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        shape, deriv, gpa, myknots, weights, DRT::Element::nurbs27);

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SONURBS27, NUMDIM_SONURBS27> jac;
    jac.Multiply(deriv, xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SONURBS27; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SONURBS27; ++nodid)
          xrefegp(dim) += shape(nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    const double fac = intpoints.qwgt[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SONURBS27; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0) ? DRT::Problem::Instance()
                                 ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(functnum - 1)
                                 .Evaluate(xrefegp.A(), time, dim)
                           : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SONURBS27; ++nodid)
        {
          elevec1[nodid * NUMDIM_SONURBS27 + dim] += shape(nodid) * dim_fac;
        }
      }
    }

  } /* end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_nurbs27::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::InitJacobianMapping(DRT::Discretization& dis)
{
  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

  if (nurbsdis == nullptr)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele)
  {
    return;
  }

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    auto* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  const static std::vector<LINALG::Matrix<3, 27>> derivs = sonurbs27_derivs(myknots, weights);
  LINALG::Matrix<27, 3> xrefe;
  for (int i = 0; i < 27; ++i)
  {
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }

  const int numgp = 27;

  invJ_.resize(numgp);
  detJ_.resize(numgp);
  for (int gp = 0; gp < numgp; ++gp)
  {
    invJ_[gp].Multiply(derivs[gp], xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT %12.5e IN ELEMENT ID %d, gauss point %d", detJ_[gp],
          Id(), gp);
  }
  return;
}  // DRT::ELEMENTS::So_nurbs27::InitJacobianMapping()

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_nlnstiffmass(
    std::vector<int>& lm,                 // location matrix
    DRT::Discretization& discretization,  // discretisation to extract knot vector
    std::vector<double>& disp,            // current displacements
    std::vector<double>& residual,        // current residual displ
    LINALG::Matrix<81, 81>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<81, 81>* massmatrix,   // element mass matrix
    LINALG::Matrix<81, 1>* force,         // element internal force vector
    Teuchos::ParameterList& params)       // strain output option
{
  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if (nurbsdis == nullptr)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele)
  {
    return;
  }

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    auto* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // update element geometry
  LINALG::Matrix<27, 3> xrefe;  // material coord. of element
  LINALG::Matrix<27, 3> xcurr;  // current  coord. of element
  for (int i = 0; i < 27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * 3];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * 3 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * 3 + 2];
  }

  /*------------------------------------------------------------------*/
  /*                    Loop over Gauss Points                        */
  /*------------------------------------------------------------------*/
  const int numgp = 27;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::GaussRule3D::hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  invJ_.resize(numgp);
  detJ_.resize(numgp);

  LINALG::Matrix<27, 1> funct;
  LINALG::Matrix<3, 27> deriv;

  LINALG::Matrix<3, 27> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<3, 3> defgrd(true);
  for (int gp = 0; gp < numgp; ++gp)
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = intpoints.qxg[gp][0];
    gpa(1) = intpoints.qxg[gp][1];
    gpa(2) = intpoints.qxg[gp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        funct, deriv, gpa, myknots, weights, DRT::Element::nurbs27);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<3, 3> invJac(true);

    invJac.Multiply(deriv, xrefe);
    double detJ = invJac.Invert();

    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT %12.5e IN ELEMENT ID %d, gauss point %d", detJ_[gp],
          Id(), gp);

    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJac, deriv);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<3, 3> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(6);
    LINALG::Matrix<6, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<6, 81> bop;
    for (int i = 0; i < 27; ++i)
    {
      bop(0, 3 * i) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, 3 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, 3 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, 3 * i) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, 3 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, 3 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, 3 * i) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, 3 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, 3 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, 3 * i) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, 3 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, 3 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, 3 * i) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, 3 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, 3 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, 3 * i) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, 3 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, 3 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // call material law
    LINALG::Matrix<6, 6> cmat(true);
    LINALG::Matrix<6, 1> stress(true);
    UTILS::GetTemperatureForStructuralMaterial<nurbs27>(funct, params);
    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
    // end of call material law

    double detJ_w = detJ * intpoints.qwgt[gp];
    // update internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // update stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6, 81> cb;
      cb.Multiply(cmat, bop);
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
      sfac.Scale(detJ_w);                 // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);       // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < 27; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < 27; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < 3; ++idim)
          {
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          }

          (*stiffmatrix)(3 * inod, 3 * jnod) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness
    }    // if (stiffmatrix)

    if (massmatrix != nullptr)  // evaluate mass matrix
    {
      double density = Material()->Density(gp);
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < 27; ++inod)
      {
        ifactor = funct(inod) * factor;
        for (int jnod = 0; jnod < 27; ++jnod)
        {
          massfactor = funct(jnod) * ifactor;  // intermediate factor
          (*massmatrix)(3 * inod, 3 * jnod) += massfactor;
          (*massmatrix)(3 * inod + 1, 3 * jnod + 1) += massfactor;
          (*massmatrix)(3 * inod + 2, 3 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix

  } /* end of Loop over GP */

  return;
}  // DRT::ELEMENTS::So_nurbs27::sonurbs27_nlnstiffmass

/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Shape fcts at all 27 Gauss Points                     |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<27, 1>> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_shapefcts(
    const std::vector<Epetra_SerialDenseVector>& myknots, const LINALG::Matrix<27, 1>& weights)
{
  const int numgp = 27;

  std::vector<LINALG::Matrix<27, 1>> shapefcts(numgp);
  // (r,s,t) gp-locations of fully integrated quadratic Nurbs 27
  // fill up nodal f at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::GaussRule3D::hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp)
  {
    LINALG::Matrix<3, 1> gp;
    gp(0) = intpoints.qxg[igp][0];
    gp(1) = intpoints.qxg[igp][1];
    gp(2) = intpoints.qxg[igp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct(
        shapefcts[igp], gp, myknots, weights, DRT::Element::nurbs27);
  }
  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Shape fct derivs at all 27 Gauss Points              |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<3, 27>> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_derivs(
    const std::vector<Epetra_SerialDenseVector>& myknots, const LINALG::Matrix<27, 1>& weights)
{
  const int numgp = 27;

  std::vector<LINALG::Matrix<3, 27>> derivs(numgp);
  // (r,s,t) gp-locations of fully integrated quadratic Nurbs 27
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::GaussRule3D::hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp)
  {
    LINALG::Matrix<3, 1> gp;
    gp(0) = intpoints.qxg[igp][0];
    gp(1) = intpoints.qxg[igp][1];
    gp(2) = intpoints.qxg[igp][2];

    LINALG::Matrix<27, 1> dummyfct;

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        dummyfct, derivs[igp], gp, myknots, weights, DRT::Element::nurbs27);
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Weights at all 27 Gauss Points                     |
 *----------------------------------------------------------------------*/
const std::vector<double> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_gpweights()
{
  const int numgp = 27;

  std::vector<double> gpweights(numgp);
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::GaussRule3D::hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int i = 0; i < numgp; ++i)
  {
    gpweights[i] = intpoints.qwgt[i];
  }
  return gpweights;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::NURBS::So_nurbs27*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_nurbs27* failed");
    actele->InitJacobianMapping(dis);
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate internal energy of the element (private)                  |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::NURBS::So_nurbs27::CalcIntEnergy(
    DRT::Discretization& discretization,  // discretisation to extract knot vector
    std::vector<double>& disp,            // current displacements
    Teuchos::ParameterList& params)       // strain output option
{
  double energy = 0.;

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));
  if (nurbsdis == nullptr) dserror("So_nurbs27 appeared in non-nurbs discretisation\n");

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele) return 0.;

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    auto* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // update element geometry
  LINALG::Matrix<27, 3> xrefe;  // material coord. of element
  LINALG::Matrix<27, 3> xcurr;  // current  coord. of element
  for (int i = 0; i < 27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * 3];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * 3 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * 3 + 2];
  }
  /*------------------------------------------------------------------*/
  /*                    Loop over Gauss Points                        */
  /*------------------------------------------------------------------*/
  const int numgp = 27;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::GaussRule3D::hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  invJ_.resize(numgp);
  detJ_.resize(numgp);

  LINALG::Matrix<27, 1> funct;
  LINALG::Matrix<3, 27> deriv;

  LINALG::Matrix<3, 27> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<3, 3> defgrd(true);
  for (int gp = 0; gp < numgp; ++gp)
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = intpoints.qxg[gp][0];
    gpa(1) = intpoints.qxg[gp][1];
    gpa(2) = intpoints.qxg[gp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        funct, deriv, gpa, myknots, weights, DRT::Element::nurbs27);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<3, 3> invJac(true);

    invJac.Multiply(deriv, xrefe);
    double detJ = invJac.Invert();

    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT %12.5e IN ELEMENT ID %d, gauss point %d", detJ_[gp],
          Id(), gp);

    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJac, deriv);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<3, 3> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(6);
    LINALG::Matrix<6, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    double psi = 0.0;
    SolidMaterial()->StrainEnergy(glstrain, psi, gp, Id());

    double detJ_w = detJ * intpoints.qwgt[gp];
    energy += detJ_w * psi;
  }

  return energy;
}

/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::lumpmass(
    LINALG::Matrix<NUMDOF_SONURBS27, NUMDOF_SONURBS27>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/**
 * \brief Helper function to evaluate the NURBS interpolation inside the element.
 */
template <unsigned int n_points, unsigned int n_val>
void EvalNurbs3DInterpolation(LINALG::Matrix<n_val, 1, double>& r,
    const LINALG::Matrix<n_points * n_val, 1, double>& q, const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<n_points, 1, double>& weights,
    const std::vector<Epetra_SerialDenseVector>& myknots,
    const DRT::Element::DiscretizationType& distype)
{
  // Get the shape functions.
  LINALG::Matrix<n_points, 1, double> N;
  DRT::NURBS::UTILS::nurbs_get_3D_funct(N, xi, myknots, weights, distype);

  // Multiply the shape functions with the control point values.
  r.Clear();
  for (unsigned int i_node_nurbs = 0; i_node_nurbs < n_points; i_node_nurbs++)
  {
    for (unsigned int i_dim = 0; i_dim < n_val; i_dim++)
    {
      r(i_dim) += N(i_node_nurbs) * q(i_node_nurbs * 3 + i_dim);
    }
  }
}

/**
 *
 */
unsigned int DRT::ELEMENTS::NURBS::So_nurbs27::AppendVisualizationGeometry(
    const DRT::Discretization& discret, std::vector<uint8_t>& cell_types,
    std::vector<double>& point_coordinates) const
{
  // This NURBS element will be displayed like a hex27 element in the vtk output.
  const int number_of_ouput_points = 27;
  const auto vtk_cell_info =
      DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(DRT::Element::hex27);
  const std::vector<int>& numbering = vtk_cell_info.second;

  // Add the cell type to the output.
  cell_types.push_back(vtk_cell_info.first);

  // Create the "nodes" for the hex27 visualization.
  {
    // Get the knots and weights for this element.
    LINALG::Matrix<27, 1, double> weights(true);
    std::vector<Epetra_SerialDenseVector> myknots(true);
    const bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discret, this, myknots, weights);
    if (zero_size) dserror("GetMyNurbsKnotsAndWeights has to return a non zero size.");

    // Get the control points position in the reference configuration.
    LINALG::Matrix<27 * 3, 1, double> q;
    for (unsigned int i_node = 0; i_node < (unsigned int)this->NumNode(); ++i_node)
    {
      const DRT::Node* node = this->Nodes()[i_node];
      for (int i_dim = 0; i_dim < 3; ++i_dim) q(3 * i_node + i_dim) = node->X()[i_dim];
    }

    // Loop over the "nodes" of the hex27 element.
    LINALG::Matrix<3, 1, double> r;
    LINALG::Matrix<3, 1, double> xi;
    for (unsigned int i_node_hex = 0; i_node_hex < number_of_ouput_points; i_node_hex++)
    {
      for (unsigned int i = 0; i < 3; i++)
        xi(i) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[numbering[i_node_hex]][i];

      // Get the reference position at the parameter coordinate.
      EvalNurbs3DInterpolation(r, q, xi, weights, myknots, this->Shape());

      for (unsigned int i_dim = 0; i_dim < 3; i_dim++) point_coordinates.push_back(r(i_dim));
    }
  }

  return number_of_ouput_points;
}

/**
 *
 */
unsigned int DRT::ELEMENTS::NURBS::So_nurbs27::AppendVisualizationDofBasedResultDataVector(
    const DRT::Discretization& discret, const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
    unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
    std::vector<double>& vtu_point_result_data) const
{
  if (read_result_data_from_dofindex != 0)
    dserror("Nurbs output is only implemented for read_result_data_from_dofindex == 0");

  if (result_num_dofs_per_node != 3)
    dserror("The nurbs elements can only output nodal data with dimension 3, e.g., displacements");

  // This NURBS element will be displayed like a hex27 element in the vtk output.
  const int number_of_ouput_points = 27;
  const std::vector<int>& numbering =
      DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(DRT::Element::hex27).second;

  // Add the data at the "nodes" of the hex27 visualization.
  {
    // Get the knots and weights for this element.
    LINALG::Matrix<27, 1, double> weights(true);
    std::vector<Epetra_SerialDenseVector> myknots(true);
    const bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discret, this, myknots, weights);
    if (zero_size) dserror("GetMyNurbsKnotsAndWeights has to return a non zero size.");

    // Get the element result vector.
    LINALG::Matrix<27 * 3, 1, double> q;
    std::vector<double> eledisp;
    std::vector<int> lm, lmowner, lmstride;
    this->LocationVector(discret, lm, lmowner, lmstride);
    DRT::UTILS::ExtractMyValues(*result_data_dofbased, eledisp, lm);
    q.SetView(eledisp.data());

    // Loop over the "nodes" of the hex27 element.
    LINALG::Matrix<3, 1, double> r;
    LINALG::Matrix<3, 1, double> xi;
    for (unsigned int i_node_hex = 0; i_node_hex < number_of_ouput_points; i_node_hex++)
    {
      for (unsigned int i = 0; i < 3; i++)
        xi(i) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[numbering[i_node_hex]][i];

      // Get the reference position at the parameter coordinate.
      EvalNurbs3DInterpolation(r, q, xi, weights, myknots, this->Shape());

      for (unsigned int i_dim = 0; i_dim < 3; i_dim++) vtu_point_result_data.push_back(r(i_dim));
    }
  }

  return number_of_ouput_points;
}