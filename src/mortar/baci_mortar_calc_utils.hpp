/*-----------------------------------------------------------------------*/
/*! \file
\brief utils for mortar buisiness

\level 1

*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 01/14 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_CALC_UTILS_HPP
#define FOUR_C_MORTAR_CALC_UTILS_HPP

/*----------------------------------------------------------------------*
 | Header                                                   farah 01/14 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_mortar_element.hpp"
#include "baci_mortar_node.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | Utils                                                    farah 01/14 |
 *----------------------------------------------------------------------*/
namespace MORTAR
{
  class Element;
  class Node;

  namespace UTILS
  {
    /*----------------------------------------------------------------------*
     |  Get global coords for given local coords                 farah 01/14|
     *----------------------------------------------------------------------*/
    template <CORE::FE::CellType distype>
    bool LocalToGlobal(MORTAR::Element& ele, const double* xi, double* globcoord, int inttype)
    {
      // check input
      if (!xi) dserror("ERROR: LocalToGlobal called with xi=nullptr");
      if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=nullptr");

      static constexpr int n = CORE::FE::num_nodes<distype>;
      static constexpr int ndim = CORE::FE::dim<distype> + 1;

      DRT::Node** mynodes = ele.Points();
      if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

      std::fill(globcoord, globcoord + ndim, 0.0);

      //===========================================
      // Basic value
      switch (inttype)
      {
        case 0:
        {
          CORE::LINALG::Matrix<n, 1> val;
          switch (ndim)
          {
            case 2:
            {
              switch (distype)
              {
                case CORE::FE::CellType::nurbs2:
                case CORE::FE::CellType::nurbs3:
                {
                  CORE::LINALG::SerialDenseVector auxval(n);
                  CORE::LINALG::SerialDenseMatrix auxderiv(n, 1);
                  ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

                  for (int i = 0; i < n; ++i) val(i) = auxval(i);

                  break;
                }
                default:
                {
                  CORE::FE::shape_function_1D(val, xi[0], distype);
                  break;
                }
              }
              break;
            }
            case 3:
            {
              switch (distype)
              {
                case CORE::FE::CellType::nurbs4:
                case CORE::FE::CellType::nurbs8:
                case CORE::FE::CellType::nurbs9:
                {
                  CORE::LINALG::SerialDenseVector auxval(n);
                  CORE::LINALG::SerialDenseMatrix auxderiv(n, 2);
                  ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

                  for (int i = 0; i < n; ++i) val(i) = auxval(i);

                  break;
                }
                default:
                {
                  CORE::FE::shape_function_2D(val, xi[0], xi[1], distype);
                  break;
                }
              }

              break;
            }
            default:
            {
              dserror("Wrong Dimension");
              exit(EXIT_FAILURE);
            }
          }

          for (int i = 0; i < n; ++i)
          {
            Node* mymrtrnode = static_cast<Node*>(mynodes[i]);
            dsassert(mymrtrnode, "ERROR: LocalToGlobal: Null pointer!");

            for (int j = 0; j < ndim; ++j)
            {
              // use shape function values for interpolation
              globcoord[j] += val(i) * mymrtrnode->xspatial()[j];
            }
          }
          break;
        }
        //===========================================
        // First Derivation (xi)
        case 1:
        {
          CORE::LINALG::Matrix<ndim - 1, n> deriv1;

          switch (ndim)
          {
            case 2:
            {
              switch (distype)
              {
                case CORE::FE::CellType::nurbs2:
                case CORE::FE::CellType::nurbs3:
                {
                  CORE::LINALG::SerialDenseVector auxval(n);
                  CORE::LINALG::SerialDenseMatrix auxderiv(n, 1);
                  ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

                  for (int i = 0; i < n; ++i) deriv1(0, i) = auxderiv(i, 0);

                  break;
                }
                default:
                {
                  CORE::FE::shape_function_1D_deriv1(deriv1, xi[0], distype);
                  break;
                }
              }

              break;
            }
            case 3:
            {
              switch (distype)
              {
                case CORE::FE::CellType::nurbs4:
                case CORE::FE::CellType::nurbs8:
                case CORE::FE::CellType::nurbs9:
                {
                  CORE::LINALG::SerialDenseVector auxval(n);
                  CORE::LINALG::SerialDenseMatrix auxderiv(n, 2);
                  ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

                  for (int i = 0; i < n; ++i) deriv1(0, i) = auxderiv(i, 0);

                  break;
                }
                default:
                {
                  CORE::FE::shape_function_2D_deriv1(deriv1, xi[0], xi[1], distype);
                  break;
                }
              }

              break;
            }
            default:
            {
              dserror("Wrong Dimension");
              break;
            }
          }
          for (int i = 0; i < n; ++i)
          {
            Node* mymrtrnode = static_cast<Node*>(mynodes[i]);
            dsassert(mymrtrnode, "ERROR: LocalToGlobal: Null pointer!");

            for (int j = 0; j < ndim; ++j)
            {
              // use shape function values for interpolation
              globcoord[j] += deriv1(0, i) * mymrtrnode->xspatial()[j];
            }
          }
          break;
        }
        //===========================================
        // Second Derivation (eta)
        case 2:
        {
          CORE::LINALG::Matrix<ndim - 1, n> deriv2(true);
          switch (ndim)
          {
            case 2:
            {
              switch (distype)
              {
                case CORE::FE::CellType::nurbs2:
                case CORE::FE::CellType::nurbs3:
                {
                  CORE::LINALG::SerialDenseVector auxval(n);
                  CORE::LINALG::SerialDenseMatrix auxderiv(n, 1);
                  ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

                  for (int i = 0; i < n; ++i) deriv2(0, i) = auxderiv(i, 0);

                  break;
                }
                default:
                {
                  CORE::FE::shape_function_1D_deriv1(deriv2, xi[0], distype);
                  break;
                }
              }

              break;
            }
            case 3:
            {
              switch (distype)
              {
                case CORE::FE::CellType::nurbs4:
                case CORE::FE::CellType::nurbs8:
                case CORE::FE::CellType::nurbs9:
                {
                  CORE::LINALG::SerialDenseVector auxval(n);
                  CORE::LINALG::SerialDenseMatrix auxderiv(n, 2);
                  ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

                  for (int i = 0; i < n; ++i) deriv2(1, i) = auxderiv(i, 1);

                  break;
                }
                default:
                {
                  CORE::FE::shape_function_2D_deriv1(deriv2, xi[0], xi[1], distype);
                  break;
                }
              }

              break;
            }
            default:
            {
              dserror("Wrong Dimension");
              exit(EXIT_FAILURE);
            }
          }

          for (int i = 0; i < n; ++i)
          {
            Node* mymrtrnode = static_cast<Node*>(mynodes[i]);
            dsassert(mymrtrnode, "ERROR: LocalToGlobal: Null pointer!");

            if constexpr (ndim > 2)
            {
              for (int j = 0; j < ndim; ++j)
              {
                // use shape function values for interpolation
                globcoord[j] += deriv2(1, i) * mymrtrnode->xspatial()[j];
              }
            }
            else
              dserror("Not implemented.");
          }
          break;
        }
        default:
        {
          dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");
          exit(EXIT_FAILURE);
        }
      }
      return true;
    };


    /*----------------------------------------------------------------------*
     |  Get global coords for given local coords (ref pos)       farah 01/14|
     *----------------------------------------------------------------------*/
    template <CORE::FE::CellType distype>
    bool LocalToGlobal(DRT::Element& ele, const double* xi, double* globcoord)
    {
      // check input
      if (!xi) dserror("ERROR: LocalToGlobal called with xi=nullptr");
      if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=nullptr");

      static constexpr int n = CORE::FE::num_nodes<distype>;
      static constexpr int ndim = CORE::FE::dim<distype> + 1;

      DRT::Node** mynodes = ele.Nodes();
      if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

      for (int i = 0; i < ndim; ++i) globcoord[i] = 0.0;

      CORE::LINALG::Matrix<ndim, n> coord;

      CORE::LINALG::Matrix<n, 1> val;
      if (ndim == 2)
        CORE::FE::shape_function_1D(val, xi[0], distype);
      else if (ndim == 3)
        CORE::FE::shape_function_2D(val, xi[0], xi[1], distype);
      else
        dserror("Wrong Dimension");

      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < ndim; ++j)
        {
          coord(j, i) = mynodes[i]->X()[j];

          // use shape function values for interpolation
          globcoord[j] += val(i) * coord(j, i);
        }
      }


      return true;
    };
    /*---------------------------------------------------------------------------*
     |  Get global coords for given local coords (ref pos)            rauch 04/14|
     |  in current configuration.                                                |
     |                                                                           |
     |  template parameter distype is the discr. type of input element           |
     |                                                                           |
     |  ele     (in): element to consider.                                       |
     |  globdim (in): dimension of global problem                                |
     |  xi      (in): local reference coordinate within "ele".                   |
     |  edisp   (in): nodal displacements of "ele".                              |
     |                                                                           |
     |  globcoord (out): global xyz coordinate in current configuration of "xi"  |
     *---------------------------------------------------------------------------*/
    template <CORE::FE::CellType distype>
    bool LocalToCurrentGlobal(DRT::Element& ele, const int globdim, const double* xi,
        const std::vector<double> edisp, double* globcoord)
    {
      // check input
      if (!xi) dserror("ERROR: LocalToGlobal called with xi=nullptr");
      if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=nullptr");

      static constexpr int n = CORE::FE::num_nodes<distype>;
      static constexpr int ndim = CORE::FE::dim<distype>;

      if ((int)edisp.size() != n * globdim)
        dserror("ERROR: vector of element displacements has wrong dimension (%d != %d)",
            n * globdim, edisp.size());

      DRT::Node** mynodes = ele.Nodes();
      if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

      for (int i = 0; i < globdim; ++i) globcoord[i] = 0.0;

      CORE::LINALG::Matrix<ndim + 1, n> coord;

      CORE::LINALG::Matrix<n, 1> val;
      if (ndim == 1)
        CORE::FE::shape_function_1D(val, xi[0], distype);
      else if (ndim == 2)
        CORE::FE::shape_function_2D(val, xi[0], xi[1], distype);
      else if (ndim == 3)
        CORE::FE::shape_function_3D(val, xi[0], xi[1], xi[2], distype);

      else
        dserror("Wrong Dimension");

      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < globdim; ++j)
        {
          coord(j, i) = mynodes[i]->X()[j] + edisp[i * globdim + j];

          // use shape function values for interpolation
          globcoord[j] += val(i) * coord(j, i);
        }
      }

      return true;
    };
    /*----------------------------------------------------------------------*
     |  Get local coords for given global coords (ref position)  farah 01/14|
     *----------------------------------------------------------------------*/
    template <CORE::FE::CellType distype>
    void GlobalToLocal(DRT::Element& ele,  // element (input)
        double* xgl,                       // global position (input)
        double* xi)                        // local position  (output)
    {
      static constexpr int numnod = CORE::FE::num_nodes<distype>;
      static constexpr int ndim = CORE::FE::dim<distype>;

      CORE::LINALG::Matrix<numnod, 1> funct(true);
      CORE::LINALG::Matrix<ndim, numnod> xref(true);
      CORE::LINALG::Matrix<ndim, ndim> xjm(true);
      CORE::LINALG::Matrix<ndim, numnod> deriv(true);

      // spatial configuration of this element!
      for (int k = 0; k < numnod; ++k)
        for (int j = 0; j < ndim; ++j) xref(j, k) = ele.Nodes()[k]->X()[j];

      // first estimation for parameter space coordinates
      for (int p = 0; p < ndim; ++p)
      {
        if (distype == CORE::FE::CellType::quad4 or distype == CORE::FE::CellType::quad8 or
            distype == CORE::FE::CellType::quad9 or distype == CORE::FE::CellType::hex8 or
            distype == CORE::FE::CellType::hex20 or distype == CORE::FE::CellType::hex27)
        {
          xi[p] = 0.0;
        }
        else if (distype == CORE::FE::CellType::tri3 or distype == CORE::FE::CellType::tri6 or
                 distype == CORE::FE::CellType::tet4 or distype == CORE::FE::CellType::tet10)
        {
          xi[p] = 1.0 / 3.0;
        }
        else if (distype == CORE::FE::CellType::pyramid5)
        {
          if (p < 2)
            xi[p] = 0.0;
          else
            xi[p] = 1.0 / 3.0;
        }
        else
        {
          dserror("ERROR: Element type not supported for parameter space mapping!");
        }
      }

      double rhs[ndim];
      for (int p = 0; p < ndim; ++p) rhs[p] = 0.0;

      // converged
      bool converged = false;
      int j = 0;

      //************************************************
      // loop
      while (!converged and j < 10)
      {
        // reset matriced
        xjm.Clear();
        deriv.Clear();

        if (ndim == 2)
        {
          CORE::FE::shape_function_2D(funct, xi[0], xi[1], distype);
          CORE::FE::shape_function_2D_deriv1(deriv, xi[0], xi[1], distype);
        }
        else if (ndim == 3)
        {
          CORE::FE::shape_function_3D(funct, xi[0], xi[1], xi[2], distype);
          CORE::FE::shape_function_3D_deriv1(deriv, xi[0], xi[1], xi[2], distype);
        }
        else
          dserror("ERROR");

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p)
            for (int l = 0; l < ndim; ++l) xjm(p, l) += deriv(l, k) * xref(p, k);

        // rhs of (linearized equation)
        for (int p = 0; p < ndim; ++p) rhs[p] = 0.0;

        for (int p = 0; p < ndim; ++p) rhs[p] = -xgl[p];

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p) rhs[p] += funct(k) * xref(p, k);

        double norm = 0.0;
        for (int p = 0; p < ndim; ++p) norm += rhs[p] * rhs[p];

        if (sqrt(norm) < 1e-12)
        {
          converged = true;
        }
        //     else
        //     {
        //       std::cout << "norm= " << norm << " xi[0]=  " << xi[0]<< " xi[1]=  " << xi[1]<< "
        //       xi[2]=  " << xi[2] << std::endl;
        //     }

        // solve equation
        if (abs(xjm.Determinant()) < 1e-15)
        {
          dserror("*** WARNING: jacobi singular ***");
        }

        double xjm_invert = xjm.Invert();
        if (abs(xjm_invert) < 1e-12) dserror("ERROR: Singular Jacobian for advection map");

        double deltaxi[ndim];
        for (int p = 0; p < ndim; ++p) deltaxi[p] = 0.0;

        for (int z = 0; z < ndim; ++z)
          for (int p = 0; p < ndim; ++p) deltaxi[z] -= xjm(z, p) * rhs[p];

        // incremental update
        for (int p = 0; p < ndim; ++p) xi[p] += deltaxi[p];

        j = j + 1;
      }  // end loop
      //************************************************

      if (!converged) dserror("Evaluation of element coordinates not converged!");

      return;
    };

    /*----------------------------------------------------------------------*
     |  Get local coords for given global coords (ref position)  farah 01/14|
     *----------------------------------------------------------------------*/
    template <CORE::FE::CellType distype>
    void GlobalToLocal(DRT::Element& ele,  // element (input)
        double* xgl,                       // global position (input)
        double* xi,
        bool& converged)  // converged solution ?
    {
      static constexpr int numnod = CORE::FE::num_nodes<distype>;
      static constexpr int ndim = CORE::FE::dim<distype>;

      // converged
      converged = false;

      CORE::LINALG::Matrix<numnod, 1> funct;
      CORE::LINALG::Matrix<ndim, numnod> xref;
      CORE::LINALG::Matrix<ndim, ndim> xjm;
      CORE::LINALG::Matrix<ndim, numnod> deriv;

      // spatial configuration of this element!
      for (int k = 0; k < numnod; ++k)
        for (int j = 0; j < ndim; ++j) xref(j, k) = ele.Nodes()[k]->X()[j];

      // first estimation for parameter space coordinates
      for (int p = 0; p < ndim; ++p)
      {
        if (distype == CORE::FE::CellType::quad4 or distype == CORE::FE::CellType::quad8 or
            distype == CORE::FE::CellType::quad9 or distype == CORE::FE::CellType::hex8 or
            distype == CORE::FE::CellType::hex20 or distype == CORE::FE::CellType::hex27)
        {
          xi[p] = 0.0;
        }
        else if (distype == CORE::FE::CellType::tri3 or distype == CORE::FE::CellType::tri6 or
                 distype == CORE::FE::CellType::tet4 or distype == CORE::FE::CellType::tet10)
        {
          xi[p] = 1.0 / 3.0;
        }
        else if (distype == CORE::FE::CellType::pyramid5)
        {
          if (p < 2)
            xi[p] = 0.0;
          else
            xi[p] = 1.0 / 3.0;
        }
        else
        {
          dserror("ERROR: Element type not supported for parameter space mapping!");
        }
      }
      double rhs[ndim];

      int j = 0;

      //************************************************
      // loop
      while (!converged and j < 10)
      {
        // reset matriced
        xjm.Clear();
        deriv.Clear();

        if (ndim == 2)
        {
          CORE::FE::shape_function_2D(funct, xi[0], xi[1], distype);
          CORE::FE::shape_function_2D_deriv1(deriv, xi[0], xi[1], distype);
        }
        else if (ndim == 3)
        {
          CORE::FE::shape_function_3D(funct, xi[0], xi[1], xi[2], distype);
          CORE::FE::shape_function_3D_deriv1(deriv, xi[0], xi[1], xi[2], distype);
        }
        else
          dserror("ERROR");

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p)
            for (int l = 0; l < ndim; ++l) xjm(p, l) += deriv(l, k) * xref(p, k);

        // rhs of (linearized equation)
        for (int p = 0; p < ndim; ++p) rhs[p] = 0.0;

        for (int p = 0; p < ndim; ++p) rhs[p] = -xgl[p];

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p) rhs[p] += funct(k) * xref(p, k);

        double norm = 0.0;
        for (int p = 0; p < ndim; ++p) norm += rhs[p] * rhs[p];

        if (sqrt(norm) < 1e-12) converged = true;

        if (converged == true) break;

        // solve equation
        if (abs(xjm.Determinant()) < 1e-15)
        {
          std::cout << "WARNING !!! jacobi determinant singular! In GlobalToLocal(...)"
                    << std::endl;
          std::cout << "JAC= " << xjm.Determinant() << std::endl;
          std::cout << "CONVERGED= " << converged << std::endl;
          //       dserror("*** WARNING: jacobi singular ***");
          converged = false;
          break;
        }

        double xjm_invert = xjm.Invert();
        if (abs(xjm_invert) < 1e-15) dserror("ERROR: Singular Jacobian");

        double deltaxi[3];
        for (int p = 0; p < ndim; ++p) deltaxi[p] = 0.0;

        for (int z = 0; z < ndim; ++z)
          for (int p = 0; p < ndim; ++p) deltaxi[z] -= xjm(z, p) * rhs[p];

        // incremental update
        for (int p = 0; p < ndim; ++p) xi[p] += deltaxi[p];

        j = j + 1;
      }  // end loop
      //************************************************

      return;
    };

    /*----------------------------------------------------------------------*
     |  Get local coords for given global coords (curr position) rauch 08/14|
     *----------------------------------------------------------------------*/
    template <CORE::FE::CellType distype>
    void GlobalToCurrentLocal(DRT::Element& ele,  // element (input)
        const double* targetdisp,
        double* xgl,  // global position (input)
        double* xi, bool& converged,
        double& residual)  // converged solution ?
    {
      static constexpr int numnod = CORE::FE::num_nodes<distype>;
      static constexpr int ndim = CORE::FE::dim<distype>;

      // converged
      converged = false;

      CORE::LINALG::Matrix<numnod, 1> funct;
      CORE::LINALG::Matrix<ndim, numnod> xref;
      CORE::LINALG::Matrix<ndim, ndim> xjm;
      CORE::LINALG::Matrix<ndim, numnod> deriv;

      // spatial configuration of this element!
      for (int k = 0; k < numnod; ++k)
        for (int j = 0; j < ndim; ++j)
          xref(j, k) = (ele.Nodes()[k]->X()[j]) + targetdisp[k * ndim + j];

      // first estimation for parameter space coordinates
      for (int p = 0; p < 3; ++p) xi[p] = 0.0;

      double rhs[ndim];

      int j = 0;

      //************************************************
      // loop
      while (!converged and j < 50)
      {
        // reset matrices
        xjm.Clear();
        deriv.Clear();

        if (ndim == 2)
        {
          CORE::FE::shape_function_2D(funct, xi[0], xi[1], distype);
          CORE::FE::shape_function_2D_deriv1(deriv, xi[0], xi[1], distype);
        }
        else if (ndim == 3)
        {
          CORE::FE::shape_function_3D(funct, xi[0], xi[1], xi[2], distype);
          CORE::FE::shape_function_3D_deriv1(deriv, xi[0], xi[1], xi[2], distype);
        }
        else
          dserror("ERROR");

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p)
            for (int l = 0; l < ndim; ++l) xjm(p, l) += deriv(l, k) * xref(p, k);

        for (int p = 0; p < ndim; ++p) rhs[p] = -xgl[p];

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p) rhs[p] += funct(k) * xref(p, k);

        double norm = 0.0;
        for (int p = 0; p < ndim; ++p) norm += rhs[p] * rhs[p];
        residual = sqrt(norm);
        if (residual < 1e-13) converged = true;

        if (converged == true) break;

        // solve equation
        if (abs(xjm.Determinant()) < 1e-15)
        {
          std::cout << "WARNING !!! jacobi determinant singular! In GlobalToCurrentLocal(...)"
                    << std::endl;
          std::cout << "JAC= " << xjm.Determinant() << std::endl;
          std::cout << "CONVERGED= " << converged << std::endl;
          //       dserror("*** WARNING: jacobi singular ***");
          converged = false;
          break;
        }

        double xjm_invert = xjm.Invert();
        if (abs(xjm_invert) < 1e-15) dserror("ERROR: Singular Jacobian");

        double deltaxi[3];
        for (int p = 0; p < ndim; ++p) deltaxi[p] = 0.0;

        for (int z = 0; z < ndim; ++z)
          for (int p = 0; p < ndim; ++p) deltaxi[z] -= xjm(z, p) * rhs[p];

        // incremental update
        for (int p = 0; p < ndim; ++p) xi[p] += deltaxi[p];

        j = j + 1;
      }  // end loop
      //************************************************

      return;
    };

  }  // namespace UTILS
}  // namespace MORTAR



BACI_NAMESPACE_CLOSE

#endif