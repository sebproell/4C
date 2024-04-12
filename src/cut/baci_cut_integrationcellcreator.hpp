/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells for the Tessellation routine

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_INTEGRATIONCELLCREATOR_HPP
#define FOUR_C_CUT_INTEGRATIONCELLCREATOR_HPP

#include "baci_config.hpp"

#include "baci_cut_volumecell.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class Mesh;

    /*!
    \brief Special cases for integration cell creation: a library of cuts.

    Some cuts lead to known shapes those meshes are known. Here we collect
    these cases. Feel free to add your own.
   */
    class IntegrationCellCreator
    {
     public:
      static bool CreateCells(Mesh& mesh, Element* element, const plain_volumecell_set& cells);

      static bool CreateCell(Mesh& mesh, CORE::FE::CellType shape, VolumeCell* cell);

     private:
      /** \brief loop over all volume cells and initiate the volume and boundary integration cell
       *  creation process
       *
       *  This is done after a successful initiation of the volume cells.*/
      void Execute(Mesh& mesh);

      /// fill the point1 boundary cell
      bool CreatePoint1Cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the line2 volume cell
      bool CreateLine2Cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /** \brief fill the 2-D volume cells
       *
       *  \author hiermeier \date 01/17 */
      template <CORE::FE::CellType celltype,
          CORE::FE::CellType facetype = CORE::FE::DisTypeToFaceShapeType<celltype>::shape,
          unsigned numfaces = CORE::FE::num_faces<celltype>>
      bool Create2DCell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the tet4 volume cell
      bool CreateTet4Cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the hex8 volume cell
      bool CreateHex8Cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the wedge6 volume cell
      bool CreateWedge6Cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the pyramid5 volume cell
      bool CreatePyramid5Cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      bool CreateSpecialCases(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      bool Hex8HorizontalCut(Mesh& mesh, Element* element, VolumeCell* cell,
          const plain_facet_set& facets, int axis, double r);

      /// add the volume cell information during the Create<Shape>Cell calls
      void Add(VolumeCell* vc, CORE::FE::CellType shape, const std::vector<Point*>& points)
      {
        volume& v = cells_[vc];
        std::vector<ic>& cells = v.domain_;
        cells.push_back(ic());
        ic& cell = cells.back();
        cell.shape_ = shape;
        cell.points_.reserve(points.size());
        cell.points_.assign(points.begin(), points.end());
      }

      /** \brief add the side (boundary cell) information during Create<Shape>Cell calls
       *         for the desired boundary cell positions
       *
       *  \author hiermeier \date 01/17   */
      void AddSide(INPAR::CUT::BoundaryCellPosition bcell_position, VolumeCell* vc, Facet* facet,
          CORE::FE::CellType shape, const std::vector<Point*>& side);

      /// add the side (boundary cell) information during the Create<Shape>Cell calls
      void AddSide(
          VolumeCell* vc, Facet* facet, CORE::FE::CellType shape, const std::vector<Point*>& side)
      {
        volume& cell = cells_[vc];
        std::vector<bc>& bcells = cell.boundary_;
        bcells.push_back(bc());
        bc& bcell = bcells.back();
        bcell.shape_ = shape;
        bcell.facet_ = facet;
        bcell.side_.reserve(side.size());
        bcell.side_.assign(side.begin(), side.end());
      }

      /*-------------------------------------------------------------------------*/
      /// construction of a boundary integration cell
      struct bc
      {
        CORE::FE::CellType shape_;
        std::vector<Point*> side_;
        Facet* facet_;

        /// actual creation of the desired boundary integration cell
        void Execute(Mesh& mesh, VolumeCell* vc)
        {
          vc->NewBoundaryCell(mesh, shape_, facet_, side_);
        }
      };

      /*-------------------------------------------------------------------------*/
      /// construction of a volume integration cell
      struct ic
      {
        CORE::FE::CellType shape_;
        std::vector<Point*> points_;

        /// actual creation of the desired volume integration cell
        void Execute(Mesh& mesh, VolumeCell* vc) { vc->NewIntegrationCell(mesh, shape_, points_); }
      };

      /*-------------------------------------------------------------------------*/
      /** \brief Create all volume and boundary integration cells                */
      struct volume
      {
        std::vector<ic> domain_;
        std::vector<bc> boundary_;

        /** \brief This routine initiates the volume and boundary cell creation.
         *
         *  At this point the necessary information for the creation must have been
         *  already added by the Add() and AddSide() methods. */
        void Execute(Mesh& mesh, VolumeCell* vc)
        {
          for (std::vector<ic>::iterator i = domain_.begin(); i != domain_.end(); ++i)
          {
            ic& cell = *i;
            cell.Execute(mesh, vc);
          }
          for (std::vector<bc>::iterator i = boundary_.begin(); i != boundary_.end(); ++i)
          {
            bc& bcell = *i;
            bcell.Execute(mesh, vc);
          }
        }
      };  // struct volume

      std::map<VolumeCell*, volume> cells_;
    };

  }  // namespace CUT
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif