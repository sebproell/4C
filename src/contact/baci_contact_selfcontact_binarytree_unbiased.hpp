/*-----------------------------------------------------------------------*/
/*! \file
\brief Search tree for unbiased self-contact problems

\level 2

*/
/*-----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_SELFCONTACT_BINARYTREE_UNBIASED_HPP
#define FOUR_C_CONTACT_SELFCONTACT_BINARYTREE_UNBIASED_HPP

#include "baci_config.hpp"

#include "baci_contact_selfcontact_binarytree.hpp"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
  \brief A class for performing unbiased self contact search in 2D / 3D based on a binary search
  tree and processor specific dual graphs to enhance parallel search

  */

  class UnbiasedSelfBinaryTree : public SelfBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    \param discret (in):    The contact interface discretization
    \param iparams (in):    interface specific parameter list
    \param elements (in):   All elements on self contact interface (fully overlapping map)
    \param dim (in):        The problem dimension

    */
    UnbiasedSelfBinaryTree(DRT::Discretization& discret, const Teuchos::ParameterList& iparams,
        Teuchos::RCP<Epetra_Map> elements, int dim, double eps);


    /*!
    \brief Initialize the unbiased self binary tree

    */
    void Init() final;

   private:
    //! @name Evaluation methods

    /*!
    \brief Decide whether tree nodes need to be added to contact pairs

    \param [in]  treenode1:  first tree node
    \param [in]  treenode2:  second tree node

    */
    void AddTreeNodesToContactPairs(Teuchos::RCP<SelfBinaryTreeNode> treenode1,
        Teuchos::RCP<SelfBinaryTreeNode> treenode2) final;

    /*!
    \brief Calculate the processor specific dual graph

    \param [out] dualgraph:  construction of binary tree is based on this data
    \param [in]  elelist:    list (gids) of all contact elements of the surface
    \param [in]  p:          number of current (not necessarily calling) processor

    */
    void CalculateProcSpecificDualGraph(
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualGraph,
        const std::vector<int>& elelist, const int p);

    /*!
    \brief Use list of possible contact pairs to define the search elements

    */
    void DefineSearchElements();

    /*!
    \brief Get the mortar element specific number of first order nodes

    \param [in]      contractedEdge:  dual edge that is contracted
    \param [in,out]  contractedNode:  node that consists of both nodes of contracted edge

    */
    void GetContractedNode(Teuchos::RCP<SelfDualEdge>& contractedEdge,
        Teuchos::RCP<SelfBinaryTreeNode>& contractedNode) final;

    /*!
    \brief Initialize Tree in a bottom up way based on dual graph

    \param [in] procdualgraph:  processor specific dual graph, construction of binary tree is based
    on this data

    */
    void InitializeTreeBottomUp(std::map<int,
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>>*
            procdualGraph);

    /*!
    \brief Checks roughly whether self contact of two elements shall be evaluated (3D)

    This method checks if the normal at the slave element center and the vector connecting this
    slave element center and the center of the element it is projected to (master element) point in
    the same direction. All these vectors are evaluated in the reference configuration to check the
    initial state. Only if both vectors point in the same direction integration shall be performed.

    */
    bool RoughCheckRefConfig(int ele1gid, int ele2gid);

    /*!
    \brief Evaluate unbiased binary search tree

    */
    void SearchContact() final;

    /*!
    \brief Communicate the Search Elements to all processors

    */
    void CommunicateSearchElementsAllProcs();
    //@}

    // don't want = operator and cctor
    UnbiasedSelfBinaryTree operator=(const SelfBinaryTree& old) = delete;
    UnbiasedSelfBinaryTree(const SelfBinaryTree& old) = delete;

    //! use two half pass approach
    const bool Two_half_pass_;
    //! perform reference configuration check for non-smooth self contact
    const bool Check_nonsmooth_selfcontactsurface_;
    //! the contact pairs are communicated to all processors
    const bool Searchele_AllProc_;
  };  // class UnbiasedSelfBinaryTree
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif