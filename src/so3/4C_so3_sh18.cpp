/*----------------------------------------------------------------------*/
/*! \file

\brief ToDo Add meaningful comment.

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_sh18.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::SoSh18Type DRT::ELEMENTS::SoSh18Type::instance_;

DRT::ELEMENTS::SoSh18Type& DRT::ELEMENTS::SoSh18Type::Instance() { return instance_; }
namespace
{
  const std::string name = DRT::ELEMENTS::SoSh18Type::Instance().Name();
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoSh18Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoSh18(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoSh18Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoSh18(id, owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoSh18Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoSh18(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoSh18Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX18"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("HEX18", 18)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedString("TSL")
                      .AddNamedString("MEL")
                      .AddNamedString("CTL")
                      .AddNamedString("VOL")
                      .AddOptionalNamedDoubleVector("RAD", 3)
                      .AddOptionalNamedDoubleVector("AXI", 3)
                      .AddOptionalNamedDoubleVector("CIR", 3)
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .AddOptionalNamedDoubleVector("FIBER2", 3)
                      .AddOptionalNamedDoubleVector("FIBER3", 3)
                      .AddOptionalNamedDouble("STRENGTH")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoSh18::SoSh18(int id, int owner) : SoBase(id, owner), SoHex18(id, owner)
{
  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoSh18::SoSh18(const DRT::ELEMENTS::SoSh18& old)
    : SoBase(old),
      SoHex18(old),
      dsg_shear_(old.dsg_shear_),
      dsg_membrane_(old.dsg_membrane_),
      dsg_ctl_(old.dsg_ctl_),
      eas_(old.eas_)
{
  SetupDSG();
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoSh18::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoSh18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh18::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  SoBase::Pack(data);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  // element technology bools
  AddtoPack(data, (int)dsg_shear_);
  AddtoPack(data, (int)dsg_membrane_);
  AddtoPack(data, (int)dsg_ctl_);
  AddtoPack(data, (int)eas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh18::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  SoBase::Unpack(basedata);

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, CORE::LINALG::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  // element technology bools
  dsg_shear_ = ExtractInt(position, data);
  dsg_membrane_ = ExtractInt(position, data);
  dsg_ctl_ = ExtractInt(position, data);
  eas_ = ExtractInt(position, data);
  SetupDSG();

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh18::Print(std::ostream& os) const
{
  os << "So_sh18 ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE