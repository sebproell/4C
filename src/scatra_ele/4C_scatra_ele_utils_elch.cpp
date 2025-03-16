// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_utils_elch.hpp"

#include "4C_fem_condition.hpp"
#include "4C_mat_ion.hpp"
#include "4C_scatra_ele_calc_elch.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 08/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::ScaTraEleUtilsElch<distype>*
Discret::Elements::ScaTraEleUtilsElch<distype>::instance(
    const int numdofpernode,    ///< number of degrees of freedom per node
    const int numscal,          ///< number of transported scalars per node
    const std::string& disname  ///< name of discretization
)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleUtilsElch<distype>>(
            new ScaTraEleUtilsElch<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].instance(
      Core::Utils::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 07/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::ScaTraEleUtilsElch<distype>::ScaTraEleUtilsElch(
    const int numdofpernode,    ///< number of degrees of freedom per node
    const int numscal,          ///< number of transported scalars per node
    const std::string& disname  ///< name of discretization
    )
    : numdofpernode_(numdofpernode), numscal_(numscal)
{
  return;
}


/*---------------------------------------------------------------------------------------------------------*
 | evaluation of electrochemistry kinetics at integration point on domain or boundary element   fang
 07/15 |
 *---------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::ScaTraEleUtilsElch<distype>::evaluate_elch_kinetics_at_integration_point(
    const Core::Elements::Element* ele,                        ///< current element
    Core::LinAlg::SerialDenseMatrix& emat,                     ///< element matrix
    Core::LinAlg::SerialDenseVector& erhs,                     ///< element right-hand side vector
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp,  ///< state variables at element nodes
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>&
        ehist,                                   ///< history variables at element nodes
    const double timefac,                        ///< time factor
    const double fac,                            ///< Gauss integration factor
    const Core::LinAlg::Matrix<nen_, 1>& funct,  ///< shape functions at int. point
    Core::Conditions::Condition& cond,           ///< condition
    const int nume,                              ///< number of transferred electrons
    const std::vector<int>& stoich,              ///< stoichiometry of the reaction
    const double valence_k,                      ///< valence of the single reactant
    const int kinetics,                          ///< desired electrode kinetics model
    const double pot0,                           ///< actual electrode potential on metal side
    const double frt,                            ///< factor F/RT
    const double fns,                            ///< factor fns = s_k / (nume * faraday * (-1))
    const double
        scalar,  ///< scaling factor for element matrix and right-hand side vector contributions
    const int k  ///< index of evaluated scalar
) const
{
  // for pre-multiplication of i0 with 1/(F z_k)
  const double faraday = Discret::Elements::ScaTraEleParameterElch::instance("scatra")->faraday();

  // concentration of active species at integration point
  std::vector<double> conint(numscal_, 0.0);

  // concentration is evaluated at all GP since some reaction models depend on all concentrations
  for (int kk = 0; kk < numscal_; ++kk) conint[kk] = funct.dot(ephinp[kk]);

  // el. potential at integration point
  const double potint = funct.dot(ephinp[numscal_]);

  // history of potential on electrode boundary at integration point
  const double pothistint = funct.dot(ehist[numscal_]);

  // electrode potential difference (epd) at integration point
  const double epd = (pot0 - potint);

  // concentration-dependent Butler-Volmer law(s)
  switch (kinetics)
  {
    case Inpar::ElCh::butler_volmer:
    case Inpar::ElCh::butler_volmer_yang1997:
    {
      // read model-specific parameters
      const auto alphaa = cond.parameters().get<double>("ALPHA_A");
      const auto alphac = cond.parameters().get<double>("ALPHA_C");
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      double pot0hist = 0.0;
      if (dlcap != 0.0) pot0hist = cond.parameters().get<double>("pot0hist");
      auto i0 = cond.parameters().get<double>("I0");
      if (i0 < -1e-14)
        FOUR_C_THROW(
            "i0 is negative, \n"
            "a positive definition is necessary due to generalized reaction models: {}",
            i0);
      // add time factor
      i0 *= timefac;
      const auto gamma = cond.parameters().get<double>("GAMMA");
      const auto refcon = cond.parameters().get<double>("REFCON");
      if (refcon < 1e-12) FOUR_C_THROW("reference concentration is too small: {}", refcon);

      if (valence_k != nume)
        FOUR_C_THROW(
            "Kinetic model Butler-Volmer: The number of transferred electrons need to be  \n "
            "the same as the charge number of the reacting species {}",
            k);

      // open circuit potential is assumed to be zero
      const double ocp = 0.0;

      // overpotential based on opencircuit potential
      const double eta = epd - ocp;

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // some safety checks/ user warnings
      if ((alphaa * frt * eta) > 100.0)
        std::cout << "WARNING: Exp(alpha_a...) in Butler-Volmer law is near overflow!"
                  << exp(alphaa * frt * eta) << std::endl;
      if (((-alphac) * frt * eta) > 100.0)
        std::cout << "WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
                  << exp((-alphac) * frt * eta) << std::endl;
#endif
      double pow_conint_gamma_k = 0.0;
      if ((conint[k] / refcon) < 1e-13)
      {
        pow_conint_gamma_k = std::pow(1e-13, gamma);
#ifdef FOUR_C_ENABLE_ASSERTIONS
        std::cout << "WARNING: Rel. Conc. in Butler-Volmer formula is zero/negative: "
                  << (conint[k] / refcon) << std::endl;
        std::cout << "-> Replacement value: pow(EPS,gamma) = " << pow_conint_gamma_k << std::endl;
#endif
      }
      else
        pow_conint_gamma_k = std::pow(conint[k] / refcon, gamma);

      if (kinetics == Inpar::ElCh::butler_volmer)
      {
        // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
        const double expterm = exp(alphaa * frt * eta) - exp((-alphac) * frt * eta);

        double concterm = 0.0;
        if (conint[k] > 1e-13)
          concterm = gamma * pow(conint[k], (gamma - 1.0)) / pow(refcon, gamma);
        else
          concterm = gamma * pow(1e-13, (gamma - 1.0)) / pow(refcon, gamma);

        for (int vi = 0; vi < nen_; ++vi)
        {
          const double fac_fns_i0_funct_vi = scalar * fac * fns * i0 * funct(vi);

          // ------matrix: d(R_k)/dx = d(theta*dt*(-1)*(w_k,j_k))/dx
          for (int ui = 0; ui < nen_; ++ui)
          {
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + k) +=
                -fac_fns_i0_funct_vi * concterm * funct(ui) * expterm;
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
                -fac_fns_i0_funct_vi * pow_conint_gamma_k *
                (((-alphaa) * frt * exp(alphaa * frt * eta)) +
                    ((-alphac) * frt * exp((-alphac) * frt * eta))) *
                funct(ui);
          }

          // -----right-hand-side: -R_k = -theta*dt*(-1)*(w_k,j_k)
          erhs[vi * numdofpernode_ + k] -= -fac_fns_i0_funct_vi * pow_conint_gamma_k * expterm;
        }

        if (dlcap != 0.0)
        {
          for (int vi = 0; vi < nen_; ++vi)
          {
            // TODO: Do we need epsilon here
            // add terms of double layer capacitance current density
            for (int ui = 0; ui < nen_; ++ui)
              emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
                  fac * funct(vi) * funct(ui) * dlcap / (nume * faraday);

            // -----right-hand-side: -R_k = -(theta*dt*(-1)*(w_k,j_k)
            erhs[vi * numdofpernode_ + k] += fac * funct(vi) * dlcap / (nume * faraday) *
                                             (pothistint - pot0hist - potint + pot0);
          }
        }
      }  // end if(kinetics=="Butler-Volmer")

      else if (kinetics == Inpar::ElCh::butler_volmer_yang1997)
      {
        if (dlcap != 0.0)
          FOUR_C_THROW(
              "double layer charging is not implemented for Butler-Volmer-Yang electrode kinetics");

        // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
        double concterm = 0.0;
        if ((conint[k] / refcon) > 1e-13)
          concterm = gamma * pow(conint[k], (gamma - 1.0)) / pow(refcon, gamma);
        else
          concterm = gamma * pow(1e-13, (gamma - 1.0)) / pow(refcon, gamma);

        for (int vi = 0; vi < nen_; ++vi)
        {
          const double fac_fns_i0_funct_vi = scalar * fac * fns * i0 * funct(vi);
          // ------matrix: d(R_k)/dx = d(theta*dt*(-1)*(w_k,j_k))/dx
          for (int ui = 0; ui < nen_; ++ui)
          {
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + k) +=
                -fac_fns_i0_funct_vi * funct(ui) * (-(concterm * exp((-alphac) * frt * eta)));
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
                -fac_fns_i0_funct_vi *
                (((-alphaa) * frt * exp(alphaa * frt * eta)) +
                    (pow_conint_gamma_k * (-alphac) * frt * exp((-alphac) * frt * eta))) *
                funct(ui);
          }

          // -----right-hand-side: -R_k = -theta*dt*(-1)*(w_k,j_k)
          erhs[vi * numdofpernode_ + k] -=
              -fac_fns_i0_funct_vi *
              (exp(alphaa * frt * eta) - (pow_conint_gamma_k * exp((-alphac) * frt * eta)));
        }
      }  // if (kinetics=="Butler-Volmer-Yang1997")
      else
        FOUR_C_THROW(
            "You should not be here!! Two options: Butler-Volmer-Yang1997 and "
            "Butler-Volmer-Yang1997 ");

      break;
    }

    // Tafel law (see phd-thesis Georg Bauer, pp.25):
    // implementation of cathodic path: i_n = i_0 * (-exp(-alpha * frt* eta)
    // -> cathodic reaction path: i_0 > 0 and alpha > 0
    // -> anodic reacton path:    i_0 < 0 and alpha < 0
    case Inpar::ElCh::tafel:
    {
      // read model-specific parameter
      const auto alpha = cond.parameters().get<double>("ALPHA");
      auto i0 = cond.parameters().get<double>("I0");
      i0 *= timefac;
      const auto gamma = cond.parameters().get<double>("GAMMA");
      const auto refcon = cond.parameters().get<double>("REFCON");
      if (refcon < 1e-12) FOUR_C_THROW("reference concentration is too small: {}", refcon);
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      if (dlcap != 0.0)
        FOUR_C_THROW("double layer charging is not implemented for Tafel electrode kinetics");

      if (valence_k != nume)
        FOUR_C_THROW(
            "Kinetic model Butler-Volmer: The number of transferred electrons need to be  \n "
            "the same as the charge number of the reacting species {}",
            k);
      // opencircuit potential is assumed to be zero
      const double ocp = 0.0;
      // overpotential based on opencircuit potential
      const double eta = epd - ocp;

      // concentration-dependent Tafel law
      double pow_conint_gamma_k(0.0);

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // some safety checks/ user warnings
      if (((-alpha) * frt * eta) > 100.0)
        std::cout << "WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
                  << exp((-alpha) * frt * eta) << std::endl;
#endif
      if ((conint[k] / refcon) < 1e-13)
      {
        pow_conint_gamma_k = std::pow(1e-13, gamma);
#ifdef FOUR_C_ENABLE_ASSERTIONS
        std::cout << "WARNING: Rel. Conc. in Tafel formula is zero/negative: "
                  << (conint[k] / refcon) << std::endl;
        std::cout << "-> Replacement value: pow(EPS,gamma) = " << pow_conint_gamma_k << std::endl;
#endif
      }
      else
        pow_conint_gamma_k = std::pow(conint[k] / refcon, gamma);

      const double expterm = -exp((-alpha) * frt * eta);

      double concterm = 0.0;
      // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
      if (conint[k] > 1e-13)
        concterm = gamma * pow(conint[k], (gamma - 1.0)) / pow(refcon, gamma);
      else
        concterm = gamma * pow(1e-13, (gamma - 1.0)) / pow(refcon, gamma);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double fac_fns_i0_funct_vi = scalar * fac * fns * i0 * funct(vi);
        // ---------------------matrix
        for (int ui = 0; ui < nen_; ++ui)
        {
          emat(vi * numdofpernode_ + k, ui * numdofpernode_ + k) +=
              -fac_fns_i0_funct_vi * concterm * funct(ui) * expterm;
          emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
              -fac_fns_i0_funct_vi * pow_conint_gamma_k * (-alpha) * frt *
              exp((-alpha) * frt * eta) *
              funct(ui);  // do not forget the (-1) from differentiation of eta!
        }
        // ------------right-hand-side
        erhs[vi * numdofpernode_ + k] -= -fac_fns_i0_funct_vi * pow_conint_gamma_k * expterm;
      }

      break;
    }

    // linear law:  i_n = frt*i_0*((alphaa+alpha_c)*(V_M - phi)) -> changed 10/13
    // previously implemented: i_n = i_0*(alphaa*frt*(V_M - phi) + 1.0)
    //                         -> linearization in respect to anodic branch!!
    //                         this is not the classical version of a linear electrode kinetics law
    case Inpar::ElCh::linear:
    {
      // read model-specific parameter
      const auto alphaa = cond.parameters().get<double>("ALPHA");
      auto i0 = cond.parameters().get<double>("I0");
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      double pot0hist = 0.0;
      if (dlcap != 0.0) pot0hist = cond.parameters().get<double>("pot0hist");
      i0 *= timefac;
      if (i0 < -1e-14)
        FOUR_C_THROW(
            "i0 is negative, \n"
            "a positive definition is necessary due to generalized reaction models: {}",
            i0);
      const auto gamma = cond.parameters().get<double>("GAMMA");
      const auto refcon = cond.parameters().get<double>("REFCON");
      if (refcon < 1e-12) FOUR_C_THROW("reference concentration is too small: {}", refcon);

      if (valence_k != nume)
        FOUR_C_THROW(
            "Kinetic model Butler-Volmer: The number of transferred electrons need to be  \n "
            "the same as the charge number of the reacting species {}",
            k);
      // opencircuit potential is assumed to be zero
      const double ocp = 0.0;
      // overpotential based on opencircuit potential
      const double eta = epd - ocp;

      double pow_conint_gamma_k = 0.0;
      if ((conint[k] / refcon) < 1e-13)
      {
        pow_conint_gamma_k = std::pow(1e-13, gamma);
#ifdef FOUR_C_ENABLE_ASSERTIONS
        std::cout << "WARNING: Rel. Conc. in Tafel formula is zero/negative: "
                  << (conint[k] / refcon) << std::endl;
        std::cout << "-> Replacement value: pow(EPS,gamma) = " << pow_conint_gamma_k << std::endl;
#endif
      }
      else
        pow_conint_gamma_k = std::pow(conint[k] / refcon, gamma);
      const double linearfunct = (alphaa * frt * eta);
      // note: gamma==0 deactivates concentration dependency
      double concterm = 0.0;
      if (conint[k] > 1e-13)
        concterm = gamma * pow(conint[k], (gamma - 1.0)) / pow(refcon, gamma);
      else
        FOUR_C_THROW("Better stop here!");

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double fac_fns_i0_funct_vi = scalar * fac * fns * i0 * funct(vi);
        const int fvi = vi * numdofpernode_ + k;
        // ---------------------matrix
        for (int ui = 0; ui < nen_; ++ui)
        {
          emat(fvi, ui * numdofpernode_ + k) +=
              -fac_fns_i0_funct_vi * concterm * funct(ui) * linearfunct;
          emat(fvi, ui * numdofpernode_ + numscal_) +=
              -fac_fns_i0_funct_vi * pow_conint_gamma_k * (-alphaa) * frt * funct(ui);
        }
        // ------------right-hand-side
        erhs[fvi] -= -fac_fns_i0_funct_vi * pow_conint_gamma_k * linearfunct;
      }

      if (dlcap != 0.0)
      {
        for (int vi = 0; vi < nen_; ++vi)
        {  // add terms of double layer capacitance current density
          for (int ui = 0; ui < nen_; ++ui)
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
                fac * funct(vi) * funct(ui) * dlcap / (nume * faraday);

          // -----right-hand-side: -R_k = -(theta*dt*(-1)*(w_k,j_k)
          erhs[vi * numdofpernode_ + k] +=
              fac * funct(vi) * dlcap / (nume * faraday) * (pothistint - pot0hist - potint + pot0);
        }
      }

      break;
    }

    case Inpar::ElCh::butler_volmer_newman:
    {
      // "Electrochemical systems"
      // Newman and Thomas-Alyea, 2004
      // General stoichiometry: pp. 212-213, e.q. 8.26
      // consideration of a elementary step of the form:
      // Sum_i s_i M_i ->  ne-
      // n is one if charge transfer is involved, multiple electron transfers "being unlikely in
      // an elementary step

      const auto k_a = cond.parameters().get<double>("K_A");
      const auto k_c = cond.parameters().get<double>("K_C");
      const auto beta = cond.parameters().get<double>("BETA");
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      if (dlcap != 0.0)
        FOUR_C_THROW(
            "double layer charging is not implemented for Butler-Volmer-Newman electrode kinetics");

      // reaction order of the cathodic and anodic reactants of ionic species k
      std::vector<int> q(numscal_, 0);
      std::vector<int> p(numscal_, 0);

      for (int ii = 0; ii < numscal_; ii++)
      {
        // according to the convention: anodic reactant is positive
        if (stoich[ii] > 0)
        {
          q[ii] = 0;
          p[ii] = stoich[ii];
        }
        // according to the convention: cathodic reactant is negative
        else
        {
          q[ii] = -stoich[ii];
          p[ii] = 0;
        }
      }

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // some safety checks/ user warnings
      if (((1 - beta) * frt * epd) > 100.0)
        std::cout << "WARNING: Exp((1-beta)...) in Butler-Volmer law is near overflow!"
                  << exp((1 - beta) * frt * epd) << std::endl;
      if (((-beta) * frt * epd) > 100.0)
        std::cout << "WARNING: Exp(-beta...) in Butler-Volmer law is near overflow!"
                  << exp((-beta) * frt * epd) << std::endl;
#endif

      double pow_conint_p(1.0);  // product over i (c_i)^(p_i)
      double pow_conint_q(1.0);  // product over i (c_i)^(q_i)
      std::vector<double> pow_conint_p_derivative(
          numscal_, 1.0);  // pow_conint_p derived after conint[nspec]
      std::vector<double> pow_conint_q_derivative(
          numscal_, 1.0);  // pow_conint_q derived after conint[nspec]

      // concentration term (product of cathodic and anodic species)
      for (int kk = 0; kk < numscal_; ++kk)
      {
        if ((conint[kk]) < 1e-13)  // 1.0E-16)
        {
          pow_conint_p *= std::pow(1e-13, p[kk]);
          pow_conint_q *= std::pow(1e-13, q[kk]);
#ifdef FOUR_C_ENABLE_ASSERTIONS
          std::cout << "WARNING: Rel. Conc. of species" << k
                    << " in Butler-Volmer formula is zero/negative: " << (conint[k]) << std::endl;
          std::cout << "-> Replacement value: pow(1.0E-16,p[ispec]) = " << pow(1e-13, p[k])
                    << " pow(1.0E-13,q[k]) = " << pow(1e-13, q[k]) << std::endl;
#endif
        }
        else
        {
          pow_conint_p *= std::pow((conint[kk]), p[kk]);
          pow_conint_q *= std::pow((conint[kk]), q[kk]);
        }
      }

      // derivation of concentration term  with respect to ionic species kk
      for (int kk = 0; kk < numscal_; ++kk)
      {
        pow_conint_p_derivative[kk] = pow_conint_p * p[kk] / conint[kk];
        pow_conint_q_derivative[kk] = pow_conint_q * q[kk] / conint[kk];
      }

      // loop over reacting species; determines the line of the matrix
      const double expterma = exp((1 - beta) * nume * frt * epd);
      const double exptermc = exp((-beta) * nume * frt * epd);

      for (int vi = 0; vi < nen_; ++vi)
      {
        // see Wittmann, Erweiterte Reaktionsmodelle fuer die numerische Simulation von
        // elektrochemischen Systemen, p.20, equ. 3.4
        const double fac_fns_funct_vi = scalar * faraday * nume * fac * fns * funct(vi);
        for (int ui = 0; ui < nen_; ++ui)
        {
          // loop over the columns of the matrix, makes sure that the linearisation w.r.t the first
          // concentration is added to the first column
          for (int kk = 0; kk < numscal_; ++kk)
          {
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + kk) +=
                -fac_fns_funct_vi *
                ((k_a * expterma * pow_conint_p_derivative[kk]) -
                    (k_c * exptermc * pow_conint_q_derivative[kk])) *
                funct(ui) * timefac;
          }
          // linearisation w.r.t the potential
          emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
              -fac_fns_funct_vi *
              ((-k_a * (1 - beta) * nume * frt * expterma * pow_conint_p) -
                  (k_c * beta * nume * frt * exptermc * pow_conint_q)) *
              funct(ui) * timefac;
        }
        // ------------right-hand-side
        erhs[vi * numdofpernode_ + k] -=
            -(fac_fns_funct_vi *
                ((k_a * expterma * pow_conint_p) - (k_c * exptermc * pow_conint_q))) *
            timefac;
      }

      break;
    }

    case Inpar::ElCh::butler_volmer_bard:
    {
      // "Electrochemical Methods Fundamentals and Applications"
      // Bard and Faulkner, 2001, pp. 94 ff; pp. 99 eq. 3.4.10
      // reaction model for a one-step, one-electron process (elementar step)
      // O + e -> R
      const auto e0 = cond.parameters().get<double>("E0");
      const auto k0 = cond.parameters().get<double>("K0");
      const auto beta = cond.parameters().get<double>("BETA");
      const auto c_c0 = cond.parameters().get<double>("C_C0");
      const auto c_a0 = cond.parameters().get<double>("C_A0");
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      if (dlcap != 0.0)
        FOUR_C_THROW(
            "double layer charging is not implemented for Butler-Volmer-Bard electrode kinetics");

      if (nume != 1)
        FOUR_C_THROW(
            "electron != 1; \n "
            "this Butler-Volmer-equation (Bard/Faulkner) works for elementary steps (one electron) "
            "only!");

      // only one reactant and product are supported by the basic model
      // only stoichiometry of 1
      {
        int check1 = 0;
        for (int kk = 0; kk < numscal_; kk++)
        {
          if (abs(stoich[kk]) > 1)
            FOUR_C_THROW(
                "Stoichiometry is larger than 1!! \n"
                "This is not supported by the reaction model based on Bard");

          check1 += abs(stoich[kk]);
        }
        if (check1 > 2 or check1 == 0)
          FOUR_C_THROW(
              "More than one reactant or product defined!! \n"
              "This is not supported by the reaction model based on Bard");

        // At the moment it is not checked if two products or reactants are defined
      }

      // equilibrium potential (equilpot):
      // defined in Bard, 2001, p.98, eq. 3.4.3
      const double equilpot = e0 + (log(c_c0 / c_a0)) / (frt * nume);
      // overpotential based on equilibrium potential
      const double eta_equilpot = epd - equilpot;

      // negative sign: we look at electron flow
      const double i0 = k0 * pow(c_c0, 1 - beta) * pow(c_a0, beta) * nume * faraday;

      // reactant or product not a species in the electrolyte
      // -> concentration = 1.0
      double conctermc = 1.0;
      double concterma = 1.0;
      double conctermc_der = 1.0;
      double concterma_der = 1.0;
      // species id of the anodic and cathodic reactant
      int anodic = 0;
      int cathodic = 0;
      bool checkc = 0;
      bool checka = 0;

      // concentration terms for anodic and cathodic reaction
      // only one reactant and product are supported by the basic model
      // only stoichiometry of 1
      for (int kk = 0; kk < numscal_; kk++)
      {
        if (stoich[kk] == 1)
        {
          concterma = conint[kk] / c_a0;
          concterma_der = 1.0 / c_a0;
          anodic = kk;
          checka = true;
        }
        else if (stoich[kk] == -1)
        {
          conctermc = conint[kk] / c_c0;
          conctermc_der = 1.0 / c_c0;
          cathodic = kk;
          checkc = true;
        }
      }

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // some safety checks/ user warnings
      if (((1 - beta) * (frt * nume) * eta_equilpot) > 100.0)
        std::cout << "WARNING: Exp((1-beta)...) in Butler-Volmer law is near overflow!"
                  << exp((1 - beta) * (frt * nume) * eta_equilpot) << std::endl;
      if (((-beta) * (frt * nume) * eta_equilpot) > 100.0)
        std::cout << "WARNING: Exp(-beta...) in Butler-Volmer law is near overflow!"
                  << exp((-beta) * (frt * nume) * eta_equilpot) << std::endl;
#endif

      const double expterma = exp((1 - beta) * (frt * nume) * eta_equilpot);
      const double exptermc = exp((-beta) * (frt * nume) * eta_equilpot);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const double fac_i0_funct_vi = scalar * fac * fns * i0 * funct(vi);
        // ---------------------matrix
        for (int ui = 0; ui < nen_; ++ui)
        {
          // derivation wrt concentration
          if (checkc == true)
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + cathodic) +=
                fac_i0_funct_vi * conctermc_der * exptermc * funct(ui) * timefac;
          if (checka == true)
            emat(vi * numdofpernode_ + k, ui * numdofpernode_ + anodic) +=
                -fac_i0_funct_vi * concterma_der * expterma * funct(ui) * timefac;
          // derivation wrt potential
          emat(vi * numdofpernode_ + k, ui * numdofpernode_ + numscal_) +=
              -fac_i0_funct_vi *
              (-concterma * (1 - beta) * frt * nume * expterma -
                  conctermc * beta * nume * frt * exptermc) *
              funct(ui) * timefac;
        }
        // ------------right-hand-side
        erhs[vi * numdofpernode_ + k] -=
            -fac_i0_funct_vi * (concterma * expterma - conctermc * exptermc) * timefac;
      }

      break;
    }

    case Inpar::ElCh::nernst:
      break;

    default:
    {
      FOUR_C_THROW("Kinetic model not implemented!");
      break;
    }
  }

  return;
}  // Discret::Elements::ScaTraEleUtilsElch<distype>::evaluate_elch_kinetics_at_integration_point


/*----------------------------------------------------------------------------------------------------------------*
 | evaluate electrode kinetics status information at integration point on domain or boundary element
 fang 07/15 |
 *----------------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::ScaTraEleUtilsElch<distype>::evaluate_electrode_status_at_integration_point(
    const Core::Elements::Element* ele,        ///< current element
    Core::LinAlg::SerialDenseVector& scalars,  ///< scalars to be computed
    const Teuchos::ParameterList& params,      ///< parameter list
    Core::Conditions::Condition& cond,         ///< condition
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephidtnp,  ///< nodal time derivative vector
    const Core::LinAlg::Matrix<nen_, 1>& funct,  ///< shape functions at integration point
    const int zerocur,                           ///< flag for zero current
    const int kinetics,                          ///< desired electrode kinetics model
    const std::vector<int>& stoich,              ///< stoichiometry of the reaction
    const int nume,                              ///< number of transferred electrons
    const double pot0,     ///< actual electrode potential on metal side at t_{n+1}
    const double frt,      ///< factor F/RT
    const double timefac,  ///< factor due to time discretization
    const double fac,      ///< integration factor
    const double scalar,   ///< scaling factor for current related quantities
    const int k            ///< index of evaluated scalar
) const
{
  // get Faraday constant
  const double faraday = Discret::Elements::ScaTraEleParameterElch::instance("scatra")->faraday();

  // get variables with their current values
  // current integrals: (i = epsilon i^E ) is calculated in case of porous media
  double currentintegral(0.);
  double currentdlintegral(0.);
  double boundaryint(0.);
  double electpotentialint(0.);
  double overpotentialint(0.);
  double electdiffpotint(0.);
  double opencircuitpotint(0.);
  double concentrationint(0.);
  double currderiv(0.);
  double currentresidual(0.);
  double boundaryint_porous(0.);

  // concentration of active species at integration point
  std::vector<double> conint(numscal_, 0.0);

  // elch-specific values at integration point
  for (int kk = 0; kk < numscal_; ++kk) conint[kk] = funct.dot(ephinp[kk]);

  // el. potential at integration point
  const double potint = funct.dot(ephinp[numscal_]);

  // history term of el. potential at integration point
  const double potdtnpint = funct.dot(ephidtnp[numscal_]);

  // electrode potential difference epd at integration point
  double epd = (pot0 - potint);

  // linearization of current w.r.t applied electrode potential "pot0"
  double linea(0.0);

  // concentration-dependent Butler-Volmer law(s)
  switch (kinetics)
  {
    case Inpar::ElCh::butler_volmer:
    case Inpar::ElCh::butler_volmer_yang1997:
    {
      // read model-specific parameter
      const auto alphaa = cond.parameters().get<double>("ALPHA_A");
      const auto alphac = cond.parameters().get<double>("ALPHA_C");
      double i0 = cond.parameters().get<double>("I0");
      if (i0 < -1e-14)
        FOUR_C_THROW(
            "i0 is negative, \n"
            "a positive definition is necessary due to generalized reaction models: {}",
            i0);
      const auto gamma = cond.parameters().get<double>("GAMMA");
      const auto refcon = cond.parameters().get<double>("REFCON");
      if (refcon < 1e-12) FOUR_C_THROW("reference concentration is too small: {}", refcon);

      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      double pot0dtnp = 0.0;
      double pot0hist = 0.0;
      if (dlcap != 0.0)
      {
        pot0dtnp = cond.parameters().get<double>("pot0dtnp");
        pot0hist = cond.parameters().get<double>("pot0hist");
      }

      // opencircuit potential is assumed to be zero here
      double ocp = 0.0;
      // surface overpotential based on opencircuit potential
      double eta = 0.0;
      // electrode potential
      double elepot = 0.0;

      if (zerocur == 0)
      {
        eta = epd - ocp;
        elepot = pot0;
      }
      else if (zerocur == 1)
      {
        elepot = potint + ocp;
        epd = ocp;
      }
      else
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

      double expterm(0.0);
      if (kinetics == Inpar::ElCh::butler_volmer)
      {
        // general Butler-Volmer
        expterm = std::pow(conint[k] / refcon, gamma) *
                  (exp(alphaa * frt * eta) - exp((-alphac) * frt * eta));
        linea = std::pow(conint[k] / refcon, gamma) * frt *
                ((alphaa * exp(alphaa * frt * eta)) + (alphac * exp((-alphac) * frt * eta)));
      }
      if (kinetics == Inpar::ElCh::butler_volmer_yang1997)
      {
        if (((conint[k] / refcon) < 1e-13) && (gamma < 1.0))
        {  // prevents NaN's in the current density evaluation
          expterm =
              (exp(alphaa * frt * eta) - (pow(1e-13 / refcon, gamma) * exp((-alphac) * frt * eta)));
          linea = ((alphaa)*frt * exp(alphaa * frt * eta)) +
                  (pow(1e-13 / refcon, gamma) * alphac * frt * exp((-alphac) * frt * eta));
        }
        else
        {
          expterm = (exp(alphaa * frt * eta) -
                     (pow(conint[k] / refcon, gamma) * exp((-alphac) * frt * eta)));
          linea = ((alphaa)*frt * exp(alphaa * frt * eta)) +
                  (pow(conint[k] / refcon, gamma) * alphac * frt * exp((-alphac) * frt * eta));
        }
      }

      // scan for NaNs due to negative concentrations under exponent gamma
      if (std::isnan(expterm) or std::isnan(linea))
        FOUR_C_THROW("NaN detected in electrode status calculation");

      // compute integrals
      electpotentialint += elepot * fac;
      overpotentialint += eta * fac;
      electdiffpotint += epd * fac;
      opencircuitpotint += ocp * fac;
      currentintegral += scalar * i0 * expterm * fac;  // the negative(!) normal flux density
      boundaryint += fac;
      boundaryint_porous += fac * scalar;
      concentrationint += scalar * conint[k] * fac;

      // tangent and rhs (= negative residual) for galvanostatic equation
      currderiv += scalar * i0 * linea * timefac * fac;
      currentresidual += scalar * i0 * expterm * timefac * fac;

      if (dlcap != 0.0)
      {
        currentdlintegral += fac * dlcap * (pot0dtnp - potdtnpint);

        // add contributions due to double-layer capacitance
        // positive due to redefinition of the exchange current density
        currderiv += scalar * fac * dlcap;
        currentresidual += scalar * fac * dlcap * (pot0 - pot0hist - (timefac * potdtnpint));
      }
      break;
    }

    // Tafel law:
    // implementation of cathodic path: i_n = i_0 * (-exp(-alpha * frt* eta)
    // -> cathodic reaction path: i_0 > 0 and alpha > 0
    // -> anodic reacton path:    i_0 < 0 and alpha < 0
    case Inpar::ElCh::tafel:
    {
      // read model-specific parameter
      const auto alpha = cond.parameters().get<double>("ALPHA");
      auto i0 = cond.parameters().get<double>("I0");

      const auto gamma = cond.parameters().get<double>("GAMMA");
      const auto refcon = cond.parameters().get<double>("REFCON");
      if (refcon < 1e-12) FOUR_C_THROW("reference concentration is too small: {}", refcon);
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      if (dlcap != 0.0)
        FOUR_C_THROW("double layer charging is not implemented for Tafel electrode kinetics");

      // opencircuit potential is assumed to be zero here
      double ocp = 0.0;
      // surface overpotential based on opencircuit potential
      double eta = 0.0;
      // electrode potential
      double elepot = 0.0;

      if (zerocur == 0)
      {
        eta = epd - ocp;
        elepot = pot0;
      }
      else if (zerocur == 1)
      {
        elepot = potint + ocp;
        epd = ocp;
      }
      else
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

      const double expterm = std::pow(conint[k] / refcon, gamma) * (-exp((-alpha) * frt * eta));
      linea = std::pow(conint[k] / refcon, gamma) * frt * (alpha * exp((-alpha) * frt * eta));
      // compute integrals
      electpotentialint += elepot * fac;
      overpotentialint += eta * fac;
      electdiffpotint += epd * fac;
      opencircuitpotint += ocp * fac;
      currentintegral += scalar * i0 * expterm * fac;  // the negative(!) normal flux density
      boundaryint += fac;
      boundaryint_porous += fac * scalar;
      concentrationint += conint[k] * fac;

      // tangent and rhs (= negative residual) for galvanostatic equation
      currderiv += scalar * i0 * linea * timefac * fac;
      currentresidual += scalar * i0 * expterm * timefac * fac;

      break;
    }

    // linear law:  i_n = frt*i_0*((alpha_a+alpha_c)*(V_M - phi)) -> changed 10/13
    // previously implemented: i_n = i_0*(alphaa*frt*(V_M - phi) + 1.0)
    //                         -> linearization in respect to anodic branch!!
    //                         this is not the classical version of a linear electrode kinetics law
    case Inpar::ElCh::linear:
    {
      // read model-specific parameter
      const auto alphaa = cond.parameters().get<double>("ALPHA");
      auto i0 = cond.parameters().get<double>("I0");
      if (i0 < -1e-14)
        FOUR_C_THROW(
            "i0 is negative, \n"
            "a positive definition is necessary due to generalized reaction models: {}",
            i0);
      const auto gamma = cond.parameters().get<double>("GAMMA");
      const auto refcon = cond.parameters().get<double>("REFCON");
      if (refcon < 1e-12) FOUR_C_THROW("reference concentration is too small: {}", refcon);
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      double pot0dtnp = 0.0;
      double pot0hist = 0.0;
      if (dlcap != 0.0)
      {
        pot0dtnp = cond.parameters().get<double>("pot0dtnp");
        pot0hist = cond.parameters().get<double>("pot0hist");
      }

      // opencircuit potential is assumed to be zero here
      double ocp = 0.0;
      // surface overpotential based on opencircuit potential
      double eta = 0.0;
      // electrode potential
      double elepot = 0.0;

      if (zerocur == 0)
      {
        eta = epd - ocp;
        elepot = pot0;
      }
      else if (zerocur == 1)
      {
        elepot = potint + ocp;
        epd = ocp;
      }
      else
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

      // compute integrals
      electpotentialint += elepot * fac;
      overpotentialint += eta * fac;
      electdiffpotint += epd * fac;
      opencircuitpotint += ocp * fac;
      currentintegral += scalar * i0 * pow(conint[k] / refcon, gamma) * (alphaa * frt * eta) *
                         fac;  // the negative(!) normal flux density
      boundaryint += fac;
      boundaryint_porous += fac * scalar;
      concentrationint += conint[k] * fac;

      // tangent and rhs (= negative residual) for galvanostatic equation
      linea = std::pow(conint[k] / refcon, gamma) * (alphaa * frt);
      currderiv += scalar * i0 * linea * timefac * fac;
      currentresidual +=
          scalar * i0 * pow(conint[k] / refcon, gamma) * (alphaa * frt * eta) * timefac * fac;

      if (dlcap != 0.0)
      {
        currentdlintegral += fac * dlcap * (pot0dtnp - potdtnpint);

        // add contributions due to double-layer capacitance
        // positive due to redefinition of the exchange current density
        currderiv += fac * dlcap;
        currentresidual += fac * dlcap * (pot0 - pot0hist - (timefac * potdtnpint));
      }

      break;
    }

    case Inpar::ElCh::butler_volmer_newman:
    {
      // "Electrochemical systems"
      // Newman ad Thomas-Alyea, 2004
      // General stoichiometry: pp. 212-213, e.q. 8.26
      // consideration of a elementary step of the form:
      // Sum_i s_i M_i ->  ne-
      // n is one if charge transfer is involved, multiple electron transfers "being unlikely in
      // an elementary step

      const auto k_a = cond.parameters().get<double>("K_A");
      const auto k_c = cond.parameters().get<double>("K_C");
      const auto beta = cond.parameters().get<double>("BETA");
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      if (dlcap != 0.0)
        FOUR_C_THROW(
            "double layer charging is not implemented for Butler-Volmer-Newman electrode kinetics");
      if (zerocur != 0)
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur is not implemented for this specific kinetic "
            "model.");

      // reaction order of the cathodic and anodic reactants of ionic species k
      std::vector<int> q(numscal_, 0);
      std::vector<int> p(numscal_, 0);

      for (int kk = 0; kk < numscal_; kk++)
      {
        // according to the convention: anodic reactant is positive
        if (stoich[kk] > 0)
        {
          q[kk] = 0;
          p[kk] = stoich[kk];
        }
        // according to the convention: cathodic reactant is negative
        else
        {
          q[kk] = -stoich[kk];
          p[kk] = 0;
        }
      }

      // linearization of current w.r.t applied electrode potential "pot0"
      double linea(0.0);
      double expterma(0.0);
      double exptermc(0.0);
      double expterm(0.0);
      double pow_conint_p = 1.0;  // product over i (c_i)^(p_i)
      double pow_conint_q = 1.0;  // product over i (c_i)^(q_i)

      // overpotential based on opencircuit potential
      double eta = 0.0;
      // electrode potential
      double elepot = 0.0;

      // open circuit potential (ocp): time dependent electrode surface concentrations
      // defined in Newman, 2004, p. 211, eq. 8.20
      double ocp = 1 / frt / nume * log(k_c / k_a);
      for (int kk = 0; kk < numscal_; ++kk)
      {
        ocp += 1 / frt / nume * (q[kk] - p[kk]) * log(conint[kk]);
        // safety check
        if ((q[kk] - p[kk]) != -stoich[kk])
          FOUR_C_THROW("stoichiometry factors and the factors q,p do not correlate!!");
      }

      if (zerocur == 0)
      {
        // overpotential based on open circuit potential
        eta = epd - ocp;
        elepot = pot0;
      }
      else if (zerocur == 1)
      {
        elepot = potint + ocp;
        epd = ocp;
      }
      else
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

      for (int kk = 0; kk < numscal_; ++kk)
      {
        if ((conint[kk]) < 1e-13)
        {
          pow_conint_p *= std::pow(1e-13, p[kk]);
          pow_conint_q *= std::pow(1e-13, q[kk]);
#ifdef FOUR_C_ENABLE_ASSERTIONS
          std::cout << "WARNING: Rel. Conc. of species" << kk
                    << " in Butler-Volmer formula is zero/negative: " << (conint[kk]) << std::endl;
          std::cout << "-> Replacement value: pow(EPS,p[ispec]) = " << pow(1e-13, p[kk])
                    << " pow(1.0E-16,q[i]) = " << pow(1e-13, q[kk]) << std::endl;
#endif
        }
        else
        {
          pow_conint_p *= std::pow((conint[kk]), p[kk]);
          pow_conint_q *= std::pow((conint[kk]), q[kk]);
        }
      }
      expterma = exp((1 - beta) * nume * frt * epd);
      exptermc = exp(-beta * nume * frt * epd);
      linea = nume * faraday *
              (frt * nume *
                  ((k_a * (1 - beta) * expterma * pow_conint_p) -
                      (k_c * (-1) * beta * exptermc * pow_conint_q)));

      // scan for NaNs due to negative concentrations under exponent gamma
      if (std::isnan(expterm) or std::isnan(linea))
        FOUR_C_THROW("NaN detected in electrode status calculation");

      // compute integrals
      currentintegral += scalar * nume * faraday *
                         ((k_a * expterma * pow_conint_p) - (k_c * exptermc * pow_conint_q)) * fac;
      boundaryint += fac;
      boundaryint_porous += fac * scalar;
      electpotentialint += elepot * fac;
      overpotentialint += eta * fac;
      electdiffpotint += epd * fac;
      opencircuitpotint += ocp * fac;
      concentrationint += conint[k] * fac;

      // tangent and rhs (= negative residual) for galvanostatic equation
      currderiv += scalar * linea * fac * timefac;
      currentresidual += scalar * nume * faraday *
                         ((k_a * expterma * pow_conint_p) - (k_c * exptermc * pow_conint_q)) *
                         timefac * fac;

      break;
    }

    case Inpar::ElCh::butler_volmer_bard:
    {
      // "Electrochemical Methods Fundamentals and Applications"
      // Bard and Faulkner, 2001, pp. 94 ff; pp. 99 eq. 3.4.10
      // reaction model for a one-step, one-electron process
      // O + e -> R
      const auto e0 = cond.parameters().get<double>("E0");
      const auto k0 = cond.parameters().get<double>("K0");
      const auto beta = cond.parameters().get<double>("BETA");
      const auto c_c0 = cond.parameters().get<double>("C_C0");
      const auto c_a0 = cond.parameters().get<double>("C_A0");
      const auto dlcap = cond.parameters().get<double>("DL_SPEC_CAP");
      if (dlcap != 0.0)
        FOUR_C_THROW(
            "double layer charging is not implemented for Butler-Volmer-Bard electrode kinetics");
      if (zerocur != 0)
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur is not implemented for this specific kinetic "
            "model.");

      if (nume != 1)
        FOUR_C_THROW(
            "electron != 1; \n "
            "this Butler-Volmer-equation (Bard/Faulkner) works for elementary steps (one electron) "
            "only!");

      // only one reactant and product are supported by the basic model
      // only stoichiometry of 1
      {
        int check1 = 0;
        for (int kk = 0; kk < numscal_; kk++)
        {
          if (abs(stoich[kk]) > 1)
            FOUR_C_THROW(
                "Stoichiometry is larger than 1!! \n"
                "This is not supported by the reaction model based on Bard");

          check1 += abs(stoich[kk]);
        }
        if (check1 > 2 or check1 == 0)
          FOUR_C_THROW(
              "More than one reactant or product defined!! \n"
              "This is not supported by the reaction model based on Bard");

        // In the moment it is not checked if two products (and no reactants) and vis versa are
        // defined
      }

      // reactant or product not a species in the electrolyte
      // -> concentration = 1.0
      double conctermc = 1.0;
      double concterma = 1.0;

      // concentration terms for anodic and cathodic reaction
      // only one reactant and product are supported by the basic model
      // only stoichiometry of 1
      for (int kk = 0; kk < numscal_; kk++)
      {
        if (stoich[kk] == 1)
          concterma = conint[kk] / c_a0;
        else if (stoich[kk] == -1)
          conctermc = conint[kk] / c_c0;
      }

      // equilibrium potential (equilpot):
      // defined in Bard, 2001, p.98, eq. 3.4.3
      const double equilpot = e0 + (log(c_c0 / c_a0)) / (frt * nume);
      // overpotential based on equilibrium potential
      const double eta_equilpot = epd - equilpot;
      // difference between equilibrium potential and open circuit potential:
      // -> equilpot: depends on initial electrode surface concentration
      // -> ocp:      depends on actual electrode surface concentration

      // open circuit potential (ocp): time dependent electrode surface concentrations
      // defined in Newman, 2004, p. 211, eq. 8.20
      const double ocp = e0 + 1 / frt / nume * log(conctermc / concterma);
      // overpotential based on open circuit potential
      const double eta = epd - ocp;

      const double expterma = exp((1 - beta) * (frt * nume) * eta_equilpot);
      const double exptermc = exp(-beta * (frt * nume) * eta_equilpot);
      const double linea = concterma * (1 - beta) * (frt * nume) * expterma +
                           conctermc * beta * (frt * nume) * exptermc;

      // scan for NaNs due to negative concentrations under exponent gamma
      if (std::isnan(expterma) or std::isnan(exptermc) or std::isnan(linea))
        FOUR_C_THROW("NaN detected in electrode status calculation");

      const double i0 = faraday * k0 * pow(c_c0, 1 - beta) * pow(c_a0, beta);

      // compute integrals
      overpotentialint += eta * fac;
      electdiffpotint += epd * fac;
      opencircuitpotint += ocp * fac;
      currentintegral += scalar * i0 * (concterma * expterma - conctermc * exptermc) *
                         fac;  // the negative(!) normal flux density
      boundaryint += fac;
      boundaryint_porous += fac * scalar;
      concentrationint += conint[k] * fac;  // concentration-output for the first species only

      // tangent and rhs (= negative residual) for galvanostatic equation
      currderiv += scalar * i0 * linea * timefac * fac;
      currentresidual +=
          scalar * i0 * (concterma * expterma - conctermc * exptermc) * timefac * fac;

      break;
    }  // end Butler-Volmer-Bard

    case Inpar::ElCh::nernst:
    {
      const auto e0 = cond.parameters().get<double>("E0");
      const auto c0 = cond.parameters().get<double>("C0");
      if (zerocur != 0)
        FOUR_C_THROW(
            "The electrode kinetics flag zero_cur is not implemented for this specific kinetic "
            "model.");

      // compute integrals
      overpotentialint += potint * fac;
      boundaryint += fac;
      boundaryint_porous += fac * scalar;
      concentrationint += conint[k] * fac;

      opencircuitpotint += e0 + (log(concentrationint / boundaryint / c0)) / (frt * nume);
      opencircuitpotint *= boundaryint;
      break;
    }

    default:
    {
      FOUR_C_THROW("Kinetic model not implemented");
      break;
    }
  }

  // add contributions from current integration point into result vector
  scalars(0) += currentintegral;
  scalars(1) += currentdlintegral;
  scalars(2) += boundaryint;
  scalars(3) += electpotentialint;
  scalars(4) += overpotentialint;
  scalars(5) += electdiffpotint;
  scalars(6) += opencircuitpotint;
  scalars(7) += concentrationint;
  scalars(8) += currderiv;
  scalars(9) += currentresidual;
  scalars(10) += boundaryint_porous;
}  // Discret::Elements::ScaTraEleUtilsElch<distype>::evaluate_electrode_status_at_integration_point


/*----------------------------------------------------------------------*
 | evaluate ion material                                     fang 07/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::ScaTraEleUtilsElch<distype>::mat_ion(
    const std::shared_ptr<const Core::Mat::Material> material,  //!< ion material
    const int k,                                                //!< ID of ion material
    const Inpar::ElCh::EquPot equpot,  //!< type of closing equation for electric potential
    const std::shared_ptr<ScaTraEleDiffManagerElch>& diffmanager  //!< diffusion manager
)
{
  // cast material to ion material
  const std::shared_ptr<const Mat::Ion> mation = std::static_pointer_cast<const Mat::Ion>(material);

  // valence of ionic species
  diffmanager->set_valence(mation->valence(), k);

  // concentration depending diffusion coefficient
  diffmanager->set_isotropic_diff(mation->diffusivity(), k);

  // Loop over materials is finished - now all material parameter are set
  if (k == numscal_ - 1)
  {
    // Material data of eliminated ion species is read from the LAST ion material
    // in the matlist!
    if (equpot == Inpar::ElCh::equpot_enc_pde_elim)
    {
      diffmanager->increase_length_vector(k, numscal_);

      // valence of ionic species
      diffmanager->set_valence(mation->elim_valence(), numscal_);

      // concentration depending diffusion coefficient
      diffmanager->set_isotropic_diff(mation->elim_diffusivity(), numscal_);
    }
  }
}


// template classes
// 1D elements
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::line2>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::line3>;

// 2D elements
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::quad4>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::quad8>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::quad9>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::tri3>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::tri6>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::nurbs3>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::hex8>;
// template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::hex20>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::hex27>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::tet4>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::tet10>;
// template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::wedge6>;
template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::pyramid5>;
// template class Discret::Elements::ScaTraEleUtilsElch<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
