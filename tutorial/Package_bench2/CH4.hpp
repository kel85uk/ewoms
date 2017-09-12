// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \ingroup Components
 *
 * \brief Properties of methane \f$CH_4\f$.
 */
#ifndef OPM_CH4_HPP
#define OPM_CH4_HPP

#include <opm/material/IdealGas.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/components/Component.hpp>

#include <cmath>

namespace Opm
{

/*!
 * \ingroup Components
 *
 * \brief Properties of pure molecular methane \f$CH_4\f$.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class CH4 : public Component<Scalar, CH4<Scalar> >
{
    typedef Opm::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for methane.
     */
    static const char *name()
    { return "CH4"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of molecular methane.
     */
    static Scalar molarMass()
    { return 16.043e-3;}

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of molecular methane
     */
    static Scalar criticalTemperature()
    { return 190.4; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of molecular methane
     */
    static Scalar criticalPressure()
    { return 46e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at molecular methane's triple point.
     */
    static Scalar tripleTemperature()
    { return 90.7; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at molecular methane's triple point.
     */
    static Scalar triplePressure()
    { return 0; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure molecular methane
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     */
    //static Scalar vaporPressure(Scalar T)
    //{ DUNE_THROW(Dune::NotImplemented, "vaporPressure for CH4"); }
    template <class Evaluation = Scalar>
    static Evaluation vaporPressure(const Evaluation& temperature)
    {
        typedef MathToolbox<Evaluation> Toolbox;

        if (temperature > criticalTemperature())
            return criticalPressure();
        if (temperature < tripleTemperature())
            return 0; // N2 is solid: We don't take sublimation into
                      // account

        // note: this is the ancillary equation given on page 1368
        const Evaluation& sigma = 1.0 - temperature/criticalTemperature();
        const Evaluation& sqrtSigma = Toolbox::sqrt(sigma);
        const Scalar N1 = -6.12445284;
        const Scalar N2 = 1.26327220;
        const Scalar N3 = -0.765910082;
        const Scalar N4 = -1.77570564;
        return
            criticalPressure() *
            Toolbox::exp(criticalTemperature()/temperature*
                         (sigma*(N1 +
                                 sqrtSigma*N2 +
                                 sigma*(sqrtSigma*N3 +
                                        sigma*sigma*sigma*N4))));
    }
    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of \f$CH_4\f$ gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation = Scalar>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(Evaluation(molarMass()), temperature, pressure);
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of gaseous \f$CH_4\f$ in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation = Scalar>
    static Evaluation gasPressure(const Evaluation& temperature, const Evaluation& density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure methane gas.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 154, 657, 671
     */
    template <class Evaluation = Scalar>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& pressure)
    {
        // method of Joback
        const Evaluation& cpVapA = 19.25;
        const Evaluation& cpVapB = 0.05213;
        const Evaluation& cpVapC = 1.197e-5;
        const Evaluation& cpVapD = -1.132e-8;

        //Scalar cp =
        //    cpVapA + T*(cpVapB + T*(cpVapC + T*cpVapD));

        // calculate: \int_0^T c_p dT
        return
            1/molarMass()* // conversion from [J/mol] to [J/kg]

            temperature*(cpVapA + temperature*
               (cpVapB/2 + temperature*
                (cpVapC/3 + temperature*
                 (cpVapD/4))));
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure methane gas.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation = Scalar>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {

        return
            gasEnthalpy(temperature, pressure) -
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]
            IdealGas::R*temperature; // = pressure * spec. volume for an ideal gas
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$CH_4\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids,
     * 4th edition, McGraw-Hill, 1987, pp 396-397, 670
     * 5th edition, McGraw-Hill, 2001  pp 9.7-9.8 (omega and V_c taken from p. A.5)
     *
     */
    template <class Evaluation = Scalar>
    static Evaluation gasViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        typedef MathToolbox<Evaluation> Toolbox;
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 98.6; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.011; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        const Scalar dipole = 0.0; // dipole moment [debye]

        Scalar mu_r4 = 131.3 * dipole / std::sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;

        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        const Evaluation& Tstar = 1.2593 * temperature/Tc;
        const Evaluation& Omega_v =
            1.16145*Toolbox::pow(Tstar, -0.14874) +
            0.52487*Toolbox::exp(- 0.77320*Tstar) +
            2.16178*Toolbox::exp(- 2.43787*Tstar);
        const Evaluation& mu = 40.785*Fc*Toolbox::sqrt(M*temperature)/(Toolbox::pow(Vc, 2./3)*Omega_v);

        // convertion from micro poise to Pa s
        return mu/1e6 / 10;
    }
    
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& temperature,
                                             const Evaluation& pressure)
    { return 0.035; }// Value liberally taken from Engineering toolbox... This should never be used anyways!!!!!! - KKL REF:http://www.engineeringtoolbox.com/methane-d_1420.html
};

} // end namespace

#endif

