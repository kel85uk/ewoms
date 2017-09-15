// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::FluidSystems::TNOBRINE_CH4
 */
#ifndef OPM_TNO_BRINE_CH4_FLUID_SYSTEM_HPP
#define OPM_TNO_BRINE_CH4_FLUID_SYSTEM_HPP

#include <opm/material/IdealGas.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>

#include "BrineVarSalinity.hpp"
#include "NaCl.hpp"
#include "CH4.hpp"
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/TabulatedComponent.hpp>
#include "brine_ch4.hpp"
#include <opm/common/Valgrind.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <iostream>
#include <cassert>

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with water and nitrogen as components.
 */
template <class Scalar, bool useComplexRelations = true>
class TNOBRINE_CH4
    : public BaseFluidSystem<Scalar, TNOBRINE_CH4<Scalar, useComplexRelations> >
{
    typedef TNOBRINE_CH4<Scalar, useComplexRelations> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    // convenience typedefs
    typedef Opm::IdealGas<Scalar> IdealGas;
    typedef Opm::H2O<Scalar> IapwsH2O;
    typedef Opm::TabulatedComponent<Scalar, IapwsH2O > TabulatedH2O;
    typedef Opm::CH4<Scalar> SimpleCH4;
    typedef Opm::NaCl<Scalar> SimpleNACL;
    typedef Opm::BrineVarSalinity<Scalar, TabulatedH2O> Brine;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = NullParameterCache<Evaluation>;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the liquid phase
    static const int liquidPhaseIdx = 0;
    //! Index of the gas phase
    static const int gasPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {
            "liquid",
            "gas"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        return
            (phaseIdx == gasPhaseIdx)
            ? true
            :H2O::liquidIsCompressible();// the water component decides for the liquid phase...
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return
            (phaseIdx == gasPhaseIdx)
            ? H2O::gasIsIdeal() && CH4::gasIsIdeal() // let the components decide
            : false; // not a gas
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 3;

    //! The component index of water
    static const int H2OIdx = 0;
    //! The component index of molecular nitrogen
    static const int CH4Idx = 1;
    //! The component index of NaCl
    static const int NACLIdx = 2;
    

    //! The component for pure water
    typedef TabulatedH2O H2O;

    //! The component for pure nitrogen
    typedef SimpleCH4 CH4;
    
    //! The component for NACL
    typedef SimpleNACL NACL;
    
    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        static const char* name[] = {
            H2O::name(),
            CH4::name(),
            NACL::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == CH4Idx)
            ? CH4::molarMass()
            : (compIdx == NACLIdx)
            ? NACL::molarMass()
            : 1e30;
    }


    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/50,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/50);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress, Scalar salinity)
    {
        if (H2O::isTabulated) {
            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
            Brine::salinity = salinity;
        }
    }

    /*!
     * \copydoc BaseFluidSystem::density
     *
     * If useComplexRelations == true, we apply Formula (2.6) from S.O.Ochs: "Development
     * of a multiphase multicomponent model for PEMFC - Technical report: IRTG-NUPUS",
     * University of Stuttgart, 2008
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        typedef Opm::MathToolbox<LhsEval> LhsToolbox;
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& T = FsToolbox::template decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = FsToolbox::template decay<LhsEval>(fluidState.pressure(phaseIdx));
        

        switch (phaseIdx)
        {
            case liquidPhaseIdx:
                {
                    const auto& xlNACL = FsToolbox::template decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, NACLIdx));
                    const auto& XlNACL = moleTomassFrac_(xlNACL);
                    return Brine::liquidDensity(T,p,XlNACL);
                }
            case gasPhaseIdx:
                {
                    const auto& xgH2O = FsToolbox::template decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx));
                    return gasDensity_(T, p, xgH2O);
                }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
        }
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        typedef Opm::MathToolbox<LhsEval> LhsToolbox;
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& T = FsToolbox::template decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = FsToolbox::template decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            {
                const auto& xNaCl = FsToolbox::template decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, NACLIdx));
                const auto& XNaCl = moleTomassFrac_(xNaCl);
                auto result = Brine::liquidViscosity(T, p, XNaCl);
                // assume pure water for the liquid phase
                return result;
            }

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        if (!useComplexRelations)
            // assume pure nitrogen for the gas phase
            DUNE_THROW(Dune::InvalidStateException, "Gas phase must use complex relations \n");
        else {
            /* Wilke method. See:
             *
             * See: R. Reid, et al.: The Properties of Gases and Liquids,
             * 4th edition, McGraw-Hill, 1987, 407-410
             * 5th edition, McGraw-Hill, 20001, p. 9.21/22
             */
            LhsEval muResult = 0;
            const auto& XgH2O = FsToolbox::template decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, H2OIdx));
			muResult = (1-XgH2O) * CH4::gasViscosity(T, p) + XgH2O * H2O::gasViscosity(T, p);
			Valgrind::CheckDefined(muResult);
            return muResult;
        }
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& T = FsToolbox::template decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& pg = FsToolbox::template decay<LhsEval>(fluidState.pressure(gasPhaseIdx)); 
        auto xlNaCl =  FsToolbox::template decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, NACLIdx));
        
        if (compIdx == NACLIdx)
        {
        	//Salts can only be present in the liquid phase
        	if(phaseIdx == liquidPhaseIdx){
                return 1e3/p;
            }
        	else
        		return 1e6;
        }

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            /*auto xgH2O = H2O::vaporPressure(T)/(pg);
            auto xgCH4 = 1.0 - xgH2O;
            auto xlCH4 = xgCH4*pg/(Opm::BinaryCoeff::Brine_CH4<Scalar>::henry(T));
            Scalar eps = 1e-6;
            if(FsToolbox::isnan(xlCH4))  xlCH4 = eps;
            xlCH4 = FsToolbox::max(eps,FsToolbox::min(1.0-eps,xlCH4));
            xgH2O = FsToolbox::max(eps,FsToolbox::min(1.0-eps,xgH2O));
            xlNaCl = FsToolbox::max(eps,FsToolbox::min(1.0-eps,xlNaCl));
            const auto xlH2O = 1.0 - xlCH4 - xlNaCl;
            xgCH4 = 1.0 - xgH2O;
            if (compIdx == H2OIdx)
                return 1.0*xgH2O/(xlH2O);
            return 1.0*xgCH4/(xlCH4);*/
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == CH4Idx)
                return Opm::BinaryCoeff::Brine_CH4<Scalar>::henry(T)/p;
        }

        assert(phaseIdx == gasPhaseIdx);

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned /*compIdx*/)

    {
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& T = FsToolbox::template decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = FsToolbox::template decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return BinaryCoeff::Brine_CH4<Scalar>::liquidDiffCoeff(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);
        return BinaryCoeff::Brine_CH4<Scalar>::gasDiffCoeff(T, p);
    }
	
	template <class Evaluation=Scalar>
    static Evaluation salinityTomoleFrac(const Evaluation& sal)
    {
    	const Evaluation Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
    	const Evaluation Ms = NACL::molarMass(); /* molecular weight of NaCl  [kg/mol] */
		
		const auto& y = sal;
    	/* XlNaCl: conversion from mass fraction to mol fraction */
		const auto& xNaCl = y*(Mw * Ms)/(Ms * ((Mw - Ms)*y + Ms) );
    	return xNaCl;
    }


private:
    static Scalar liquidDensity_(Scalar T, Scalar pl, Scalar xlCH4, Scalar xlH2O, Scalar XlNaCl)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(XlNaCl);
        Valgrind::CheckDefined(xlCH4);

        if(T < 273.15 || T > 623.15) {
            DUNE_THROW(Dune::InvalidStateException,
                       "Liquid density for Brine and air is only defined between 273.15K and 623.15K (is " << T << ")");
        }
        if(pl >= 1.0e8) {
            DUNE_THROW(Dune::InvalidStateException,
                       "Liquid density for Brine and air is only defined below 100MPa (is " << pl << ")");
        }

        //Scalar rho_brine = Brine::liquidDensity(T, pl, XlNaCl); 	// WARNING: Here we just neglect the influence of dissolved CH4 in Brine
        Scalar rho_brine = Brine::liquidDensity(T, pl); 	// WARNING: Here we just neglect the influence of dissolved CH4 in Brine
        return rho_brine;
    }
    template <class Evaluation=Scalar>
    static Evaluation moleTomassFrac_(const Evaluation& xlNaCl)
    {
    	const Evaluation Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
    	const Evaluation Ms = NACL::molarMass(); /* molecular weight of NaCl  [kg/mol] */

    	const auto& x_NaCl = xlNaCl;
    	/* XlNaCl: conversion from mass fraction to mol fraction */
    	const auto& XlNaCl = Ms * x_NaCl / ((Ms - Mw) * x_NaCl + Mw);
    	return XlNaCl;
    }
    static Scalar liquidEnthalpyBrine_(Scalar T, Scalar p, Scalar XlNaCl)
    {
        /* XlCH4 : mass fraction of CH4 in brine */
        /* same function as enthalpy_brine, only extended by CH4 content */
        /*Numerical coefficents from PALLISER (The values for F [] are given in ADAMS and BACHO)*/

        static const Scalar f[] = {2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10};

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] =
        {{ 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }};

        Scalar theta, h_NaCl;
        Scalar m, h_ls, d_h, hw;
        Scalar S_lSAT, delta_h;
        int i, j;
        XlNaCl = std::abs(XlNaCl); // Added by Vishal

        theta = T - 273.15;

        S_lSAT = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta; // This is NaCl specific (S_lSAT = 0.265234)
       // std::cout<<"saturation limit : " << S_lSAT<< std::endl;

        /*Regularization*/
        if (XlNaCl > S_lSAT)
            XlNaCl = S_lSAT;

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T // enthalpy of Halite (Confirm with Dr. Melani)
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        m = (1E3/58.44)*(XlNaCl/(1-XlNaCl));
        i = 0;
        j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
            }
        }
        /* heat of dissolution for halite according to Michaelides 1971 */
        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without CH4 */
        h_ls = (1-XlNaCl)*hw + XlNaCl*h_NaCl + XlNaCl*delta_h; /* kJ/kg */
        return (h_ls);
    };
    
    template <class Evaluation>
    static Evaluation gasDensity_(const Evaluation& T, const Evaluation& pg, const Evaluation& xgH2O)
    {
    	auto pH2O = xgH2O*pg; //Dalton' Law
        auto pCH4 = pg - pH2O;
        auto gasDensityCH4 = CH4::gasDensity(T, pCH4);
        auto gasDensityH2O = H2O::gasDensity(T, pH2O);
        auto gasDensity = gasDensityCH4 + gasDensityH2O;
        return gasDensity;
    }
};

} // namespace FluidSystems

} // namespace Opm

#endif
