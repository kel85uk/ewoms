/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Binary coefficients for water and methane.
 */
#ifndef OPM_BINARY_COEFF_BRINE_CH4_HH
#define OPM_BINARY_COEFF_BRINE_CH4_HH

#include <opm/material/binarycoefficients/HenryIapws.hpp>
#include <opm/material/binarycoefficients/FullerMethod.hpp>

#include "CH4.hpp"
#include <opm/material/components/H2O.hpp>

namespace Opm
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and methane.
 */
template<class Scalar>
class Brine_CH4
{
    typedef Opm::H2O<Scalar> H2O;
    typedef Opm::CH4<Scalar> CH4;
public:


    static constexpr Scalar Ms_ = 58.8e-3; //Molar Mass NaCl
    static constexpr Scalar eps_ = 1.0E-15;



    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for molecular methane in liquid water.
     *
     * See:
     *
     * IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid
     * Distribution Constant for Gases in H2O and D2O at High
     * Temperatures"
     * http://www.iapws.org/relguide/HenGuide.pdf
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& temperature)
    {
        const Scalar E = 2215.6977;
        const Scalar F = -0.1089;
        const Scalar G = -6.6240;
        const Scalar H = 4.6789;
        
        return henryIAPWS(E, F, G, H, temperature);
    };

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and methane.
     *
     * \copybody fullerMethod()
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        typedef Opm::H2O<Scalar> H2O;
        typedef Opm::CH4<Scalar> CH4;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  25.14 /* CH4, estimated from the increments */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*1e3, CH4::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    };

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular methane in liquid water.
     *
     * The empirical equations for estimating the diffusion
     * coefficient in infinite solution which are presented in Reid,
     * 1987 all show a linear dependency on temperature. We thus
     * simply scale the experimentally obtained diffusion coefficient
     * given in Kobayashi by the temperature.
     *
     * See:
     *
     * R. Reid et al.: "The properties of Gases and Liquids", 4th edition,
     * pp. 599, McGraw-Hill, 1987
     *
     * K. Kobayashi: "Optimization Methods for Multiphase Systems in
     * the Subsurface - Application to Methane Migration in Coal
     * Mining Arenas", PhD thesis, University of Stuttgart, Institute
     * of Hydraulic Engineering (Mitteilungsheft 139), p 58, 2004
     */
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& Texp = 273.15 + 25; // [K]
        const Evaluation& Dexp = 3.55e-9; // [m^2/s]
        return Dexp * temperature/Texp;
    };

    //Calculation of the CH4 solubility in brine according to Duan et al. The prediction of methane solubility in natural waters to high ionic strength
    //from 0 to 250°C and from 0 to 1600 bar (1992)
    template <class Evaluation = Scalar>
    static void calculateMoleFraction(const Evaluation& Temp,
            const Evaluation& pg, const Evaluation& XlNaCl,
            const Evaluation& ygCH4, Evaluation &xlCH4)
    {
        Evaluation MH2O = H2O::molarMass();
        Evaluation xlNaCl   = -MH2O*XlNaCl/((Ms_-MH2O)*XlNaCl - MH2O); /* salinity: conversion from mass fraction to mole fraction */
        Evaluation mNaCl = -55.56*xlNaCl/(xlNaCl-1);         /* salinity: conversion from mole fraction to molality */

        Evaluation pg_bar = pg/1.0E5;

        Evaluation A = mComputeA(Temp,pg);    /* mu_{CH4}^{l(0)}/RT */
        Evaluation B = mComputeB(Temp,pg);    /* lambda_{CH4-Na+} */
        Evaluation C = mComputeC(Temp,pg);    /* Xi_{CH4-Na+-Cl-} */
        Evaluation phiCH4 = mFugacityCoeffCH4(Temp,pg);

        Evaluation exponent = A - log(phiCH4) + 2*B*mNaCl + C*mNaCl*mNaCl;

        Evaluation mCH4 = ygCH4 * pg_bar / exp(exponent);  /* paper: equation (6) */

        xlCH4 = mCH4/(mCH4 + 55.56);              /* conversion: molality to mole fraction */
//        XlCH4 = xlCH4w*MCH4/(xlCH4*MCH4 + (1-xlCH4)*MH2O);   /* conversion: mole fraction to mass fraction */


    }
private:

    /**********************************************************/
    /*                                                        */
    /* mComputeA: computation of mu_{CH4}^{l(0)}/RT            */
    /*                                                        */
    /**********************************************************/
    template <class Evaluation = Scalar>
    static Evaluation mComputeA(Evaluation& Temp, Evaluation& pg)
    {
        Evaluation T, pg_bar, p;
        Evaluation c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
        Evaluation A;

        T      = Temp;
        pg_bar = pg/1.0E5;   /* conversion from Pa to bar */
        p      = pg_bar;

            c1 = 43.0210345;
            c2 = -0.0683277221;
            c3 = -5687.1873;
            c4 = 3.56636281e-05;
            c5 = -57.913379;
            c6 = 6.11616662E-03;
            c7 = -7.85528103E-04;
            c8 = -9.42540759E-02;
            c9 = 1.92132040E-02;
            c10 = -9.17186899E-06;

        A = c1+c2*T+c3/T+c4*T*T+c5/(680.0-T)+c6*p+c7*p*log(T)+c8*p/T+c9*p/(680.0-T)+c10*p*p/T;

        return(A);
    }

    /**********************************************************/
    /*                                                        */
    /* mComputeB: computation of lambda_{CH4-Na}               */
    /*                                                        */
    /**********************************************************/
    template <class Evaluation = Scalar>
    static Evaluation mComputeB(Evaluation Temp, Evaluation pg)
    {
        Evaluation c1,c2,c8,c10;
        Evaluation T,pg_bar,p;
        Evaluation B;

        c1 = 9.92230792E-02;
        c2 = 2.57906811E-05;
        c8 = 1.83451402E-02;
        c10= -8.07196716E-06;

        T      = Temp;
        pg_bar = pg/1.0E5;   /* conversion from Pa to bar */
        p      = pg_bar;

        B = c1+c2*T+c8*p/T+c10*p*p/T;

        return(B);
    }

    /**********************************************************/
    /*                                                        */
    /* mComputeB: computation of xi_{CH4-Na-Cl}               */
    /*                                                        */
    /**********************************************************/
    template <class Evaluation = Scalar>
    static Evaluation mComputeC(Evaluation Temp, Evaluation pg)
    {
        Evaluation c1;
        Evaluation C;

        c1 = -6.23943799E-03;
        C = c1;

        return(C);
    }

    /************************************************************************/
    /*                                                                      */
    /* Computation of real gas factor Z as suggested in Duan,Moeller and    */
    /* Weare (1992)                                                         */
    /*                                                                      */
    /************************************************************************/

    template <class Evaluation = Scalar>
    static Evaluation mComputeZ(Evaluation Temp, Evaluation pg)
    {
        Evaluation pr, pcrit, Tr, Tcrit, Vr, R, Z;

        pcrit = CH4::criticalPressure();
        Tcrit = CH4::criticalTemperature();
        R   = 8.314467;   /* universal gas constant [Pa m^3/(K mol)] */

        pr = pg/pcrit;    /* computation of reduced pressure */
        Tr = Temp/Tcrit;  /* computation of reduced temperature */

        Vr = m_iterateVr(Temp, pg);

        Z = pr*Vr/Tr;

    /*    printf("p = %.2f bar T = %.2f C \t Vr = %.4f \t Z = %.4f \n", pg/1.0E5, Temp-273.15, Vr, Z);
    */

        return Z;
    }
    template <class Evaluation = Scalar>
    static Evaluation mFugacityCoeffCH4(Evaluation Temp, Evaluation pg)
    {
        Evaluation A, B, C, D;
        Evaluation Tcrit, Tr, pcrit, pr, Vr;
        Evaluation a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
        Evaluation Z;
        Evaluation lnphiCH4, phiCH4;

        Z = 1.; //mComputeZ(Temp, pg);

            a1 = 0.08725553928;
            a2 = -0.752599476;
            a3 = 0.375419887;
            a4 = 0.0107291342;
            a5 = 0.0054962636;
            a6 = -0.0184772802;
            a7 = 0.000318993183;
            a8 = 0.000211079375;
            a9 = 2.01682801e-05;
            a10 = -1.65606189e-05;
            a11 = 0.000119614546;
            a12 = -0.000108087289;
            a13 = 0.0448262295;
            a14 = 0.75397;
            a15 = 0.077167;

        pcrit = CH4::criticalPressure();     /* critical pressure [Pa] */
        Tcrit = CH4::criticalTemperature();     /* critical temperature [K] */

        Tr = Temp/Tcrit;    /* reduced temperature */
        pr = pg/pcrit;    /* reduced pressure */

        Vr = Z*Tr/pr;    /* compute Vr backwards */

        A = a1 + a2/(Tr*Tr) + a3/(Tr*Tr*Tr);
        B = a4 + a5/(Tr*Tr) + a6/(Tr*Tr*Tr);
        C = a7 + a8/(Tr*Tr) + a9/(Tr*Tr*Tr);
        D = a10 + a11/(Tr*Tr) + a12/(Tr*Tr*Tr);

        lnphiCH4 =  Z - 1 - log(Z) + A/Vr + B/(2*Vr*Vr) + C/(4*Vr*Vr*Vr*Vr) + D/(5*Vr*Vr*Vr*Vr*Vr)
                        + a13/(2*Tr*Tr*Tr*a15) * ( a14 + 1 - (a14+1+a15/(Vr*Vr))*exp(-a15/(Vr*Vr)) );

        phiCH4 = exp(lnphiCH4);

        return(phiCH4);
    }
    template <class Evaluation = Scalar>
    static Evaluation m_iterateVr(Evaluation T, Evaluation p)
    {
        Evaluation A, B, C, D;
        Evaluation Tcrit, pcrit, Tr, pr;
        Evaluation Vr, linkeSeite;
        Evaluation a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
        Evaluation Vc, Vd, Ve, c, d, e, i;
        int n;

            a1 = 0.08725553928;
            a2 = -0.752599476;
            a3 = 0.375419887;
            a4 = 0.0107291342;
            a5 = 0.0054962636;
            a6 = -0.0184772802;
            a7 = 0.000318993183;
            a8 = 0.000211079375;
            a9 = 2.01682801e-05;
            a10 = -1.65606189e-05;
            a11 = 0.000119614546;
            a12 = -0.000108087289;
            a13 = 0.0448262295;
            a14 = 0.75397;
            a15 = 0.077167;

        pcrit = CH4::criticalPressure();
        Tcrit = CH4::criticalTemperature();

        Tr = T/Tcrit;    /* reduced temperature */
        pr = p/pcrit;    /* reduced pressure */

        A = a1 + a2/(Tr*Tr) + a3/(Tr*Tr*Tr);
        B = a4 + a5/(Tr*Tr) + a6/(Tr*Tr*Tr);
        C = a7 + a8/(Tr*Tr) + a9/(Tr*Tr*Tr);
        D = a10 + a11/(Tr*Tr) + a12/(Tr*Tr*Tr);

            Vc = 1.0; Vd = Vc+0.0001;i=-1;
            for (n=1; n<=20; n++)
            {
                    c = Tr/(pr*Vc) * (1 + A/Vc + B/(Vc*Vc) + C/(Vc*Vc*Vc*Vc) + D/(Vc*Vc*Vc*Vc*Vc) + a13/(Tr*Tr*Tr*Vc*Vc)*(a14+a15/(Vc*Vc))*exp(-a15/(Vc*Vc))) - 1.0;
                    d = Tr/(pr*Vd) * (1 + A/Vd + B/(Vd*Vd) + C/(Vd*Vd*Vd*Vd) + D/(Vd*Vd*Vd*Vd*Vd) + a13/(Tr*Tr*Tr*Vd*Vd)*(a14+a15/(Vd*Vd))*exp(-a15/(Vd*Vd))) - 1.0;
                    Ve = Vc - (c*(Vd-Vc))/(d-c);
                    e = Tr/(pr*Ve) * (1 + A/Ve + B/(Ve*Ve) + C/(Ve*Ve*Ve*Ve) + D/(Ve*Ve*Ve*Ve*Ve) + a13/(Tr*Tr*Tr*Ve*Ve)*(a14+a15/(Ve*Ve))*exp(-a15/(Ve*Ve))) - 1.0;
                    if(fabs(e)<1.0E-4) break;
                    Vc=Ve;
                    Vd = Vc+0.0001;
            }
            Vr=Ve;

            if (isnan(Vr))
            {
                    Vc = 0.001; Vd = 85.0; i=-1;
            for (n=1; n<=100; n++)
            {
                    c = Tr/(pr*Vc) * (1 + A/Vc + B/(Vc*Vc) + C/(Vc*Vc*Vc*Vc) + D/(Vc*Vc*Vc*Vc*Vc) + a13/(Tr*Tr*Tr*Vc*Vc)*(a14+a15/(Vc*Vc))*exp(-a15/(Vc*Vc))) - 1.0;
                    d = Tr/(pr*Vd) * (1 + A/Vd + B/(Vd*Vd) + C/(Vd*Vd*Vd*Vd) + D/(Vd*Vd*Vd*Vd*Vd) + a13/(Tr*Tr*Tr*Vd*Vd)*(a14+a15/(Vd*Vd))*exp(-a15/(Vd*Vd))) - 1.0;
                    if (i<n) Ve = (Vc+Vd)/2;
                    else Ve = Vd-((Vd-Vc)*d)/(d-c);
                    e = Tr/(pr*Ve) * (1 + A/Ve + B/(Ve*Ve) + C/(Ve*Ve*Ve*Ve) + D/(Ve*Ve*Ve*Ve*Ve) + a13/(Tr*Tr*Tr*Ve*Ve)*(a14+a15/(Ve*Ve))*exp(-a15/(Ve*Ve))) - 1.0;
                    if(fabs(e)<1.0E-4) break;
                    if ((c*e)<0) Vd=Ve;
                    if ((d*e)<0) Vc=Ve;
                    i=-10.; //-1.2*i;
            }

            }
            Vr=Vc;
            if (fabs(e) > 1.0E-4)
            {
                printf("Vorsicht: keine Loesung fuer Vr gefunden!!! \t p = %.2f bar T = %.2f K \n", p/1.0E5, T);
            }
        return(Vr);
    }

};

}
} // end namespace

#endif

