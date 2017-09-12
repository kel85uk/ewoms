#ifndef NACL_HH
#define NACL_HH


/*!
 * \file
 *
 * \brief A class for the NaCl properties
 */
#ifndef OPM_NACL_HPP
#define OPM_NACL_HPP

#include <opm/common/Exceptions.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/common/MathToolbox.hpp>
//#include <cmath>
#include <iostream>


namespace Opm
{
/*!
 * \brief A class for the NACL fluid properties for tracer purpose
 */
template <class Scalar>
class NaCl : public Component<Scalar, NaCl<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the NaCl.
     */
    static const char *name()
    { return "NaCl"; }

    /*!
     * \brief The mass in [kg] of one mole of NaCl.
     */
    static Scalar molarMass()
    { return 58.4428e-3 ; } // kg/mol

    /*!
     * \brief The diffusion Coefficient of NaCl in water.
     */
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    { return 2e-9; }

    static Scalar charge()
    {
        return 0.0;
    }

    static Scalar Density()
    {
        return (2165); /* 2165 kg/m³*/
    }
};

} // end namespace

#endif



#endif // NACL_HH

