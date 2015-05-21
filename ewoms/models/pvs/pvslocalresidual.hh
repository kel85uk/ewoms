/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2009-2012 by Klaus Mosthaf

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::PvsLocalResidual
 */
#ifndef EWOMS_PVS_LOCAL_RESIDUAL_HH
#define EWOMS_PVS_LOCAL_RESIDUAL_HH

#include "pvsproperties.hh"

#include <ewoms/models/common/diffusionmodule.hh>
#include <ewoms/models/common/energymodule.hh>

namespace Ewoms {

/*!
 * \ingroup PvsModel
 *
 * \brief Element-wise calculation of the local residual for the
 *        compositional multi-phase primary variable switching model.
 */
template <class TypeTag>
class PvsLocalResidual : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef Ewoms::DiffusionModule<TypeTag, enableDiffusion> DiffusionModule;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    /*!
     * \copydoc ImmiscibleLocalResidual::addPhaseStorage
     */
    template <class LhsEval>
    void addPhaseStorage(Dune::FieldVector<LhsEval, numEq> &storage,
                         const ElementContext &elemCtx,
                         int dofIdx,
                         int timeIdx,
                         int phaseIdx) const
    {
        const IntensiveQuantities &intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto &fs = intQuants.fluidState();

        // compute storage term of all components within all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            int eqIdx = conti0EqIdx + compIdx;
            storage[eqIdx] +=
                Toolbox::template toLhs<LhsEval>(fs.molarity(phaseIdx, compIdx))
                * Toolbox::template toLhs<LhsEval>(fs.saturation(phaseIdx))
                * Toolbox::template toLhs<LhsEval>(intQuants.porosity());
        }

        EnergyModule::addPhaseStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx), phaseIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        storage = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorage(storage, elemCtx, dofIdx, timeIdx, phaseIdx);

        EnergyModule::addSolidHeatStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx));
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux, const ElementContext &elemCtx, int scvfIdx, int timeIdx) const
    {
        flux = Toolbox::createConstant(0.0);
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);

        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          int scvfIdx, int timeIdx) const
    {
        const auto &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        int interiorIdx = extQuants.interiorIndex();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // data attached to upstream and the downstream DOFs
            // of the current phase
            int upIdx = extQuants.upstreamIndex(phaseIdx);
            const IntensiveQuantities &up = elemCtx.intensiveQuantities(upIdx, timeIdx);

            // this is a bit hacky because it is specific to the element-centered
            // finite volume scheme. (N.B. that if finite differences are used to
            // linearize the system of equations, it does not matter.)
            if (upIdx == interiorIdx) {
                Evaluation tmp =
                    up.fluidState().molarDensity(phaseIdx)
                    * extQuants.volumeFlux(phaseIdx);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    flux[conti0EqIdx + compIdx] +=
                        tmp*up.fluidState().moleFraction(phaseIdx, compIdx);
                }
            }
            else {
                Evaluation tmp =
                    Toolbox::value(up.fluidState().molarDensity(phaseIdx))
                    * extQuants.volumeFlux(phaseIdx);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    flux[conti0EqIdx + compIdx] +=
                        tmp*Toolbox::value(up.fluidState().moleFraction(phaseIdx, compIdx));
                }
            }
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          int scvfIdx, int timeIdx) const
    {
        DiffusionModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeSource
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int dofIdx,
                       int timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(source);
    }
};

} // namespace Ewoms

#endif
