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
/**
 * \file
 *
 * \copydoc Ewoms::EclCarterTracyAquifer
 */
#ifndef EWOMS_ECL_AQUIFER_CARTER_TRACY_HH
#define EWOMS_ECL_AQUIFER_CARTER_TRACY_HH

#include <ewoms/models/blackoil/blackoilproperties.hh>
#include <ewoms/aux/baseauxiliarymodule.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/alignedallocator.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/common/Valgrind.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>

#include <map>

namespace Ewoms {

template <class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief The aquifer model of Peaceman.
 *
 * This class is tailored for the element centered finite volume
 * discretization, assumes a vertical borehole and is intended to be
 * used by the EclAquiferManager.
 *
 * See:
 *
 * Z. Chen, G. Huan, Y. Ma: Computational Methods for Multiphase
 * Flows in Porous Media, 1st edition, SIAM, 2006, pp. 445-446
 *
 * and
 *
 * D. W. Peaceman: Interpretation of aquifer-block pressures in numerical
 * reservoir simulation, The 52nd Annual SPE Fall Technical Conference
 * and Exhibition, Denver, CO., 1977
 */
template <class TypeTag>
class EclCarterTracyAquifer : public BaseAuxiliaryModule<TypeTag>
{
    typedef BaseAuxiliaryModule<TypeTag> AuxModule;

    typedef typename AuxModule::NeighborSet NeighborSet;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef Element  ElementStorage;

    // the dimension of the simulator's world
    static const int dimWorld = GridView::dimensionworld;

    // convenient access to the number of phases and the number of
    // components
    static const unsigned numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static const unsigned numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

    // convenient access to the phase and component indices. If the compiler bails out
    // here, you're probably using an incompatible fluid system. This class has only been
    // tested with Opm::FluidSystems::BlackOil...
    static const unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static const unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static const unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static const unsigned oilCompIdx = FluidSystem::oilCompIdx;
    static const unsigned waterCompIdx = FluidSystem::waterCompIdx;
    static const unsigned gasCompIdx = FluidSystem::gasCompIdx;

    static const unsigned numModelEq = GET_PROP_VALUE(TypeTag, NumEq);
    static const unsigned conti0EqIdx = GET_PROP_TYPE(TypeTag, Indices)::conti0EqIdx;

    typedef Opm::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> FluidState;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    // all quantities that need to be stored per degree of freedom that intersects the
    // aquifer.
    struct DofVariables {
        DofVariables() = default;
        DofVariables(const DofVariables&) = default;

        // retrieve the solution dependent quantities that are only updated at the
        // beginning of a time step from the IntensiveQuantities of the model
        void updateBeginTimestep(const IntensiveQuantities& intQuants OPM_UNUSED)
        {}

        // retrieve the solution dependent quantities from the IntensiveQuantities of the
        // model
        void update(const IntensiveQuantities& intQuants)
        {
            const auto& fs = intQuants.fluidState();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                pressure[phaseIdx] = fs.pressure(phaseIdx);
                density[phaseIdx] = fs.density(phaseIdx);
                mobility[phaseIdx] = intQuants.mobility(phaseIdx);
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                oilMassFraction[compIdx] = fs.massFraction(oilPhaseIdx, compIdx);
                gasMassFraction[compIdx] = fs.massFraction(gasPhaseIdx, compIdx);
            }
        }

        // the depth of the centroid of the DOF
        Scalar depth;

        // the volume in m^3 of the DOF
        Scalar totalVolume;

        // the effective size of an element in each direction. This is defined as the
        // distance of the face centers along the respective axis.
        std::array<Scalar, dimWorld> effectiveSize;

        // the intrinsic permeability matrix for the degree of freedom
        DimMatrix permeability;

        // the effective permeability of the connection. usually that's the geometric
        // mean of the X and Y permeabilities of the DOF times the DOF's height
        Scalar effectivePermeability;

        // The connection transmissibility factor to be used for a given DOF. this is
        // usually computed from the values above but it can be explicitly specified by
        // the user...
        Scalar connectionTransmissibilityFactor;

        // the radius of the aquifer for the given degree of freedom
        Scalar boreholeRadius;

        // The skin factor of the aquifer at the given degree of freedom
        Scalar skinFactor;

        //////////////
        // the following quantities depend on the considered solution and are thus updated
        // at the beginning of each Newton-Raphson iteration.
        //////////////

        // the phase pressures inside a DOF
        std::array<Evaluation, numPhases> pressure;

        // the phase densities at the DOF
        std::array<Evaluation, numPhases> density;

        // the phase mobilities of the DOF
        std::array<Evaluation, numPhases> mobility;

        // the composition of the oil phase at the DOF
        std::array<Evaluation, numComponents> oilMassFraction;

        // the composition of the gas phase at the DOF
        std::array<Evaluation, numComponents> gasMassFraction;

        ElementStorage element;
        unsigned pvtRegionIdx;
        unsigned localDofIdx;
    };

    // some safety checks/caveats
    static_assert(std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value,
                  "The Carter Tracy aquifer model is only implemented for the "
                  "element-centered finite volume discretization!");
    static_assert(dimWorld == 3,
                  "The Carter Tracy aquifer model is only implemented for 3D grids!");

public:
    enum ControlMode {
        BottomHolePressure,
        TubingHeadPressure,
        VolumetricSurfaceRate,
        VolumetricReservoirRate
    };

    enum AquiferType {
        Undefined,
        Injector,
        Producer
    };

    enum AquiferStatus {
        // production/injection is ongoing
        Open,

        // no production/injection, but aquifer is only closed above the reservoir, so cross
        // flow is possible
        Closed,

        // aquifer is completely separated from the reservoir, e.g. by filling it with
        // concrete.
        Shut
    };

    EclCarterTracyAquifer(const Simulator& simulator)
        : simulator_(simulator)
    {
        // set the initial status of the aquifer
        aquiferType_ = Undefined;
        aquiferStatus_ = Shut;
        controlMode_ = BottomHolePressure;

        aquiferTotalVolume_ = 0.0;

        bhpLimit_ = 0.0;
        thpLimit_ = 0.0;

        targetBottomHolePressure_ = 0.0;
        actualBottomHolePressure_ = 0.0;
        maximumSurfaceRate_ = 0.0;
        maximumReservoirRate_ = 0.0;

        actualWeightedSurfaceRate_ = 0.0;
        actualWeightedResvRate_ = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            actualSurfaceRates_[phaseIdx] = 0.0;
            actualResvRates_[phaseIdx] = 0.0;

            volumetricWeight_[phaseIdx] = 0.0;
        }

        refDepth_ = 0.0;

        // set the composition of the injected fluids based. If
        // somebody is stupid enough to inject oil, we assume he wants
        // to loose his fortune on dry oil...
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
                injectionFluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);
        injectionFluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
        injectionFluidState_.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
        injectionFluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1.0);

        // set the temperature to 25 deg C, just so that it is set
        injectionFluidState_.setTemperature(273.15 + 25);

        injectedPhaseIdx_ = oilPhaseIdx;
    }

    /*!
     * \copydoc Ewoms::BaseAuxiliaryModule::numDofs()
     */
    virtual unsigned numDofs() const
    { return 1; }

    /*!
     * \copydoc Ewoms::BaseAuxiliaryModule::addNeighbors()
     */
    virtual void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        int aquiferGlobalDof = AuxModule::localToGlobalDof(/*localDofIdx=*/0);

        // the aquifer's bottom hole pressure always affects itself...
        neighbors[aquiferGlobalDof].insert(aquiferGlobalDof);

        // add the grid DOFs which are influenced by the aquifer, and add the aquifer dof to
        // the ones neighboring the grid ones
        auto aquiferDofIt = dofVariables_.begin();
        const auto& aquiferDofEndIt = dofVariables_.end();
        for (; aquiferDofIt != aquiferDofEndIt; ++ aquiferDofIt) {
            neighbors[aquiferGlobalDof].insert(aquiferDofIt->first);
            neighbors[aquiferDofIt->first].insert(aquiferGlobalDof);
        }
    }

    /*!
     * \copydoc Ewoms::BaseAuxiliaryModule::addNeighbors()
     */
    virtual void applyInitial()
    {
        auto& sol = const_cast<SolutionVector&>(simulator_.model().solution(/*timeIdx=*/0));

        int aquiferGlobalDof = AuxModule::localToGlobalDof(/*localDofIdx=*/0);
        sol[aquiferGlobalDof] = 0.0;
    }

    /*!
     * \copydoc Ewoms::BaseAuxiliaryModule::linearize()
     */
    virtual void linearize(JacobianMatrix& matrix, GlobalEqVector& residual)
    {
        const SolutionVector& curSol = simulator_.model().solution(/*timeIdx=*/0);

        unsigned aquiferGlobalDofIdx = AuxModule::localToGlobalDof(/*localDofIdx=*/0);
        residual[aquiferGlobalDofIdx] = 0.0;

        auto& diagBlock = matrix[aquiferGlobalDofIdx][aquiferGlobalDofIdx];
        diagBlock = 0.0;
        for (unsigned i = 0; i < numModelEq; ++ i)
            diagBlock[i][i] = 1.0;

        if (aquiferStatus() == Shut) {
            // if the aquifer is shut, make the auxiliary DOFs a trivial equation in the
            // matrix: the main diagonal is already set to the identity matrix, the
            // off-diagonal matrix entries must be set to 0.
            auto aquiferDofIt = dofVariables_.begin();
            const auto& aquiferDofEndIt = dofVariables_.end();
            for (; aquiferDofIt != aquiferDofEndIt; ++ aquiferDofIt) {
                matrix[aquiferGlobalDofIdx][aquiferDofIt->first] = 0.0;
                matrix[aquiferDofIt->first][aquiferGlobalDofIdx] = 0.0;
                residual[aquiferGlobalDofIdx] = 0.0;
            }
            return;
        }

        Scalar aquiferResid = aquiferResidual_(actualBottomHolePressure_);
        residual[aquiferGlobalDofIdx][0] = aquiferResid;

        // account for the effect of the grid DOFs which are influenced by the aquifer on
        // the aquifer equation and the effect of the aquifer on the grid DOFs
        auto aquiferDofIt = dofVariables_.begin();
        const auto& aquiferDofEndIt = dofVariables_.end();

        ElementContext elemCtx(simulator_);
        for (; aquiferDofIt != aquiferDofEndIt; ++ aquiferDofIt) {
            unsigned gridDofIdx = aquiferDofIt->first;
            const auto& dofVars = *dofVariables_[gridDofIdx];
            DofVariables tmpDofVars(dofVars);
            auto priVars(curSol[gridDofIdx]);

            /////////////
            // influence of grid on aquifer
            auto& curBlock = matrix[aquiferGlobalDofIdx][gridDofIdx];

            elemCtx.updateStencil( dofVars.element );
            curBlock = 0.0;
            for (unsigned priVarIdx = 0; priVarIdx < numModelEq; ++priVarIdx) {
                // calculate the derivative of the aquifer equation w.r.t. the current
                // primary variable using forward differences
                Scalar eps =
                    1e3
                    *std::numeric_limits<Scalar>::epsilon()
                    *std::max<Scalar>(1.0, priVars[priVarIdx]);
                priVars[priVarIdx] += eps;

                elemCtx.updateIntensiveQuantities(priVars, dofVars.localDofIdx, /*timeIdx=*/0);
                tmpDofVars.update(elemCtx.intensiveQuantities(dofVars.localDofIdx, /*timeIdx=*/0));

                Scalar dAquiferEq_dPV =
                    (aquiferResidual_(actualBottomHolePressure_, &tmpDofVars, gridDofIdx) - aquiferResid)
                    / eps;
                curBlock[0][priVarIdx] = dAquiferEq_dPV;

                // go back to the original primary variables
                priVars[priVarIdx] -= eps;
            }
            //
            /////////////

            /////////////
            // influence of aquifer on grid:
            RateVector q(0.0);
            RateVector modelRate;
            std::array<Scalar, numPhases> resvRates;

            elemCtx.updateIntensiveQuantities(priVars, dofVars.localDofIdx, /*timeIdx=*/0);

            const auto& fluidState = elemCtx.intensiveQuantities(dofVars.localDofIdx, /*timeIdx=*/0).fluidState();

            // first, we need the source term of the grid for the slightly disturbed aquifer.
            Scalar eps =
                1e3
                *std::numeric_limits<Scalar>::epsilon()
                *std::max<Scalar>(1e5, actualBottomHolePressure_);
            computeVolumetricDofRates_(resvRates, actualBottomHolePressure_ + eps, *dofVariables_[gridDofIdx]);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                modelRate.setVolumetricRate(fluidState, phaseIdx, resvRates[phaseIdx]);
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    q[compIdx] += modelRate[compIdx];
            }

            // then, we subtract the source rates for a undisturbed aquifer.
            computeVolumetricDofRates_(resvRates, actualBottomHolePressure_, *dofVariables_[gridDofIdx]);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                modelRate.setVolumetricRate(fluidState, phaseIdx, resvRates[phaseIdx]);
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    q[compIdx] -= modelRate[compIdx];
            }

            // and finally, we divide by the epsilon to get the derivative
            for (unsigned eqIdx = 0; eqIdx < numModelEq; ++eqIdx)
                q[eqIdx] /= eps;

            // now we put this derivative into the right place in the Jacobian
            // matrix. This is a bit hacky because it assumes that the model uses a mass
            // rate for each component as its first conservation equation, but we require
            // the black-oil model for now anyway, so this should not be too much of a
            // problem...
            Opm::Valgrind::CheckDefined(q);
            auto& matrixEntry = matrix[gridDofIdx][aquiferGlobalDofIdx];
            matrixEntry = 0.0;
            for (unsigned eqIdx = 0; eqIdx < numModelEq; ++ eqIdx)
                matrixEntry[eqIdx][0] = - Toolbox::value(q[eqIdx])/dofVars.totalVolume;
            //
            /////////////
        }

        // effect of changing the aquifer's bottom hole pressure on the aquifer equation
        Scalar eps =
            1e3
            *std::numeric_limits<Scalar>::epsilon()
            *std::max<Scalar>(1e7, targetBottomHolePressure_);
        Scalar aquiferResidStar = aquiferResidual_(actualBottomHolePressure_ + eps);
        diagBlock[0][0] = (aquiferResidStar - aquiferResid)/eps;
    }


    // reset the aquifer to the initial state, i.e. remove all degrees of freedom...
    void clear()
    {
        dofVarsStore_.clear();
        dofVariables_.clear();
    }

    /*!
     * \brief Begin the specification of the aquifer.
     *
     * The specification process is the following:
     *
     * beginSpec()
     * setName("FOO");
     * // add degrees of freedom to the aquifer
     * for (dof in aquiferDofs)
     *    addDof(dof);
     * endSpec()
     *
     * // set the radius of the aquifer at the dof [m].
     * // optional, if not specified, it is assumed to be 0.1524m
     * setRadius(dof, someRadius);
     *
     * // set the skin factor of the aquifer.
     * // optional, if not specified, it is assumed to be 0
     * setSkinFactor(dof, someSkinFactor);
     *
     * // specify the phase which is supposed to be injected. (Optional,
     * // if unspecified, the aquifer will throw an
     * // exception if it would inject something.)
     * setInjectedPhaseIndex(phaseIdx);
     *
     * // set maximum production rate at reservoir conditions
     * // (kg/s, optional, if not specified, the aquifer is assumed to be
     * // shut for production)
     * setMaximumReservoirRate(someMassRate);
     *
     * // set maximum injection rate at reservoir conditions
     * // (kg/s, optional, if not specified, the aquifer is assumed to be
     * // shut for injection)
     * setMinmumReservoirRate(someMassRate);
     *
     * // set the relative weight of the mass rate of a fluid phase.
     * // (Optional, if unspecified each phase exhibits a weight of 1)
     * setPhaseWeight(phaseIdx, someWeight);
     *
     * // set maximum production rate at surface conditions
     * // (kg/s, optional, if not specified, the aquifer is assumed to be
     * // not limited by the surface rate)
     * setMaximumSurfaceRate(someMassRate);
     *
     * // set maximum production rate at surface conditions
     * // (kg/s, optional, if not specified, the aquifer is assumed to be
     * // not limited by the surface rate)
     * setMinimumSurfaceRate(someMassRate);
     *
     * // set the minimum pressure at the bottom of the aquifer (Pa,
     * // optional, if not specified, the aquifer is assumes it estimates
     * // the bottom hole pressure based on the tubing head pressure
     * // assuming hydrostatic conditions.)
     * setMinimumBottomHolePressure(somePressure);
     *
     * // set the pressure at the top of the aquifer (Pa,
     * // optional, if not specified, the tubing head pressure is
     * // assumed to be 1 bar)
     * setTubingHeadPressure(somePressure);
     *
     * // set the control mode of the aquifer [m].
     * // optional, if not specified, it is assumed to be "BottomHolePressure"
     * setControlMode(Aquifer::TubingHeadPressure);
     *
     * // set the tubing head pressure of the aquifer [Pa]
     * // only require  if the control mode is "TubingHeadPressure"
     * setTubingHeadPressure(1e5);
     */
    void beginSpec()
    {
        // this is going to be set to a real value by any realistic grid. Shall we bet?
        refDepth_ = 1e100;

        // By default, take the bottom hole pressure as a given
        controlMode_ = ControlMode::BottomHolePressure;

        // use one bar for the default bottom hole and tubing head
        // pressures. For the bottom hole pressure, this is probably
        // off by at least one magnitude...
        bhpLimit_ = 1e5;
        thpLimit_ = 1e5;

        // reset the actually observed bottom hole pressure
        actualBottomHolePressure_ = 0.0;

        // By default, all fluids exhibit the weight 1.0
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            volumetricWeight_[phaseIdx] = 1.0;

        aquiferType_ = Undefined;

        aquiferTotalVolume_ = 0.0;
    }

    /*!
     * \brief Set the relative weight of the volumetric phase rates.
     */
    void setVolumetricPhaseWeights(Scalar oilWeight, Scalar gasWeight, Scalar waterWeight)
    {
        volumetricWeight_[oilPhaseIdx] = oilWeight;
        volumetricWeight_[gasPhaseIdx] = gasWeight;
        volumetricWeight_[waterPhaseIdx] = waterWeight;
    }

    /*!
     * \brief Return the human-readable name of the aquifer
     *
     * Aquifer, let's say "readable by some humans".
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Set the human-readable name of the aquifer
     */
    void setName(const std::string& newName)
    { name_ = newName; }

    /*!
     * \brief Add a degree of freedom to the aquifer.
     */
    template <class Context>
    void addDof(const Context& context, unsigned dofIdx)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        if (applies(globalDofIdx))
            // we already have this DOF in the aquifer!
            return;

        const auto& dofPos = context.pos(dofIdx, /*timeIdx=*/0);

        dofVarsStore_.push_back(DofVariables());
        dofVariables_[globalDofIdx] = &dofVarsStore_.back();
        DofVariables& dofVars = *dofVariables_[globalDofIdx];
        aquiferTotalVolume_ += context.model().dofTotalVolume(globalDofIdx);

        dofVars.element = context.element();

        dofVars.localDofIdx = dofIdx;
        dofVars.pvtRegionIdx = context.problem().pvtRegionIndex(context, dofIdx, /*timeIdx=*/0);
        assert(dofVars.pvtRegionIdx == 0);

        // determine the size of the element
        dofVars.effectiveSize.fill(0.0);

        // we assume all elements to be hexahedrons!
        assert(context.element().subEntities(/*codim=*/dimWorld) == 8);

        const auto& refElem = Dune::ReferenceElements<Scalar, /*dim=*/3>::cube();

        // determine the current element's effective size
        const auto& elem = context.element();
        unsigned faceIdx = 0;
        unsigned numFaces = refElem.size(/*codim=*/1);
        for (; faceIdx < numFaces; ++faceIdx) {
            const auto& faceCenterLocal = refElem.position(faceIdx, /*codim=*/1);
            const auto& faceCenter = elem.geometry().global(faceCenterLocal);

            switch (faceIdx) {
            case 0:
                dofVars.effectiveSize[0] -= faceCenter[0];
                break;
            case 1:
                dofVars.effectiveSize[0] += faceCenter[0];
                break;
            case 2:
                dofVars.effectiveSize[1] -= faceCenter[1];
                break;
            case 3:
                dofVars.effectiveSize[1] += faceCenter[1];
                break;
            case 4:
                dofVars.depth += faceCenter[2];
                dofVars.effectiveSize[2] -= faceCenter[2];
                break;
            case 5:
                dofVars.depth += faceCenter[2];
                dofVars.effectiveSize[2] += faceCenter[2];
                break;
            }
        }

        // the volume associated with the DOF
        dofVars.totalVolume = context.model().dofTotalVolume(globalDofIdx);

        // the depth of the degree of freedom
        dofVars.depth /= 2;

        // default borehole radius: 1/2 foot
        dofVars.boreholeRadius = 0.3048/2;

        // default skin factor: 0
        dofVars.skinFactor = 0;

        // the intrinsic permeability tensor of the DOF
        const auto& K = context.problem().intrinsicPermeability(context, dofIdx, /*timeIdx=*/0);
        dofVars.permeability = K;

        // default the effective permeability: Geometric mean of the x and y components
        // of the intrinsic permeability of DOF times the DOF's height.
        assert(K[0][0] > 0);
        assert(K[1][1] > 0);
        dofVars.effectivePermeability =
            std::sqrt(K[0][0]*K[1][1])*dofVars.effectiveSize[2];

        // from that, compute the default connection transmissibility factor
        computeConnectionTransmissibilityFactor_(globalDofIdx);

        // we assume that the z-coordinate represents depth (and not
        // height) here...
        if (dofPos[2] < refDepth_)
            refDepth_ = dofPos[2];
    }

    /*!
     * \brief Finalize the specification of the borehole.
     */
    void endSpec()
    {
        const auto& comm = simulator_.gridView().comm();

        if (dofVariables_.size() == 0) {
            std::cout << "Aquifer " << name() << " does not penetrate any active cell."
                      << " Assuming it to be shut!\n";
            setAquiferStatus(AquiferStatus::Shut);
            return;
        }

        // determine the maximum depth of the aquifer over all processes
        refDepth_ = comm.min(refDepth_);

        // the total volume of the aquifer must also be summed over all processes
        aquiferTotalVolume_ = comm.sum(aquiferTotalVolume_);
    }

    /*!
     * \brief Set the control mode of the aquifer.
     *
     * This specifies which quantities are assumed to be externally
     * given and which must be calculated based on those.
     */
    void setControlMode(ControlMode controlMode)
    { controlMode_ = controlMode; }

    /*!
     * \brief Set the connection transmissibility factor for a given degree of freedom.
     */
    template <class Context>
    void setConnectionTransmissibilityFactor(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx]->connectionTransmissibilityFactor = value;
    }

    /*!
     * \brief Set the effective permeability Kh to be used for a given degree of freedom.
     *
     * By default, Kh is sqrt(K_xx * K_yy) * h, where K_xx and K_yy is the permeability
     * for the DOF in X and Y directions and h is the height associated with the degree
     * of freedom.
     *
     * Note: The connection transmissibility factor is updated after calling this method,
     *       so if setConnectionTransmissibilityFactor() is to have any effect, it should
     *       be called after setEffectivePermeability()!
     */
    template <class Context>
    void setEffectivePermeability(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx].effectivePermeability = value;

        computeConnectionTransmissibilityFactor_(globalDofIdx);
    }

    /*!
     * \brief Set the type of the aquifer (i.e., injector or producer).
     */
    void setAquiferType(AquiferType aquiferType)
    { aquiferType_ = aquiferType; }

    /*!
     * \brief Returns the type of the aquifer (i.e., injector or producer).
     */
    AquiferType aquiferType() const
    { return aquiferType_; }

    /*!
     * \brief Set the index of fluid phase to be injected.
     *
     * This is only relevant if the aquifer type is an injector.
     */
    void setInjectedPhaseIndex(unsigned injPhaseIdx)
    { injectedPhaseIdx_ = injPhaseIdx; }

    /*!
     * \brief Sets the reference depth for the bottom hole pressure [m]
     */
    void setReferenceDepth(Scalar value)
    { refDepth_ = value; }

    /*!
     * \brief The reference depth for the bottom hole pressure [m]
     */
    Scalar referenceDepth() const
    { return refDepth_; }

    /*!
     * \brief Set whether the aquifer is open,closed or shut
     */
    void setAquiferStatus(AquiferStatus status)
    { aquiferStatus_ = status; }

    /*!
     * \brief Return whether the aquifer is open,closed or shut
     */
    AquiferStatus aquiferStatus() const
    { return aquiferStatus_; }

    /*!
     * \brief Return true iff a degree of freedom is directly affected
     *        by the aquifer
     */
    bool applies(unsigned globalDofIdx) const
    { return dofVariables_.count(globalDofIdx) > 0; }

    /*!
     * \brief Set the maximum/minimum bottom hole pressure [Pa] of the aquifer.
     */
    void setTargetBottomHolePressure(Scalar val)
    { bhpLimit_ = val; }

    /*!
     * \brief Return the maximum/minimum bottom hole pressure [Pa] of the aquifer.
     *
     * For injectors, this is the maximum, for producers it's the minimum.
     */
    Scalar targetBottomHolePressure() const
    { return bhpLimit_; }

    /*!
     * \brief Return the maximum/minimum bottom hole pressure [Pa] of the aquifer.
     */
    Scalar bottomHolePressure() const
    { return actualBottomHolePressure_; }

    /*!
     * \brief Set the tubing head pressure [Pa] of the aquifer.
     */
    void setTargetTubingHeadPressure(Scalar val)
    { thpLimit_ = val; }

    /*!
     * \brief Return the maximum/minimum tubing head pressure [Pa] of the aquifer.
     *
     * For injectors, this is the maximum, for producers it's the minimum.
     */
    Scalar targetTubingHeadPressure() const
    { return thpLimit_; }

    /*!
     * \brief Return the maximum/minimum tubing head pressure [Pa] of the aquifer.
     */
    Scalar tubingHeadPressure() const
    {
        // warning: this is a bit hacky...
        Scalar rho = 650; // kg/m^3
        Scalar g = 9.81; // m/s^2
        return actualBottomHolePressure_ + rho*refDepth_*g;
    }

    /*!
     * \brief Set the maximum combined rate of the fluids at the surface.
     */
    void setMaximumSurfaceRate(Scalar value)
    { maximumSurfaceRate_ = value; }

    /*!
     * \brief Return the weighted maximum surface rate [m^3/s] of the aquifer.
     */
    Scalar maximumSurfaceRate() const
    { return maximumSurfaceRate_; }

    /*!
     * \brief Set the maximum combined rate of the fluids at the surface.
     */
    void setMaximumReservoirRate(Scalar value)
    { maximumReservoirRate_ = value; }

    /*!
     * \brief Return the weighted maximum reservoir rate [m^3/s] of the aquifer.
     */
    Scalar maximumReservoirRate() const
    { return maximumReservoirRate_; }

    /*!
     * \brief Return the reservoir rate [m^3/s] actually seen by the aquifer in the current time
     *        step.
     */
    Scalar reservoirRate() const
    { return actualWeightedResvRate_; }

    /*!
     * \brief Return the weighted surface rate [m^3/s] actually seen by the aquifer in the current time
     *        step.
     */
    Scalar surfaceRate() const
    { return actualWeightedSurfaceRate_; }

    /*!
     * \brief Return the reservoir rate [m^3/s] of a given fluid which is actually seen
     *        by the aquifer in the current time step.
     */
    Scalar reservoirRate(unsigned phaseIdx) const
    { return actualResvRates_[phaseIdx]; }

    /*!
     * \brief Return the weighted surface rate [m^3/s] of a given fluid which is actually
     *        seen by the aquifer in the current time step.
     */
    Scalar surfaceRate(unsigned phaseIdx) const
    { return actualSurfaceRates_[phaseIdx]; }

    /*!
     * \brief Set the skin factor of the aquifer
     *
     * Note: The connection transmissibility factor is updated after calling this method,
     *       so if setConnectionTransmissibilityFactor() is to have any effect, it should
     *       be called after setSkinFactor()!
     */
    template <class Context>
    void setSkinFactor(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx].skinFactor = value;

        computeConnectionTransmissibilityFactor_(globalDofIdx);
    }

    /*!
     * \brief Return the aquifer's skin factor at a DOF [-].
     */
    Scalar skinFactor(unsigned gridDofIdx) const
    { return dofVariables_.at(gridDofIdx).skinFactor_; }

    /*!
     * \brief Set the borehole radius of the aquifer
     *
     * Note: The connection transmissibility factor is updated after calling this method,
     *       so if setConnectionTransmissibilityFactor() is to have any effect, it should
     *       be called after setRadius()!
     */
    template <class Context>
    void setRadius(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx]->boreholeRadius = value;

        computeConnectionTransmissibilityFactor_(globalDofIdx);
    }

    /*!
     * \brief Return the aquifer's radius at a cell [m].
     */
    Scalar radius(unsigned gridDofIdx) const
    { return dofVariables_.at(gridDofIdx)->radius_; }

    /*!
     * \brief Informs the aquifer that a time step has just begun.
     */
    void beginTimeStep()
    {
        if (aquiferStatus() == Shut)
            return;

        // calculate the bottom hole pressure to be actually used
        if (controlMode_ == ControlMode::TubingHeadPressure) {
            // assume a density of 650 kg/m^3 for the bottom hole pressure
            // calculation
            Scalar rho = 650.0;
            targetBottomHolePressure_ = thpLimit_ + rho*refDepth_;
        }
        else if (controlMode_ == ControlMode::BottomHolePressure)
            targetBottomHolePressure_ = bhpLimit_;
        else
            // TODO: also take the tubing head pressure limit into account...
            targetBottomHolePressure_ = bhpLimit_;

        // make it very likely that we screw up if we control for {surface,reservoir}
        // rate, but depend on the {reservoir,surface} rate somewhere...
        if (controlMode_ == ControlMode::VolumetricSurfaceRate)
            maximumReservoirRate_ = 1e100;
        else if (controlMode_ == ControlMode::VolumetricReservoirRate)
            maximumSurfaceRate_ = 1e100;

        // reset the iteration index
        iterationIdx_ = 0;
    }

    /*!
     * \brief Informs the aquifer that an iteration has just begun.
     *
     * The beginIteration*() methods, the aquifer calculates the bottom
     * and tubing head pressures, the actual unconstraint production and
     * injection rates, etc. The callback is split into three parts as
     * this arrangement avoids iterating over the whole grid and to
     * re-calculate the volume variables for each aquifer.
     *
     * This is supposed to prepare the aquifer object to do the
     * computations which are required to do the DOF specific
     * things.
     */
    void beginIterationPreProcess()
    { }

    /*!
     * \brief Do the DOF specific part at the beginning of each iteration
     */
    template <class Context>
    void beginIterationAccumulate(Context& context, unsigned timeIdx)
    {
        if (aquiferStatus() == Shut)
            return;

        for (unsigned dofIdx = 0; dofIdx < context.numPrimaryDof(timeIdx); ++dofIdx) {
            unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
            if (!applies(globalDofIdx))
                continue;

            DofVariables& dofVars = *dofVariables_.at(globalDofIdx);
            const auto& intQuants = context.intensiveQuantities(dofIdx, timeIdx);

            if (iterationIdx_ == 0)
                dofVars.updateBeginTimestep(intQuants);

            dofVars.update(intQuants);
        }
    }

    /*!
     * \brief Informs the aquifer that an iteration has just begun.
     *
     * This is the post-processing part which uses the results of the
     * accumulation callback.
     */
    void beginIterationPostProcess()
    {
        if (aquiferStatus() == Shut)
            return;

        auto& sol = const_cast<SolutionVector&>(simulator_.model().solution(/*timeIdx=*/0));
        int aquiferGlobalDof = AuxModule::localToGlobalDof(/*localDofIdx=*/0);

        // retrieve the bottom hole pressure from the global system of equations
        actualBottomHolePressure_ = Toolbox::value(dofVariables_.begin()->second->pressure[0]);
        actualBottomHolePressure_ = computeRateEquivalentBhp_();

        sol[aquiferGlobalDof][0] = actualBottomHolePressure_;

        computeOverallRates_(actualBottomHolePressure_,
                             actualResvRates_,
                             actualSurfaceRates_);

        actualWeightedResvRate_ = computeWeightedRate_(actualResvRates_);
        actualWeightedSurfaceRate_ = computeWeightedRate_(actualSurfaceRates_);
    }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    { ++ iterationIdx_; }

    /*!
     * \brief Called by the simulator after each time step.
     */
    void endTimeStep()
    {
        if (aquiferStatus() == Shut)
            return;

        // we use a condition that is always false here to prevent the code below from
        // bitrotting. (i.e., at least it stays compileable)
        if (false && simulator_.gridView().comm().rank() == 0) {
            std::cout << "Aquifer '" << name() << "':\n";
            std::cout << " Control mode: " << controlMode_ << "\n";
            std::cout << " BHP limit: " << bhpLimit_/1e5 << " bar\n";
            std::cout << " Observed BHP: " << actualBottomHolePressure_/1e5 << " bar\n";
            std::cout << " Weighted surface rate limit: " << maximumSurfaceRate_ << "\n";
            std::cout << " Weighted surface rate: " << std::abs(actualWeightedSurfaceRate_) << " (="
                      << 100*std::abs(actualWeightedSurfaceRate_)/maximumSurfaceRate_ << "%)\n";

            std::cout << " Surface rates:\n";
            std::cout << "  oil: "
                      << actualSurfaceRates_[oilPhaseIdx] << " m^3/s = "
                      << actualSurfaceRates_[oilPhaseIdx]*(24*60*60) << " m^3/day = "
                      << actualSurfaceRates_[oilPhaseIdx]*(24*60*60)/0.15898729 << " STB/day = "
                      << actualSurfaceRates_[oilPhaseIdx]*(24*60*60)
                         *FluidSystem::referenceDensity(oilPhaseIdx, /*pvtRegionIdx=*/0) << " kg/day"
                      << "\n";
            std::cout << "  gas: "
                      << actualSurfaceRates_[gasPhaseIdx] << " m^3/s = "
                      << actualSurfaceRates_[gasPhaseIdx]*(24*60*60) << " m^3/day = "
                      << actualSurfaceRates_[gasPhaseIdx]*(24*60*60)/28.316847 << " MCF/day = "
                      << actualSurfaceRates_[gasPhaseIdx]*(24*60*60)
                         *FluidSystem::referenceDensity(gasPhaseIdx, /*pvtRegionIdx=*/0) << " kg/day"
                      << "\n";
            std::cout << "  water: "
                      << actualSurfaceRates_[waterPhaseIdx] << " m^3/s = "
                      << actualSurfaceRates_[waterPhaseIdx]*(24*60*60) << " m^3/day = "
                      << actualSurfaceRates_[waterPhaseIdx]*(24*60*60)/0.15898729 << " STB/day = "
                      << actualSurfaceRates_[waterPhaseIdx]*(24*60*60)
                         *FluidSystem::referenceDensity(waterPhaseIdx, /*pvtRegionIdx=*/0) << " kg/day"
                      << "\n";
        }
    }

    /*!
     * \brief Computes the source term for a degree of freedom.
     */
    template <class Context>
    void computeTotalRatesForDof(RateVector& q,
                                 const Context& context,
                                 unsigned dofIdx,
                                 unsigned timeIdx) const
    {
        q = 0.0;

        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        if (aquiferStatus() == Shut || !applies(globalDofIdx))
            return;

        // create a DofVariables object for the current evaluation point
        DofVariables tmp(*dofVariables_.at(globalDofIdx));

        tmp.update(context.intensiveQuantities(dofIdx, timeIdx));

        std::array<Evaluation, numPhases> volumetricRates;
        computeVolumetricDofRates_(volumetricRates, actualBottomHolePressure_, tmp);

        // convert to mass rates
        RateVector modelRate;
        const auto& intQuants = context.intensiveQuantities(dofIdx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            modelRate.setVolumetricRate(intQuants.fluidState(), phaseIdx, volumetricRates[phaseIdx]);
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                q[conti0EqIdx + compIdx] += modelRate[conti0EqIdx + compIdx];
        }

        Opm::Valgrind::CheckDefined(q);
    }

protected:
    // compute the connection transmissibility factor based on the effective permeability
    // of a connection, the radius of the borehole and the skin factor.
    void computeConnectionTransmissibilityFactor_(unsigned globalDofIdx)
    {
        auto& dofVars = *dofVariables_[globalDofIdx];

        const auto& D = dofVars.effectiveSize;
        const auto& K = dofVars.permeability;
        Scalar Kh = dofVars.effectivePermeability;
        Scalar S = dofVars.skinFactor;
        Scalar rAquifer = dofVars.boreholeRadius;

        // compute the "equivalence radius" r_0 of the connection
        assert(K[0][0] > 0.0);
        assert(K[1][1] > 0.0);
        Scalar tmp1 = std::sqrt(K[1][1]/K[0][0]);
        Scalar tmp2 = 1.0 / tmp1;
        Scalar r0 = std::sqrt(D[0]*D[0]*tmp1 + D[1]*D[1]*tmp2);
        r0 /= std::sqrt(tmp1) + std::sqrt(tmp2);
        r0 *= 0.28;

        // we assume the aquifer borehole in the center of the dof and that it is vertical,
        // i.e., the area which is exposed to the flow is 2*pi*r0*h. (for non-vertical
        // aquifers this would need to be multiplied with the cosine of the angle and the
        // height must be adapted...)
        const Scalar exposureFactor = 2*M_PI;

        dofVars.connectionTransmissibilityFactor = exposureFactor*Kh/(std::log(r0 / rAquifer) + S);
    }

    template <class ResultEval, class BhpEval>
    void computeVolumetricDofRates_(std::array<ResultEval, numPhases>& volRates,
                                    const BhpEval& bottomHolePressure,
                                    const DofVariables& dofVars) const
    {
        typedef Opm::MathToolbox<Evaluation> DofVarsToolbox;
        typedef typename std::conditional<std::is_same<BhpEval, Scalar>::value,
                                          ResultEval,
                                          Scalar>::type DofEval;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            volRates[phaseIdx] = 0.0;

        // connection transmissibility factor for the current DOF.
        Scalar Twj = dofVars.connectionTransmissibilityFactor;

        // bottom hole pressure and depth of the degree of freedom
        ResultEval pbh = bottomHolePressure;
        Scalar depth = dofVars.depth;

        // gravity constant
        Scalar g = simulator_.problem().gravity()[dimWorld - 1];

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            // aquifer model due to Peaceman; see Chen et al., p. 449

            // phase pressure in grid cell
            const DofEval& p = DofVarsToolbox::template decay<DofEval>(dofVars.pressure[phaseIdx]);

            // density and mobility of fluid phase
            const DofEval& rho = DofVarsToolbox::template decay<DofEval>(dofVars.density[phaseIdx]);
            DofEval lambda;
            if (aquiferType_ == Producer) {
                //assert(p < pbh);
                lambda = DofVarsToolbox::template decay<DofEval>(dofVars.mobility[phaseIdx]);
            }
            else if (aquiferType_ == Injector) {
                //assert(p > pbh);
                if (phaseIdx != injectedPhaseIdx_)
                    continue;

                // use the total mobility, i.e. the sum of all phase mobilities at the
                // injector cell. this seems a bit weird: at the wall of the borehole,
                // there should only be injected phase present, so its mobility should be
                // 1/viscosity...
                lambda = 0.0;
                for (unsigned phase2Idx = 0; phase2Idx < numPhases; ++phase2Idx) {
                    if (!FluidSystem::phaseIsActive(phase2Idx))
                        continue;

                    lambda += DofVarsToolbox::template decay<DofEval>(dofVars.mobility[phase2Idx]);
                }
            }
            else
                OPM_THROW(std::logic_error,
                          "Type of aquifer \"" << name() << "\" is undefined");

            Opm::Valgrind::CheckDefined(pbh);
            Opm::Valgrind::CheckDefined(p);
            Opm::Valgrind::CheckDefined(g);
            Opm::Valgrind::CheckDefined(rho);
            Opm::Valgrind::CheckDefined(lambda);
            Opm::Valgrind::CheckDefined(depth);
            Opm::Valgrind::CheckDefined(refDepth_);

            // pressure in the borehole ("hole pressure") at the given location
            ResultEval ph = pbh + rho*g*(depth - refDepth_);

            // volumetric reservoir rate for the phase
            volRates[phaseIdx] = Twj*lambda*(ph - p);

            Opm::Valgrind::CheckDefined(g);
            Opm::Valgrind::CheckDefined(ph);
            Opm::Valgrind::CheckDefined(volRates[phaseIdx]);
        }
    }

    /*!
     * \brief Given the volumetric rates for all phases, return the
     *        corresponding weighted rate
     *
     * The weights are user-specified and can be set using
     * setVolumetricPhaseWeights()
     */
    template <class Eval>
    Eval computeWeightedRate_(const std::array<Eval, numPhases>& volRates) const
    {
        Eval result = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            result += volRates[phaseIdx]*volumetricWeight_[phaseIdx];
        }
        return result;
    }

    /*!
     * \brief Convert volumetric reservoir rates into volumetric volume rates.
     *
     * This requires the density and composition of the phases and
     * thus the applicable fluid state.
     */
    template <class Eval>
    void computeSurfaceRates_(std::array<Eval, numPhases>& surfaceRates,
                              const std::array<Eval, numPhases>& reservoirRate,
                              const DofVariables& dofVars) const
    {
        // the array for the surface rates and the one for the reservoir rates must not
        // be the same!
        assert(&surfaceRates != &reservoirRate);

        int regionIdx = dofVars.pvtRegionIdx;

        // If your compiler bails out here, you have not chosen the correct fluid
        // system. Currently, only Opm::FluidSystems::BlackOil is supported, sorry...
        Scalar rhoOilSurface = FluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
        Scalar rhoGasSurface = FluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rhoWaterSurface = FluidSystem::referenceDensity(waterPhaseIdx, regionIdx);

        // oil
        if (FluidSystem::phaseIsActive(oilPhaseIdx))
            surfaceRates[oilPhaseIdx] =
                // oil in gas phase
                reservoirRate[gasPhaseIdx]
                * Toolbox::value(dofVars.density[gasPhaseIdx])
                * Toolbox::value(dofVars.gasMassFraction[oilCompIdx])
                / rhoOilSurface
                +
                // oil in oil phase
                reservoirRate[oilPhaseIdx]
                * Toolbox::value(dofVars.density[oilPhaseIdx])
                * Toolbox::value(dofVars.oilMassFraction[oilCompIdx])
                / rhoOilSurface;

        // gas
        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            surfaceRates[gasPhaseIdx] =
                // gas in gas phase
                reservoirRate[gasPhaseIdx]
                * Toolbox::value(dofVars.density[gasPhaseIdx])
                * Toolbox::value(dofVars.gasMassFraction[gasCompIdx])
                / rhoGasSurface
                +
                // gas in oil phase
                reservoirRate[oilPhaseIdx]
                * Toolbox::value(dofVars.density[oilPhaseIdx])
                * Toolbox::value(dofVars.oilMassFraction[gasCompIdx])
                / rhoGasSurface;

        // water
        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            surfaceRates[waterPhaseIdx] =
                reservoirRate[waterPhaseIdx]
                * Toolbox::value(dofVars.density[waterPhaseIdx])
                / rhoWaterSurface;
    }

    /*!
     * \brief Compute the volumetric phase rate of the complete aquifer given a bottom hole
     *        pressure.
     *
     * A single degree of freedom may be different from the evaluation point.
     */
    void computeOverallRates_(Scalar bottomHolePressure,
                              std::array<Scalar, numPhases>& overallResvRates,
                              std::array<Scalar, numPhases>& overallSurfaceRates,
                              const DofVariables *evalDofVars = 0,
                              int globalEvalDofIdx = -1) const

    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            overallResvRates[phaseIdx] = 0.0;
            overallSurfaceRates[phaseIdx] = 0.0;
        }

        auto dofVarsIt = dofVariables_.begin();
        const auto& dofVarsEndIt = dofVariables_.end();
        for (; dofVarsIt != dofVarsEndIt; ++ dofVarsIt) {
            std::array<Scalar, numPhases> volumetricReservoirRates;
            const DofVariables *tmp;
            if (dofVarsIt->first == globalEvalDofIdx)
                tmp = evalDofVars;
            else
                tmp = dofVarsIt->second;

            computeVolumetricDofRates_<Scalar, Scalar>(volumetricReservoirRates, bottomHolePressure, *tmp);

            std::array<Scalar, numPhases> volumetricSurfaceRates;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                volumetricSurfaceRates[ phaseIdx ] = 0;
            }
            computeSurfaceRates_(volumetricSurfaceRates, volumetricReservoirRates, *tmp);

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                overallResvRates[phaseIdx] += volumetricReservoirRates[phaseIdx];
                overallSurfaceRates[phaseIdx] += volumetricSurfaceRates[phaseIdx];
            }
        }
    }

    /*!
     * \brief Compute the weighted volumetric rate of the complete aquifer given a bottom
     *        hole pressure.
     *
     * A single degree of freedom may be different from the evaluation point.
     */
    Scalar computeOverallWeightedSurfaceRate_(Scalar bottomHolePressure,
                                              std::array<Scalar, numPhases>& overallSurfaceRates,
                                              const DofVariables& evalDofVars,
                                              int globalEvalDofIdx) const

    {
        static std::array<Scalar, numPhases> resvRatesDummy;
        computeOverallRates_(bottomHolePressure,
                             overallSurfaceRates,
                             resvRatesDummy,
                             evalDofVars,
                             globalEvalDofIdx);
        return computeWeightedRate_(overallSurfaceRates);
    }

    // this is a more convenient version of the method above if all degrees of freedom
    // are supposed to be at their evaluation points.
    Scalar computeOverallWeightedSurfaceRate_(Scalar bottomHolePressure,
                                              std::array<Scalar, numPhases>& overallSurfaceRates) const
    {
        // create a dummy DofVariables object and call the method above using an index
        // that is guaranteed to never be part of a aquifer...
        static DofVariables dummyDofVars;
        return computeOverallWeightedSurfaceRate_(bottomHolePressure,
                                                  overallSurfaceRates,
                                                  dummyDofVars,
                                                  /*globalEvalDofIdx=*/-1);
    }

    /*!
     * \brief Compute the "rate-equivalent bottom hole pressure"
     *
     * I.e. The bottom hole pressure where the aquifer rate is exactly the one which is
     * targeted. This is zero of the "rate-equivalent bottom hole pressure" would be
     * smaller than 1 bar.
     */
    Scalar computeRateEquivalentBhp_() const
    {
        if (aquiferStatus() == Shut)
            // there is no flow happening in the aquifer, so we return 0...
            return 0.0;

        // initialize the bottom hole pressure which we would like to calculate
        Scalar bhpScalar = actualBottomHolePressure_;
        if (bhpScalar > 1e8)
            bhpScalar = 1e8;
        if (bhpScalar < 1e5)
            bhpScalar = 1e5;

        // if the BHP goes below 1 bar for the first time, we reset it to 10 bars and
        // are "on bail", i.e. if it goes below 1 bar again, we give up because the
        // final pressure would be below 1 bar...
        bool onBail = false;

        // Newton-Raphson method
        typedef Opm::DenseAd::Evaluation<Scalar, 1> BhpEval;

        BhpEval bhpEval(bhpScalar);
        bhpEval.setDerivative(0, 1.0);
        const Scalar tolerance = 1e3*std::numeric_limits<Scalar>::epsilon();
        for (int iterNum = 0; iterNum < 20; ++iterNum) {
            const auto& f = aquiferResidual_<BhpEval>(bhpEval);

            if (std::abs(f.derivative(0)) < 1e-20)
                OPM_THROW(Opm::NumericalProblem,
                          "Cannot determine the bottom hole pressure for aquifer " << name()
                          << ": Derivative of the aquifer residual is too small");
            Scalar delta = f.value()/f.derivative(0);

            bhpEval.setValue(bhpEval.value() - delta);
            if (bhpEval < 1e5) {
                bhpEval.setValue(1e5);
                if (onBail)
                    return bhpEval.value();
                else
                    onBail = true;
            }
            else
                onBail = false;

            if (std::abs(delta/bhpEval.value()) < tolerance)
                return bhpEval.value();
        }

        OPM_THROW(Opm::NumericalProblem,
                  "Could not determine the bottom hole pressure of aquifer '" << name()
                  << "' within 20 iterations.");
    }

    template <class BhpEval>
    BhpEval aquiferResidual_(const BhpEval& bhp,
                          const DofVariables *replacementDofVars = 0,
                          int replacedGridIdx = -1) const
    {
        typedef Opm::MathToolbox<BhpEval> BhpEvalToolbox;

        // compute the volumetric reservoir and surface rates for the complete aquifer
        BhpEval resvRate = 0.0;

        std::array<BhpEval, numPhases> totalSurfaceRates;
        std::fill(totalSurfaceRates.begin(), totalSurfaceRates.end(), 0.0);

        auto dofVarsIt = dofVariables_.begin();
        const auto& dofVarsEndIt = dofVariables_.end();
        for (; dofVarsIt != dofVarsEndIt; ++ dofVarsIt) {
            std::array<BhpEval, numPhases> resvRates;
            const DofVariables *dofVars = dofVarsIt->second;
            if (replacedGridIdx == dofVarsIt->first)
                dofVars = replacementDofVars;
            computeVolumetricDofRates_(resvRates, bhp, *dofVars);

            std::array<BhpEval, numPhases> surfaceRates;
            computeSurfaceRates_(surfaceRates, resvRates, *dofVars);

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                totalSurfaceRates[phaseIdx] += surfaceRates[phaseIdx];
            }

            resvRate += computeWeightedRate_(resvRates);
        }

        BhpEval surfaceRate = computeWeightedRate_(totalSurfaceRates);

        // compute the residual of aquifer equation. we currently use max(rateMax - rate,
        // bhp - targetBhp) for producers and max(rateMax - rate, bhp - targetBhp) for
        // injectors. (i.e., the target bottom hole pressure is an upper limit for
        // injectors and a lower limit for producers.) Note that with this approach, one
        // of the limits must always be reached to get the aquifer equation to zero...
        Opm::Valgrind::CheckDefined(maximumSurfaceRate_);
        Opm::Valgrind::CheckDefined(maximumReservoirRate_);
        Opm::Valgrind::CheckDefined(surfaceRate);
        Opm::Valgrind::CheckDefined(resvRate);

        BhpEval result = 1e30;

        BhpEval maxSurfaceRate = maximumSurfaceRate_;
        BhpEval maxResvRate = maximumReservoirRate_;
        if (aquiferStatus() == Closed) {
            // make the weight of the fluids on the surface equal and require that no
            // fluids are produced on the surface...
            maxSurfaceRate = 0.0;
            surfaceRate = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                surfaceRate += totalSurfaceRates[phaseIdx];
            }

            // don't care about the reservoir rate...
            maxResvRate = 1e30;
        }

        if (aquiferType_ == Injector) {
            // for injectors the computed rates are positive and the target BHP is the
            // maximum allowed pressure ...
            result = BhpEvalToolbox::min(maxSurfaceRate - surfaceRate, result);
            result = BhpEvalToolbox::min(maxResvRate - resvRate, result);
            result = BhpEvalToolbox::min(1e-7*(targetBottomHolePressure_ - bhp), result);
        }
        else {
            assert(aquiferType_ == Producer);
            // ... for producers the rates are negative and the bottom hole pressure is
            // is the minimum
            result = BhpEvalToolbox::min(maxSurfaceRate + surfaceRate, result);
            result = BhpEvalToolbox::min(maxResvRate + resvRate, result);
            result = BhpEvalToolbox::min(1e-7*(bhp - targetBottomHolePressure_), result);
        }

        const Scalar scalingFactor = 1e-3;
        return scalingFactor*result;
    }

    const Simulator& simulator_;

    std::string name_;

    std::vector<DofVariables, Ewoms::aligned_allocator<DofVariables, alignof(DofVariables)> > dofVarsStore_;
    std::map<int, DofVariables*> dofVariables_;

    // the number of times beginIteration*() was called for the current time step
    unsigned iterationIdx_;

    // the type of the aquifer (injector, producer or undefined)
    AquiferType aquiferType_;

    // Specifies whether the aquifer is currently open, closed or shut. The difference
    // between "closed" and "shut" is that for the former, the aquifer is assumed to be
    // closed above the reservoir so that cross-flow within the aquifer is possible while
    // the aquifer is completely separated from the reservoir if it is shut. (i.e., no
    // crossflow is possible in this case.)
    AquiferStatus aquiferStatus_;

    // specifies the quantities which are controlled for (i.e., which
    // should be assumed to be externally specified and which should
    // be computed based on those)
    ControlMode controlMode_;

    // the sum of the total volumes of all the degrees of freedoms that interact with the aquifer
    Scalar aquiferTotalVolume_;

    // The assumed bottom hole and tubing head pressures as specified by the user
    Scalar bhpLimit_;
    Scalar thpLimit_;

    // The bottom hole pressure to be targeted by the aquifer model. This may be computed
    // from the tubing head pressure (if the control mode is TubingHeadPressure), or it may be
    // just the user-specified bottom hole pressure if the control mode is
    // BottomHolePressure.
    Scalar targetBottomHolePressure_;

    // The bottom hole pressure which is actually observed in the aquifer
    Scalar actualBottomHolePressure_;

    // The maximum weighted volumetric surface rates specified by the
    // user. This is used to apply rate limits and it is to be read as
    // the maximum absolute value of the rate, i.e., the aquifer can
    // produce or inject the given amount.
    Scalar maximumSurfaceRate_;

    // The maximum weighted volumetric reservoir rates specified by
    // the user. This is used to apply rate limits and it is to be
    // read as the maximum absolute value of the rate, i.e., the aquifer
    // can produce or inject the given amount.
    Scalar maximumReservoirRate_;

    // The volumetric surface rate which is actually observed in the aquifer
    Scalar actualWeightedSurfaceRate_;
    std::array<Scalar, numPhases> actualSurfaceRates_;

    // The volumetric reservoir rate which is actually observed in the aquifer
    Scalar actualWeightedResvRate_;
    std::array<Scalar, numPhases> actualResvRates_;

    // The relative weight of the volumetric rate of each fluid
    Scalar volumetricWeight_[numPhases];

    // the reference depth for the bottom hole pressure. if not specified otherwise, this
    // is the position of the _highest_ DOF in the aquifer.
    Scalar refDepth_;

    // The thermodynamic state of the fluid which gets injected
    //
    // The fact that this attribute is mutable is kind of an hack
    // which can be avoided using a PressureOverlayFluidState, but
    // then performance would be slightly worse...
    mutable FluidState injectionFluidState_;

    unsigned injectedPhaseIdx_;
};
} // namespace Ewoms

#endif
