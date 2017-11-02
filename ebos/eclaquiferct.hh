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
        
    }

    /*!
     * \copydoc Ewoms::BaseAuxiliaryModule::addNeighbors()
     */
    virtual void applyInitial()
    {
        
    }

    /*!
     * \copydoc Ewoms::BaseAuxiliaryModule::linearize()
     */
    virtual void linearize(JacobianMatrix& matrix, GlobalEqVector& residual)
    {
        
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
        
    }

    /*!
     * \brief Finalize the specification of the borehole.
     */
    void endSpec()
    {
        
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
        return 0.;
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
        
    }

    /*!
     * \brief Informs the aquifer that an iteration has just begun.
     *
     * This is the post-processing part which uses the results of the
     * accumulation callback.
     */
    void beginIterationPostProcess()
    {
        
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

            modelRate.setVolumetricRate(intQuants.fluidState(), phaseIdx, 0.);
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
    }

    template <class ResultEval, class BhpEval>
    void computeVolumetricDofRates_(std::array<ResultEval, numPhases>& volRates,
                                    const BhpEval& bottomHolePressure,
                                    const DofVariables& dofVars) const
    {
        
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
        return 0.;
    }

    // this is a more convenient version of the method above if all degrees of freedom
    // are supposed to be at their evaluation points.
    Scalar computeOverallWeightedSurfaceRate_(Scalar bottomHolePressure,
                                              std::array<Scalar, numPhases>& overallSurfaceRates) const
    {
        return 0.;
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
        return 0.;
    }

    template <class BhpEval>
    BhpEval aquiferResidual_(const BhpEval& bhp,
                          const DofVariables *replacementDofVars = 0,
                          int replacedGridIdx = -1) const
    {
        typedef Opm::MathToolbox<BhpEval> BhpEvalToolbox;

        // compute the volumetric reservoir and surface rates for the complete aquifer
        BhpEval resvRate = 0.0;
        return resvRate;
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
