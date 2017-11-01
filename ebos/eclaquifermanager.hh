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
 * \copydoc Ewoms::EclAquiferManager
 */
#ifndef EWOMS_ECL_AQUIFER_MANAGER_HH
#define EWOMS_ECL_AQUIFER_MANAGER_HH

#include "eclaquiferct.hh"

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>
//#include <opm/parser/eclipse/EclipseState/Schedule/Aquifer.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/parallel/threadedentityiterator.hh>

#include <dune/grid/common/gridenums.hh>

#include <map>
#include <string>
#include <vector>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Grid);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief A class which handles aquifer controls as specified by an
 *        Eclipse deck
 */
template <class TypeTag>
class EclAquiferManager
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };

    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Ewoms::EclCarterTracyAquifer<TypeTag> Aquifer;

    typedef std::map<int, std::pair<const Opm::Completion*, std::shared_ptr<Aquifer> > > AquiferCompletionsMap;

    typedef Dune::FieldVector<Evaluation, numEq> EvalEqVector;

public:
    EclAquiferManager(Simulator& simulator)
        : simulator_(simulator)
    { }

    /*!
     * \brief This sets up the basic properties of all aquifers.
     *
     * I.e., aquifer positions, names etc...
     */
    void init(const Opm::EclipseState& eclState)
    {
        const auto& deckSchedule = eclState.getSchedule();

        // create the aquifers which intersect with the current process' grid
        for (size_t deckAquiferIdx = 0; deckAquiferIdx < deckSchedule.numAquifers(); ++deckAquiferIdx)
        {
            const Opm::Aquifer* deckAquifer = deckSchedule.getAquifers()[deckAquiferIdx];
            const std::string& aquiferName = deckAquifer->name();

            // set the name of the aquifer but not much else. (i.e., if it is not completed,
            // the aquifer primarily serves as a placeholder.) The big rest of the aquifer is
            // specified by the updateAquiferCompletions_() method
            auto aquifer = std::make_shared<Aquifer>(simulator_);
            aquifer->setName(aquiferName);
            aquifer->setAquiferStatus(Aquifer::Shut);

            aquifers_.push_back(aquifer);
            aquiferNameToIndex_[aquifer->name()] = aquifers_.size() - 1;
        }
    }

    /*!
     * \brief This should be called the problem before each simulation
     *        episode to adapt the aquifer controls.
     */
    void beginEpisode(const Opm::EclipseState& eclState, bool wasRestarted=false)
    {
        unsigned episodeIdx = simulator_.episodeIndex();

        const auto& deckSchedule = eclState.getSchedule();
        AquiferCompletionsMap aquiferCompMap;
        computeAquiferCompletionsMap_(episodeIdx, aquiferCompMap);

        if (wasRestarted || aquiferTopologyChanged_(eclState, episodeIdx))
            updateAquiferTopology_(episodeIdx, aquiferCompMap, gridDofIsPenetrated_);

        // set those parameters of the aquifers which do not change the topology of the
        // linearized system of equations
        updateAquiferParameters_(episodeIdx, aquiferCompMap);

        const std::vector<const Opm::Aquifer*>& deckAquifers = deckSchedule.getAquifers(episodeIdx);
        // set the injection data for the respective aquifers.
        for (size_t deckAquiferIdx = 0; deckAquiferIdx < deckAquifers.size(); ++deckAquiferIdx) {
            const Opm::Aquifer* deckAquifer = deckAquifers[deckAquiferIdx];

            if (!hasAquifer(deckAquifer->name()))
                continue;

            auto aquifer = this->aquifer(deckAquifer->name());

            Opm::AquiferCommon::StatusEnum deckAquiferStatus = deckAquifer->getStatus(episodeIdx);
            switch (deckAquiferStatus) {
            case Opm::AquiferCommon::AUTO:
                // TODO: for now, auto means open...
            case Opm::AquiferCommon::OPEN:
                aquifer->setAquiferStatus(Aquifer::Open);
                break;
            case Opm::AquiferCommon::STOP:
                aquifer->setAquiferStatus(Aquifer::Closed);
                break;
            case Opm::AquiferCommon::SHUT:
                aquifer->setAquiferStatus(Aquifer::Shut);
                break;
            }

            // make sure that the aquifer is either an injector or a
            // producer for the current episode. (it is not allowed to
            // be neither or to be both...)
            assert((deckAquifer->isInjector(episodeIdx)?1:0) +
                   (deckAquifer->isProducer(episodeIdx)?1:0) == 1);

            if (deckAquifer->isInjector(episodeIdx)) {
                aquifer->setAquiferType(Aquifer::Injector);

                const Opm::AquiferInjectionProperties& injectProperties =
                    deckAquifer->getInjectionProperties(episodeIdx);

                switch (injectProperties.injectorType) {
                case Opm::AquiferInjector::WATER:
                    aquifer->setInjectedPhaseIndex(FluidSystem::waterPhaseIdx);
                    break;
                case Opm::AquiferInjector::GAS:
                    aquifer->setInjectedPhaseIndex(FluidSystem::gasPhaseIdx);
                    break;
                case Opm::AquiferInjector::OIL:
                    aquifer->setInjectedPhaseIndex(FluidSystem::oilPhaseIdx);
                    break;
                case Opm::AquiferInjector::MULTI:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Multi-phase injector aquifers");
                }

                switch (injectProperties.controlMode) {
                case Opm::AquiferInjector::RATE:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricSurfaceRate);
                    break;

                case Opm::AquiferInjector::RESV:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricReservoirRate);
                    break;

                case Opm::AquiferInjector::BHP:
                    aquifer->setControlMode(Aquifer::ControlMode::BottomHolePressure);
                    break;

                case Opm::AquiferInjector::THP:
                    aquifer->setControlMode(Aquifer::ControlMode::TubingHeadPressure);
                    break;

                case Opm::AquiferInjector::GRUP:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Aquifer groups");

                case Opm::AquiferInjector::CMODE_UNDEFINED:
                    std::cout << "Warning: Control mode of injection aquifer " << aquifer->name()
                              << " is undefined. Assuming aquifer to be shut.\n";
                    aquifer->setAquiferStatus(Aquifer::AquiferStatus::Shut);
                    continue;
                }

                switch (injectProperties.injectorType) {
                case Opm::AquiferInjector::WATER:
                    aquifer->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/0.0, /*water=*/1.0);
                    break;

                case Opm::AquiferInjector::OIL:
                    aquifer->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/0.0);
                    break;

                case Opm::AquiferInjector::GAS:
                    aquifer->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/1.0, /*water=*/0.0);
                    break;

                case Opm::AquiferInjector::MULTI:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Multi-phase injection aquifers");
                }

                aquifer->setMaximumSurfaceRate(injectProperties.surfaceInjectionRate);
                aquifer->setMaximumReservoirRate(injectProperties.reservoirInjectionRate);
                aquifer->setTargetBottomHolePressure(injectProperties.BHPLimit);

                // TODO
                aquifer->setTargetTubingHeadPressure(1e30);
                //aquifer->setTargetTubingHeadPressure(injectProperties.THPLimit);
            }

            if (deckAquifer->isProducer(episodeIdx)) {
                aquifer->setAquiferType(Aquifer::Producer);

                const Opm::AquiferProductionProperties& producerProperties =
                    deckAquifer->getProductionProperties(episodeIdx);

                switch (producerProperties.controlMode) {
                case Opm::AquiferProducer::ORAT:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricSurfaceRate);
                    aquifer->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/0.0);
                    aquifer->setMaximumSurfaceRate(producerProperties.OilRate);
                    break;

                case Opm::AquiferProducer::GRAT:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricSurfaceRate);
                    aquifer->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/1.0, /*water=*/0.0);
                    aquifer->setMaximumSurfaceRate(producerProperties.GasRate);
                    break;

                case Opm::AquiferProducer::WRAT:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricSurfaceRate);
                    aquifer->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/0.0, /*water=*/1.0);
                    aquifer->setMaximumSurfaceRate(producerProperties.WaterRate);
                    break;

                case Opm::AquiferProducer::LRAT:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricSurfaceRate);
                    aquifer->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/1.0);
                    aquifer->setMaximumSurfaceRate(producerProperties.LiquidRate);
                    break;

                case Opm::AquiferProducer::CRAT:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Linearly combined rates");

                case Opm::AquiferProducer::RESV:
                    aquifer->setControlMode(Aquifer::ControlMode::VolumetricReservoirRate);
                    aquifer->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/1.0, /*water=*/1.0);
                    aquifer->setMaximumSurfaceRate(producerProperties.ResVRate);
                    break;

                case Opm::AquiferProducer::BHP:
                    aquifer->setControlMode(Aquifer::ControlMode::BottomHolePressure);
                    break;

                case Opm::AquiferProducer::THP:
                    aquifer->setControlMode(Aquifer::ControlMode::TubingHeadPressure);
                    break;

                case Opm::AquiferProducer::GRUP:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Aquifer groups");

                case Opm::AquiferProducer::NONE:
                    // fall-through
                case Opm::AquiferProducer::CMODE_UNDEFINED:
                    std::cout << "Warning: Control mode of production aquifer " << aquifer->name()
                              << " is undefined. Assuming aquifer to be shut.";
                    aquifer->setAquiferStatus(Aquifer::AquiferStatus::Shut);
                    continue;
                }

                aquifer->setTargetBottomHolePressure(producerProperties.BHPLimit);

                // TODO
                aquifer->setTargetTubingHeadPressure(-1e30);
                //aquifer->setTargetTubingHeadPressure(producerProperties.THPLimit);
            }
        }
    }

    /*!
     * \brief Return the number of aquifers considered by the EclAquiferManager.
     */
    unsigned numAquifers() const
    { return aquifers_.size(); }

    /*!
     * \brief Return if a given aquifer name is known to the aquifers manager
     */
    bool hasAquifer(const std::string& aquiferName) const
    {
        return aquiferNameToIndex_.find( aquiferName ) != aquiferNameToIndex_.end();
    }

    /*!
     * \brief Returns true iff a given degree of freedom is currently penetrated by any aquifer.
     */
    bool gridDofIsPenetrated(unsigned globalDofIdx) const
    { return gridDofIsPenetrated_[globalDofIdx]; }

    /*!
     * \brief Given a aquifer name, return the corresponding index.
     *
     * A std::runtime_error will be thrown if the aquifer name is unknown.
     */
    unsigned aquiferIndex(const std::string& aquiferName) const
    {
        assert( hasAquifer( aquiferName ) );
        const auto& it = aquiferNameToIndex_.find(aquiferName);
        if (it == aquiferNameToIndex_.end())
        {
            OPM_THROW(std::runtime_error,
                      "No aquifer called '" << aquiferName << "'found");
        }
        return it->second;
    }

    /*!
     * \brief Given a aquifer name, return the corresponding aquifer.
     *
     * A std::runtime_error will be thrown if the aquifer name is unknown.
     */
    std::shared_ptr<const Aquifer> aquifer(const std::string& aquiferName) const
    { return aquifers_[aquiferIndex(aquiferName)]; }

    /*!
     * \brief Given a aquifer name, return the corresponding aquifer.
     *
     * A std::runtime_error will be thrown if the aquifer name is unknown.
     */
    std::shared_ptr<Aquifer> aquifer(const std::string& aquiferName)
    { return aquifers_[aquiferIndex(aquiferName)]; }

    /*!
     * \brief Given a aquifer index, return the corresponding aquifer.
     */
    std::shared_ptr<const Aquifer> aquifer(size_t aquiferIdx) const
    { return aquifers_[aquiferIdx]; }

    /*!
     * \brief Given a aquifer index, return the corresponding aquifer.
     */
    std::shared_ptr<Aquifer> aquifer(size_t aquiferIdx)
    { return aquifers_[aquiferIdx]; }

    /*!
     * \brief Informs the aquifer manager that a time step has just begun.
     */
    void beginTimeStep()
    {
        // iterate over all aquifers and notify them individually
        for (size_t aquiferIdx = 0; aquiferIdx < aquifers_.size(); ++aquiferIdx)
            aquifers_[aquiferIdx]->beginTimeStep();
    }

    /*!
     * \brief Informs the aquifer that an iteration has just begun.
     *
     * In this method, the aquifer calculates the bottom hole and tubing head pressures, the
     * actual unconstraint production and injection rates, etc.
     */
    void beginIteration()
    {
        // call the preprocessing routines
        const size_t aquiferSize = aquifers_.size();
        for (size_t aquiferIdx = 0; aquiferIdx < aquiferSize; ++aquiferIdx)
            aquifers_[aquiferIdx]->beginIterationPreProcess();

        // call the accumulation routines
        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(simulator_.gridManager().gridView());
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            ElementContext elemCtx(simulator_);
            auto elemIt = threadedElemIt.beginParallel();
            for (; !threadedElemIt.isFinished(elemIt); elemIt = threadedElemIt.increment()) {
                const Element& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                for (size_t aquiferIdx = 0; aquiferIdx < aquiferSize; ++aquiferIdx)
                    aquifers_[aquiferIdx]->beginIterationAccumulate(elemCtx, /*timeIdx=*/0);
            }
        }

        // call the postprocessing routines
        for (size_t aquiferIdx = 0; aquiferIdx < aquiferSize; ++aquiferIdx)
            aquifers_[aquiferIdx]->beginIterationPostProcess();
    }

    /*!
     * \brief Informs the aquifer manager that an iteration has just been finished.
     */
    void endIteration()
    {
        // iterate over all aquifers and notify them individually
        const size_t aquiferSize = aquifers_.size();
        for (size_t aquiferIdx = 0; aquiferIdx < aquiferSize; ++aquiferIdx)
            aquifers_[aquiferIdx]->endIteration();
    }

    /*!
     * \brief Informs the aquifer manager that a time step has just been finished.
     */
    void endTimeStep()
    {
        Scalar dt = simulator_.timeStepSize();

        // iterate over all aquifers and notify them individually. also, update the
        // production/injection totals for the active aquifers.
        const size_t aquiferSize = aquifers_.size();
        for (size_t aquiferIdx = 0; aquiferIdx < aquiferSize; ++aquiferIdx) {
            auto aquifer = aquifers_[aquiferIdx];
            aquifer->endTimeStep();

            // update the surface volumes of the produced/injected fluids
            std::array<Scalar, numPhases>* injectedVolume;
            if (aquiferTotalInjectedVolume_.count(aquifer->name()) == 0) {
                injectedVolume = &aquiferTotalInjectedVolume_[aquifer->name()];
                std::fill(injectedVolume->begin(), injectedVolume->end(), 0.0);
            }
            else
                injectedVolume = &aquiferTotalInjectedVolume_[aquifer->name()];

            std::array<Scalar, numPhases>* producedVolume;
            if (aquiferTotalProducedVolume_.count(aquifer->name()) == 0) {
                producedVolume = &aquiferTotalProducedVolume_[aquifer->name()];
                std::fill(producedVolume->begin(), producedVolume->end(), 0.0);
            }
            else
                producedVolume = &aquiferTotalProducedVolume_[aquifer->name()];

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // this assumes that the implicit Euler method is used for time
                // integration. TODO: Once the time discretization becomes pluggable,
                // this integration needs to be done by the time discretization code!
                Scalar vol = dt * aquifer->surfaceRate(phaseIdx);

                if (vol < 0)
                    (*producedVolume)[phaseIdx] += -vol;
                else
                    (*injectedVolume)[phaseIdx] += vol;
            }
        }
    }

    /*!
     * \brief Informs the aquifer manager that a simulation episode has just been finished.
     */
    void endEpisode()
    { }

    /*!
     * \brief Returns the surface volume of a fluid phase produced by a aquifer.
     */
    Scalar totalProducedVolume(const std::string& aquiferName, unsigned phaseIdx) const
    {
        if (aquiferTotalProducedVolume_.count(aquiferName) == 0)
            return 0.0; // aquifer not yet seen
        return aquiferTotalProducedVolume_.at(aquiferName)[phaseIdx];
    }

    /*!
     * \brief Returns the surface volume of a fluid phase injected by a aquifer.
     */
    Scalar totalInjectedVolume(const std::string& aquiferName, unsigned phaseIdx) const
    {
        if (aquiferTotalInjectedVolume_.count(aquiferName) == 0)
            return 0.0; // aquifer not yet seen
        return aquiferTotalInjectedVolume_.at(aquiferName)[phaseIdx];
    }

    /*!
     * \brief Computes the source term due to aquifers for a degree of
     *        freedom.
     */
    template <class Context>
    void computeTotalRatesForDof(EvalEqVector& q,
                                 const Context& context,
                                 unsigned dofIdx,
                                 unsigned timeIdx) const
    {
        q = 0.0;

        if (!gridDofIsPenetrated(context.globalSpaceIndex(dofIdx, timeIdx)))
            return;

        RateVector aquiferRate;

        // iterate over all aquifers and add up their individual rates
        const size_t aquiferSize = aquifers_.size();
        for (size_t aquiferIdx = 0; aquiferIdx < aquiferSize; ++aquiferIdx) {
            aquiferRate = 0.0;
            aquifers_[aquiferIdx]->computeTotalRatesForDof(aquiferRate, context, dofIdx, timeIdx);
            for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
                q[eqIdx] += aquiferRate[eqIdx];
        }
    }

    /*!
     * \brief This method writes the complete state of all aquifers
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter& res OPM_UNUSED)
    {
        /* do nothing: Everything which we need here is provided by the deck->.. */
    }

    /*!
     * \brief This method restores the complete state of the all aquifers
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter& res OPM_UNUSED)
    {
        // initialize the aquifers for the current episode
        beginEpisode(simulator_.gridManager().eclState(), /*wasRestarted=*/true);
    }

    /*!
     * \brief Returns true if something in a aquifer changed compared to the previous report
     *        step.
     *
     * "Something" can either be the aquifer topology (i.e., which grid blocks are contained
     * in which aquifer) or it can be a aquifer parameter like the bottom hole pressure...
     */
    bool aquifersChanged(const Opm::EclipseState& eclState, unsigned reportStepIdx) const
    {
        if (aquiferTopologyChanged_(eclState, reportStepIdx))
            return true;

        const auto& schedule = eclState.getSchedule();
        if (schedule.getTimeMap().numTimesteps() <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the aquifers don't change
            return false;

        const Opm::Events& events = schedule.getEvents();
        return events.hasEvent(Opm::ScheduleEvents::PRODUCTION_UPDATE |
                               Opm::ScheduleEvents::INJECTION_UPDATE |
                               Opm::ScheduleEvents::WELL_STATUS_CHANGE,
                               reportStepIdx);
    }

protected:
    bool aquiferTopologyChanged_(const Opm::EclipseState& eclState, unsigned reportStepIdx) const
    {
        if (reportStepIdx == 0) {
            // the aquifer topology has always been changed relative to before the
            // simulation is started...
            return true;
        }

        const auto& schedule = eclState.getSchedule();
        if (schedule.getTimeMap().numTimesteps() <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the aquifers don't change
            return false;

        const Opm::Events& events = schedule.getEvents();
        return events.hasEvent(Opm::ScheduleEvents::NEW_WELL |
                               Opm::ScheduleEvents::COMPLETION_CHANGE,
                               reportStepIdx);
    }

    void updateAquiferTopology_(unsigned reportStepIdx OPM_UNUSED,
                             const AquiferCompletionsMap& aquiferCompletions,
                             std::vector<bool>& gridDofIsPenetrated) const
    {
        auto& model = simulator_.model();
        const auto& gridManager = simulator_.gridManager();

        // first, remove all aquifers from the reservoir
        model.clearAuxiliaryModules();
        auto aquiferIt = aquifers_.begin();
        const auto& aquiferEndIt = aquifers_.end();
        for (; aquiferIt != aquiferEndIt; ++aquiferIt) {
            (*aquiferIt)->clear();
            (*aquiferIt)->beginSpec();
        }

        //////
        // tell the active aquifers which DOFs they contain
        const auto gridView = simulator_.gridManager().gridView();

        gridDofIsPenetrated.resize(model.numGridDof());
        std::fill(gridDofIsPenetrated.begin(), gridDofIsPenetrated.end(), false);

        ElementContext elemCtx(simulator_);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto elemEndIt = gridView.template end</*codim=*/0>();
        std::set<std::shared_ptr<Aquifer> > aquifers;
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue; // non-local entities need to be skipped

            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                unsigned cartesianDofIdx = gridManager.cartesianIndex(globalDofIdx);

                if (aquiferCompletions.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any aquifer, so we must skip
                    // it...
                    continue;

                gridDofIsPenetrated[globalDofIdx] = true;

                auto eclAquifer = aquiferCompletions.at(cartesianDofIdx).second;
                eclAquifer->addDof(elemCtx, dofIdx);

                aquifers.insert(eclAquifer);
            }
            //////
        }

        // register all aquifers at the model as auxiliary equations
        aquiferIt = aquifers_.begin();
        for (; aquiferIt != aquiferEndIt; ++aquiferIt) {
            (*aquiferIt)->endSpec();
            model.addAuxiliaryModule(*aquiferIt);
        }
    }

    void computeAquiferCompletionsMap_(unsigned reportStepIdx OPM_UNUSED, AquiferCompletionsMap& cartesianIdxToCompletionMap)
    {
        const auto& eclState = simulator_.gridManager().eclState();
        const auto& deckSchedule = eclState.getSchedule();

#ifndef NDEBUG
        const auto& eclGrid = eclState.getInputGrid();
        assert( int(eclGrid.getNX()) == simulator_.gridManager().cartesianDimensions()[ 0 ] );
        assert( int(eclGrid.getNY()) == simulator_.gridManager().cartesianDimensions()[ 1 ] );
        assert( int(eclGrid.getNZ()) == simulator_.gridManager().cartesianDimensions()[ 2 ] );
#endif

        // compute the mapping from logically Cartesian indices to the aquifer the
        // respective completion.
        const std::vector<const Opm::Aquifer*>& deckAquifers = deckSchedule.getAquifers(reportStepIdx);
        for (size_t deckAquiferIdx = 0; deckAquiferIdx < deckAquifers.size(); ++deckAquiferIdx) {
            const Opm::Aquifer* deckAquifer = deckAquifers[deckAquiferIdx];
            const std::string& aquiferName = deckAquifer->name();

            if (!hasAquifer(aquiferName))
            {
#ifndef NDEBUG
                if( simulator_.gridManager().grid().comm().size() == 1 )
                {
                    std::cout << "Aquifer '" << aquiferName << "' suddenly appears in the completions "
                              << "for the report step, but has not been previously specified. "
                              << "Ignoring.\n";
                }
#endif
                continue;
            }

            std::array<int, 3> cartesianCoordinate;
            // set the aquifer parameters defined by the current set of completions
            const auto& completionSet = deckAquifer->getCompletions(reportStepIdx);
            for (size_t complIdx = 0; complIdx < completionSet.size(); complIdx ++) {
                const auto& completion = completionSet.get(complIdx);
                cartesianCoordinate[ 0 ] = completion.getI();
                cartesianCoordinate[ 1 ] = completion.getJ();
                cartesianCoordinate[ 2 ] = completion.getK();
                unsigned cartIdx = simulator_.gridManager().cartesianIndex( cartesianCoordinate );

                // in this code we only support each cell to be part of at most a single
                // aquifer. TODO (?) change this?
                assert(cartesianIdxToCompletionMap.count(cartIdx) == 0);

                auto eclAquifer = aquifers_[aquiferIndex(aquiferName)];
                cartesianIdxToCompletionMap[cartIdx] = std::make_pair(&completion, eclAquifer);
            }
        }
    }

    void updateAquiferParameters_(unsigned reportStepIdx, const AquiferCompletionsMap& aquiferCompletions)
    {
        const auto& eclState = simulator_.gridManager().eclState();
        const auto& deckSchedule = eclState.getSchedule();
        const std::vector<const Opm::Aquifer*>& deckAquifers = deckSchedule.getAquifers(reportStepIdx);

        // set the reference depth for all aquifers
        for (size_t deckAquiferIdx = 0; deckAquiferIdx < deckAquifers.size(); ++deckAquiferIdx) {
            const Opm::Aquifer* deckAquifer = deckAquifers[deckAquiferIdx];
            const std::string& aquiferName = deckAquifer->name();

            if( hasAquifer( aquiferName ) )
            {
                aquifers_[aquiferIndex(aquiferName)]->clear();
                aquifers_[aquiferIndex(aquiferName)]->setReferenceDepth(deckAquifer->getRefDepth());
            }
        }

        // associate the aquifer completions with grid cells and register them in the
        // Peaceman aquifer object
        const auto& gridManager = simulator_.gridManager();
        const GridView gridView = gridManager.gridView();

        ElementContext elemCtx(simulator_);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto elemEndIt = gridView.template end</*codim=*/0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue; // non-local entities need to be skipped

            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx)
            {
                assert( elemCtx.numPrimaryDof(/*timeIdx=*/0) == 1 );
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                unsigned cartesianDofIdx = gridManager.cartesianIndex(globalDofIdx);

                if (aquiferCompletions.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any aquifer, so we must skip
                    // it...
                    continue;

                const auto& compInfo = aquiferCompletions.at(cartesianDofIdx);
                const Opm::Completion* completion = compInfo.first;
                std::shared_ptr<Aquifer> eclAquifer = compInfo.second;
                eclAquifer->addDof(elemCtx, dofIdx);

                // the catch is a hack for a ideosyncrasy of opm-parser with regard to
                // defaults handling: if the deck did not specify a radius for the
                // completion, there seems to be no other way to detect this except for
                // catching the exception
                try {
                    eclAquifer->setRadius(elemCtx, dofIdx, 0.5*completion->getDiameter());
                }
                catch (const std::logic_error& e)
                {}

                // overwrite the automatically computed effective
                // permeability by the one specified in the deck-> Note: this
                // is not implemented by opm-parser yet...
                /*
                  Scalar Kh = completion->getEffectivePermeability();
                  if (std::isfinite(Kh) && Kh > 0.0)
                      eclAquifer->setEffectivePermeability(elemCtx, dofIdx, Kh);
                */

                // overwrite the automatically computed connection
                // transmissibilty factor by the one specified in the deck->
                const auto& ctf = completion->getConnectionTransmissibilityFactorAsValueObject();
                if (ctf.hasValue() && ctf.getValue() > 0.0)
                    eclAquifer->setConnectionTransmissibilityFactor(elemCtx, dofIdx, ctf.getValue());
            }
        }
    }

    Simulator& simulator_;

    std::vector<std::shared_ptr<Aquifer> > aquifers_;
    std::vector<bool> gridDofIsPenetrated_;
    std::map<std::string, int> aquiferNameToIndex_;
    std::map<std::string, std::array<Scalar, numPhases> > aquiferTotalInjectedVolume_;
    std::map<std::string, std::array<Scalar, numPhases> > aquiferTotalProducedVolume_;
};
} // namespace Ewoms

#endif
