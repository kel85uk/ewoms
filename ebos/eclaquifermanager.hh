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
        
    }

    /*!
     * \brief This should be called the problem before each simulation
     *        episode to adapt the aquifer controls.
     */
    void beginEpisode(const Opm::EclipseState& eclState, bool wasRestarted=false)
    {
        
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
        
    }

    /*!
     * \brief Informs the aquifer manager that an iteration has just been finished.
     */
    void endIteration()
    {
        
    }

    /*!
     * \brief Informs the aquifer manager that a time step has just been finished.
     */
    void endTimeStep()
    {
        
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
        return 0.;
    }

    /*!
     * \brief Returns the surface volume of a fluid phase injected by a aquifer.
     */
    Scalar totalInjectedVolume(const std::string& aquiferName, unsigned phaseIdx) const
    {
        return 0.;
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
        
    }

protected:
    bool aquiferTopologyChanged_(const Opm::EclipseState& eclState, unsigned reportStepIdx) const
    {
        
    }

    void updateAquiferTopology_(unsigned reportStepIdx OPM_UNUSED,
                             const AquiferCompletionsMap& aquiferCompletions,
                             std::vector<bool>& gridDofIsPenetrated) const
    {
        
    }

    void computeAquiferCompletionsMap_(unsigned reportStepIdx OPM_UNUSED, AquiferCompletionsMap& cartesianIdxToCompletionMap)
    {
        
    }

    void updateAquiferParameters_(unsigned reportStepIdx, const AquiferCompletionsMap& aquiferCompletions)
    {
        
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
