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
 *
 * \copydoc Ewoms::WaterCH4Problem
 * This problem adds from bench1.hh the salt component
 */
#ifndef EWOMS_WATER_CH4_PROBLEM_HH
#define EWOMS_WATER_CH4_PROBLEM_HH

#include <ewoms/models/pvs/pvsmodel.hh>
#include <ewoms/models/pvs/pvsproperties.hh>
#include <ewoms/linear/parallelistlbackend.hh>

#include "./tnoBrineCH4FluidSystem_varsalinity.hpp"
//#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>
#include <opm/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh> /*@\label{tutorial1:include-grid-manager}@*/
#include <ewoms/io/cubegridmanager.hh> 

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

namespace Ewoms {
template <class TypeTag>
class WaterCH4Problem;
}

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(WaterCH4BaseProblem);

// Set the grid type
SET_TYPE_PROP(WaterCH4BaseProblem, Grid, Dune::YaspGrid<2>);
SET_TYPE_PROP(WaterCH4BaseProblem, GridManager, Ewoms::CubeGridManager<TypeTag>);

// Set the problem property
SET_TYPE_PROP(WaterCH4BaseProblem, Problem, Ewoms::WaterCH4Problem<TypeTag>);

// Set the material Law
SET_PROP(WaterCH4BaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffMaterialLaw;

public:
    // define the material law parameterized by absolute saturations
    // which uses the two-phase API
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Set the fluid system. in this case, we use the one which describes
// CH4 and water
SET_TYPE_PROP(WaterCH4BaseProblem, FluidSystem,
              Opm::FluidSystems::TNOBRINE_CH4<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(WaterCH4BaseProblem, EnableGravity, false);

SET_BOOL_PROP(WaterCH4BaseProblem,VtkWriteMassFractions,true);

// Use forward differences instead of central differences
//SET_INT_PROP(WaterCH4BaseProblem, NumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(WaterCH4BaseProblem, NewtonWriteConvergence, true);

// The default for the end time of the simulation (1 year)
SET_SCALAR_PROP(WaterCH4BaseProblem, EndTime, 43200000);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(WaterCH4BaseProblem, InitialTimeStepSize, 10);

SET_SCALAR_PROP(WaterCH4BaseProblem, MaxTimeStepSize, 864000);

// define the physical size of the problem's domain [m]
SET_SCALAR_PROP(WaterCH4BaseProblem, DomainSizeX, 300.0);
SET_SCALAR_PROP(WaterCH4BaseProblem, DomainSizeY, 60.0);
SET_SCALAR_PROP(WaterCH4BaseProblem, DomainSizeZ, 1.0);

// // define the number of cells used for discretizing the physical domain
SET_INT_PROP(WaterCH4BaseProblem, CellsX, 100);
SET_INT_PROP(WaterCH4BaseProblem, CellsY, 1);
SET_INT_PROP(WaterCH4BaseProblem, CellsZ, 1);

// Use the restarted GMRES linear solver with the ILU-2 preconditioner from dune-istl
SET_TAG_PROP(WaterCH4BaseProblem, LinearSolverSplice, ParallelIstlLinearSolver);
SET_TYPE_PROP(WaterCH4BaseProblem, LinearSolverWrapper,
              Ewoms::Linear::SolverWrapperRestartedGMRes<TypeTag>);
SET_TYPE_PROP(WaterCH4BaseProblem, PreconditionerWrapper,
              Ewoms::Linear::PreconditionerWrapperILUn<TypeTag>);
SET_INT_PROP(WaterCH4BaseProblem, PreconditionerOrder, 1);

} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(WaterCH4Problem, INHERITS_FROM(PvsModel, WaterCH4BaseProblem));

SET_BOOL_PROP(WaterCH4Problem, EnableEnergy, false);
}}

namespace Ewoms {

template <class TypeTag >
class WaterCH4Problem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,

        // component indices
        H2OIdx = FluidSystem::H2OIdx,
        CH4Idx = FluidSystem::CH4Idx,
        NACLIdx = FluidSystem::NACLIdx,

        // phase indices
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        
        //Variable indices
        pressure0Idx = Indices::pressure0Idx,
        //switch0Idx = Indices::switch0Idx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    static const bool enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy);

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    WaterCH4Problem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();
        eps_ = 3e-6;
        
        systemperature_ = 418.15;
        
        salinity_ = 0.1;
        max_salinity_ = 0.3;

        FluidSystem::init(/*Tmin=*/417, /*Tmax=*/420, /*nT=*/3,
                          /*pmin=*/1e5, /*pmax=*/12e6, /*np=*/10000, salinity_);

        // intrinsic permeabilities
        coarseK_ = this->toDimMatrix_(0.89e-14);

        // porosities
        coarsePorosity_ = 0.11;

        // residual saturations
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.164);
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.001);

        // parameters for the Brooks-Corey law
        coarseMaterialParams_.setEntryPressure(5000);
        coarseMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.finalize();
        initFluidStates_();
        this->simulator().startNextEpisode(24.0*60.0*60.0);
        ElementContext elemCtx(this->simulator());
        salt_storage_.resize(elemCtx.numDof(0),0);
        xyCoord_.resize(elemCtx.numDof(0));
        salt_file.open("./salt_storage.dat");
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "waterCH4_" << Model::name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";

        return oss.str();
    }
    
    void endEpisode()
    {
        // in the second episode, the actual work is done (the first is "settle down"
        // episode). we need to use a pretty short initial time step here as the change
        // in conditions is quite abrupt.
        if (this->simulator().time() < 10*86400)
            this->simulator().startNextEpisode(86400.0);
        else if (this->simulator().time() < 100*86400)
            this->simulator().startNextEpisode(10*86400.0);
        else
            this->simulator().startNextEpisode(100*86400.0);
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
//#ifndef NDEBUG
        // checkConservativeness() does not include the effect of constraints, so we
        // disable it for this problem...
        //this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl << std::flush;
        }
        GlobalPosition pos;
        Scalar x,y;
        for ( auto it = 0; it < salt_storage_.size(); ++it )
        {
          if (salt_storage_.at(it) != 0){
            pos = xyCoord_.at(it);
            x = pos[0];
            y = pos[1];
            salt_file << this->simulator().time() << ";" << it << "," << x << "," << y  << "," 
                      << salt_storage_.at(it) << "\n" << std::flush;
          }
        }
//#endif // NDEBUG
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     *
     * In this problem, the upper part of the domain is sightly less
     * permeable than the lower one.
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
            return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx,
                                               unsigned timeIdx) const
    {
            return coarseMaterialParams_;
    }
    // Get Temperature of the system. Required by multiphasebaseproblem.hh
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    {
        return systemperature_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * For this problem, we inject CH4 at the inlet on the center of
     * the lower domain boundary and use a no-flow condition on the
     * top boundary and a and a free-flow condition on the left and
     * right boundaries of the domain.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        assert(onLeftBoundary_(pos) ||
               onLowerBoundary_(pos) ||
               onRightBoundary_(pos) ||
               onUpperBoundary_(pos));
        
        if(onLeftBoundary_(pos)){          
            // impose an freeflow boundary condition
            values.setOutFlow(context, spaceIdx, timeIdx, outletFluidState_);
        }
        else if(onRightBoundary_(pos)){
            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, inletFluidState_);
        }
        else
            values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * For this problem, we set the medium to be fully saturated by
     * liquid water and assume hydrostatic pressure.
     */
    template <class Context>
    void initial(PrimaryVariables& values,
                 const Context& context,
                 unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(initFluidState_, matParams, false);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(  RateVector& rate,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx        ) const
    {
		    const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        const Scalar& elemVol = context.dofVolume(spaceIdx,timeIdx);
        const unsigned& globalIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        
        const auto& fs = context.intensiveQuantities(spaceIdx,timeIdx).fluidState();
        
        Scalar moleFracNaCl_Max_lPhase = FluidSystem::salinityTomoleFrac(max_salinity_);
        
        auto NaCll_moleFraction = (fs).moleFraction(liquidPhaseIdx,NACLIdx);
        const auto& elem_Porosity = context.intensiveQuantities(spaceIdx,timeIdx).porosity();
        const auto& elem_molarDensity = (fs).molarDensity(liquidPhaseIdx);
        const auto& elem_saturation = (fs).saturation(liquidPhaseIdx);
        Scalar precip_salt = elem_Porosity * elem_molarDensity * elem_saturation
                             * (NaCll_moleFraction - moleFracNaCl_Max_lPhase);

        Scalar removal = std::max(precip_salt,0.);
        
        RateVector molar_Rate(0.0);

        molar_Rate[conti0EqIdx + NACLIdx] = -removal;
        salt_storage_.at(globalIdx) += removal*elemVol;
        xyCoord_.at(globalIdx) = pos;
        rate.setMolarRate(molar_Rate);
    }

    /*!
    * \brief Evaluate the source term for all phases within a given
    *        sub-control-volume.
    *
    * For this method, the \a values parameter stores the rate mass
    * of a component is generated or annihilate per volume
    * unit. Positive values mean that mass is created, negative ones
    * mean that it vanishes.
    */
    /*
    void solDependentSource(PrimaryVariables &source,
                                const Element &element,
                                const FVElementGeometry &fvGeometry,
                                int scvIdx, const ElementVolumeVariables &elemVolVars) const
    {
                  source = 0;
                  const  VolumeVariables &volVars = elemVolVars[scvIdx];
                  Scalar moleFracNaCl_lPhase = volVars.fluidState().moleFraction(wPhaseIdx, NaClIdx);
                  Scalar moleFracNaCl_gPhase = volVars.fluidState().moleFraction(nPhaseIdx, NaClIdx);
                  Scalar massFracNaCl_Max_lPhase = this->spatialParams().SolubilityLimit();
                  Scalar moleFracNaCl_Max_lPhase = massTomoleFrac_(massFracNaCl_Max_lPhase);
                  Scalar moleFracNaCl_Max_gPhase = moleFracNaCl_Max_lPhase / volVars.fluidState().pressure(nPhaseIdx);

                  // liquid phase
                  Scalar precipSalt = volVars.porosity() * volVars.molarDensity(wPhaseIdx)
                                                             * volVars.saturation(wPhaseIdx)
                                                             * pow(abs(moleFracNaCl_lPhase - moleFracNaCl_Max_lPhase), 1.0);
                  if (moleFracNaCl_lPhase < moleFracNaCl_Max_lPhase)
                                precipSalt *= -1;

                  // gas phase
                  if (moleFracNaCl_gPhase > moleFracNaCl_Max_gPhase)
                                precipSalt += volVars.porosity() * volVars.molarDensity(nPhaseIdx)
                                                             * volVars.saturation(nPhaseIdx)
                                                             * pow(abs(moleFracNaCl_gPhase - moleFracNaCl_Max_gPhase), 1.0);

                  // make sure we don't disolve more salt than previously precipitated
                  if (precipSalt*this->timeManager().timeStepSize() + volVars.solidity(sPhaseIdx)* volVars.molarDensity(sPhaseIdx)< 0)
                                precipSalt = - volVars.solidity(sPhaseIdx)* volVars.molarDensity(sPhaseIdx)/this->timeManager().timeStepSize();

                  if (volVars.solidity(sPhaseIdx) >= volVars.InitialPorosity() - saltPorosity_ && precipSalt > 0)
                                //if (volVars.solidity(sPhaseIdx) >= 0.1*volVars.InitialPorosity()  && precipSalt > 0)
                  precipSalt = 0;

                  source[conti0EqIdx + NaClIdx] += -precipSalt;
                  source[precipNaClEqIdx] += precipSalt;

                  Valgrind::CheckDefined(source);
    }
    */


    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }
    
    void initFluidStates_()
    {
        updateFluidState_(inletFluidState_,true,false);
        updateFluidState_(outletFluidState_,false,false);
        updateFluidState_(initFluidState_,false,true);
    }

    template <class FluidState>
    void updateFluidState_(FluidState& fs,
                            bool isInlet,
                            bool isInit) const
    {
        int b;
        Scalar xNACL;
        xNACL = FluidSystem::salinityTomoleFrac(salinity_);
        Opm::MMPCAuxConstraint<Scalar> MMPCAux;
        if (isInlet&&!isInit){
            fs.setPressure(gasPhaseIdx,10.0e6);
            fs.setSaturation(liquidPhaseIdx, 0.2);
            fs.setTemperature(systemperature_);
            fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 0.8);
            fs.setMoleFraction(liquidPhaseIdx, NACLIdx, xNACL);
            fs.setMoleFraction(liquidPhaseIdx, CH4Idx, 1.0-0.8-xNACL);
            MMPCAux.set(liquidPhaseIdx,NACLIdx,xNACL);
            b = printf("Boundary is inlet\n");
        }
        else if(!isInit){
            fs.setPressure(gasPhaseIdx,3.0e6);
            fs.setSaturation(liquidPhaseIdx, 0.1);
            fs.setTemperature(systemperature_);
            fs.setMoleFraction(gasPhaseIdx, H2OIdx, 0.1);
            fs.setMoleFraction(gasPhaseIdx, NACLIdx, 0.);
            fs.setMoleFraction(gasPhaseIdx, CH4Idx, 0.9);
            MMPCAux.set(gasPhaseIdx,NACLIdx,0.);
            b = printf("Boundary is outlet\n");
        }
        if(isInit){
			// We can copy the inlet fluid state but this provides more flexibility, in case we want to change stuff for init runs
            fs.setPressure(gasPhaseIdx,10.0e6);
            fs.setSaturation(liquidPhaseIdx, 0.2);
            fs.setTemperature(systemperature_);
            fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 0.8);
            fs.setMoleFraction(liquidPhaseIdx, NACLIdx, xNACL);
            fs.setMoleFraction(liquidPhaseIdx, CH4Idx, 1.0-0.8-xNACL);
            MMPCAux.set(liquidPhaseIdx,NACLIdx,xNACL);
            b = printf("Initial internal fluid state\n");
        }
        // set the gas saturation and pressure
        fs.setSaturation(gasPhaseIdx, 1.0 - fs.saturation(liquidPhaseIdx) );
        PhaseVector pc;
        MaterialLaw::capillaryPressures(pc, coarseMaterialParams_, fs);
        fs.setPressure(liquidPhaseIdx, fs.pressure(gasPhaseIdx) - (pc[gasPhaseIdx] - pc[liquidPhaseIdx])); // Non intuitive implementation of the index, but because of RegularizedBrooksCorey.hpp, the pc[nonWetingindex] = pc(Sat_wetting) -KKL
        
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
        MMPC::solve(fs, paramCache, 0xffffff, &MMPCAux, 1, /*setViscosity=*/false, /*setEnthalpy=*/false);
        //typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        //CFRP::solve(fs, paramCache, liquidPhaseIdx, /*setViscosity=*/false,  /*setEnthalpy=*/false);
    }
    
    DimMatrix coarseK_;

    Scalar coarsePorosity_;

    MaterialLawParams coarseMaterialParams_;

    Scalar eps_;
    
    Scalar systemperature_;
    
    Scalar salinity_, max_salinity_;

    mutable std::vector<Scalar> salt_storage_;
    mutable std::vector<GlobalPosition> xyCoord_;
    std::ofstream salt_file;
    
    Opm::CompositionalFluidState<Scalar,FluidSystem> inletFluidState_;
    Opm::CompositionalFluidState<Scalar,FluidSystem> outletFluidState_;
    Opm::CompositionalFluidState<Scalar,FluidSystem> initFluidState_;
    
};
} // namespace Ewoms

#endif
