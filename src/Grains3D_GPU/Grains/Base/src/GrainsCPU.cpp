#include "GrainsCPU.hh"
#include "Grains.hh"
#include "GrainsParameters.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
GrainsCPU<T>::GrainsCPU()
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
GrainsCPU<T>::~GrainsCPU()
{
}

// -------------------------------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
template <typename T>
void GrainsCPU<T>::simulate()
{
    using G  = Grains<T>;
    using GP = GrainsParameters<T>;
    auto& SS = GP::m_simulationState;

    Gout(std::string(80, '='));
    Gout("Starting the simulation on CPU");
    Gout(std::string(80, '='));

    // first, inserting particles on host
    auto& timer = GP::m_simTimer;
    timer.start(SimStage::Total);

    timer.start(SimStage::Insertion);
    G::m_components->insertParticles(G::m_insertion);
    G::m_components->updateSubBodyPositions();
    timer.stop(SimStage::Insertion);

    printf("%12s  %12s  %7s  %8s  %8s  %8s  %5s  %7s  %7s\n",
           "SimTime",
           "EndTime",
           "Prog%",
           "Wall(s)",
           "ETA(s)",
           "stp/s",
           "NL#",
           "Pairs",
           "CPU_MB");

    // Pre-compute forces on the initial configuration so the first
    // moveParticles call (and the KDK first half-kick) has valid accelerations.
    if(GP::m_isLeapFrog)
    {
        timer.start(SimStage::DetectCollisions);
        G::m_components->detectCollisions();
        timer.stop(SimStage::DetectCollisions);

        timer.start(SimStage::ComputeContactForces);
        G::m_components->computeContactForces(G::m_contactForce);
        timer.stop(SimStage::ComputeContactForces);
    }

    // Write initial state (t = tStart) before advancing
    SS.time = GP::m_tStart;
    timer.start(SimStage::PostProcess);
    G::postProcess(G::m_components);
    timer.stop(SimStage::PostProcess);

    // time marching
    uint stepCount  = 0;
    GP::m_loopStart = std::chrono::steady_clock::now();
    for(SS.time = GP::m_tStart + GP::m_dt; SS.time <= GP::m_tEnd; SS.time += GP::m_dt)
    {
        stepCount++;
        // Output time
        if(GP::m_verbosityFrequency > 0 && (stepCount % GP::m_verbosityFrequency == 0))
            G::printStepInfo(stepCount, G::m_components->getNumberOfPairs());

        if(GP::m_isLeapFrog)
        {
            // KDK Step 1: half-kick + drift using f_n (from pre-loop or previous step).
            timer.start(SimStage::MoveComponents);
            G::m_components->moveComponents(G::m_timeIntegrator);
            G::m_components->updateSubBodyPositions();
            timer.stop(SimStage::MoveComponents);
            timer.start(SimStage::DetectCollisions);
            G::m_components->detectCollisions();
            timer.stop(SimStage::DetectCollisions);
            timer.start(SimStage::ComputeContactForces);
            G::m_components->computeContactForces(G::m_contactForce);
            timer.stop(SimStage::ComputeContactForces);
            // KDK Step 3: second half-kick using f_{n+1}.
            timer.start(SimStage::AdvanceVelocity);
            G::m_components->advanceVelocity(G::m_timeIntegrator);
            timer.stop(SimStage::AdvanceVelocity);
        }
        else
        {
            // Single-pass scheme: compute forces at x_n, then advance.
            timer.start(SimStage::DetectCollisions);
            G::m_components->detectCollisions();
            timer.stop(SimStage::DetectCollisions);
            timer.start(SimStage::ComputeContactForces);
            G::m_components->computeContactForces(G::m_contactForce);
            timer.stop(SimStage::ComputeContactForces);
            timer.start(SimStage::MoveComponents);
            G::m_components->moveComponents(G::m_timeIntegrator);
            timer.stop(SimStage::MoveComponents);
        }

        // Post-Processing
        timer.start(SimStage::PostProcess);
        G::postProcess(G::m_components);
        timer.stop(SimStage::PostProcess);
    }

    timer.stop(SimStage::Total);
    if(timer.isEnabled())
    {
        if(GP::m_cdmTimer.isEnabled())
            GP::m_cdmTimer.printSummary();
        if(GP::m_fmTimer.isEnabled())
            GP::m_fmTimer.printSummary();
        timer.printSummary();
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class GrainsCPU<float>;
template class GrainsCPU<double>;