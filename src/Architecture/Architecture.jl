module Architecture

export AbstractCoupling,
    AbstractSubsystem,
    CapacitiveCoupling,
    CircuitCapacitiveCoupling,
    CompositeSystem,
    Resonator,
    TunableCoupler,
    TunableTransmon,
    Transmon,
    coupling_endpoints,
    couplings,
    name,
    subsystem_names,
    subsystems,
    with_coupling_parameter,
    with_subsystem_parameter

include("Validation.jl")
include("Subsystems.jl")
include("Couplings.jl")
include("CompositeSystem.jl")

end
