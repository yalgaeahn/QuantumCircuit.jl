struct SweepSeries{T}
    values::Vector{T}
    data::Vector{Float64}
    label::Symbol
end

function transition_curve(result::SweepResult; transition::Tuple{<:Integer,<:Integer} = (1, 2))
    lower, upper = transition
    data = [_transition_frequency(spectrum_result, lower, upper) for spectrum_result in result.spectra]
    label = Symbol("transition_$(lower - 1)$(upper - 1)")
    return SweepSeries(copy(result.values), data, label)
end

function transition_curve(result::EigensystemSweepResult; transition::Tuple{<:Integer,<:Integer} = (1, 2))
    lower, upper = transition
    data = [_transition_frequency(spectrum_result, lower, upper) for spectrum_result in result.spectra]
    label = Symbol("transition_$(lower - 1)$(upper - 1)")
    return SweepSeries(copy(result.values), data, label)
end

function anharmonicity_curve(result::SweepResult)
    data = [anharmonicity(spectrum_result) for spectrum_result in result.spectra]
    return SweepSeries(copy(result.values), data, :anharmonicity)
end

function anharmonicity_curve(result::EigensystemSweepResult)
    data = [anharmonicity(spectrum_result) for spectrum_result in result.spectra]
    return SweepSeries(copy(result.values), data, :anharmonicity)
end

function minimum_gap(result::SweepResult; level_pair::Tuple{<:Integer,<:Integer} = (1, 2))
    curve = transition_curve(result; transition = level_pair)
    index = argmin(curve.data)
    return (; gap = curve.data[index], sweep_value = curve.values[index], index, level_pair = level_pair)
end

function minimum_gap(result::EigensystemSweepResult; level_pair::Tuple{<:Integer,<:Integer} = (1, 2))
    curve = transition_curve(result; transition = level_pair)
    index = argmin(curve.data)
    return (; gap = curve.data[index], sweep_value = curve.values[index], index, level_pair = level_pair)
end
