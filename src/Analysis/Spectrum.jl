function transition_frequencies(result::SpectrumResult)
    return diff(result.energies)
end

function anharmonicity(result::SpectrumResult)
    length(result.energies) >= 3 || throw(ArgumentError("anharmonicity requires at least three energy levels."))
    ω01 = _transition_frequency(result, 1, 2)
    ω12 = _transition_frequency(result, 2, 3)
    return ω12 - ω01
end

function _transition_frequency(result::SpectrumResult, lower::Integer, upper::Integer)
    _validate_level_pair(lower, upper, length(result.energies))
    return result.energies[upper] - result.energies[lower]
end

function _validate_level_pair(lower::Integer, upper::Integer, available_levels::Integer)
    lower > 0 || throw(ArgumentError("Level indices must be positive integers."))
    upper > lower || throw(ArgumentError("Upper level index must be greater than lower level index."))
    upper <= available_levels || throw(
        ArgumentError(
            "Requested level pair ($(lower), $(upper)) exceeds the available spectrum size $available_levels.",
        ),
    )
    return nothing
end
