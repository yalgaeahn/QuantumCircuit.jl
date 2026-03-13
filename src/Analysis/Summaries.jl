struct SweepSummary{T}
    values::Vector{T}
    transition_01::Vector{Float64}
    anharmonicity::Vector{Float64}
end

function sweep_summary(result::SweepResult)
    transition_01 = transition_curve(result; transition = (1, 2)).data
    anharmonicity_values = anharmonicity_curve(result).data
    return SweepSummary(copy(result.values), transition_01, anharmonicity_values)
end
