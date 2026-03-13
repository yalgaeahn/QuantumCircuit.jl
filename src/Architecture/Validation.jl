function validate_identifier(id::Symbol, kind::AbstractString)
    id == Symbol("") && throw(ArgumentError("$kind name must not be empty."))
    return id
end

function validate_positive(value::Real, label::AbstractString)
    isfinite(value) || throw(ArgumentError("$label must be finite."))
    value > 0 || throw(ArgumentError("$label must be positive."))
    return Float64(value)
end

function validate_finite(value::Real, label::AbstractString)
    isfinite(value) || throw(ArgumentError("$label must be finite."))
    return Float64(value)
end

function validate_abs_leq_one(value::Real, label::AbstractString)
    isfinite(value) || throw(ArgumentError("$label must be finite."))
    abs(value) <= 1 || throw(ArgumentError("$label must satisfy |value| <= 1."))
    return Float64(value)
end

function validate_positive_integer(value::Integer, label::AbstractString)
    value > 0 || throw(ArgumentError("$label must be a positive integer."))
    return Int(value)
end
