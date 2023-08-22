using MAT
using Interpolations

include("interpolant-helpers.jl")
include(joinpath(@__DIR__, "..", "src", "coordinates.jl"))

############################################################################################

"""
    struct GriddedField{T, A, R}

A container for gridded, possible time-dependent, field data. 

### Fields
- `var_names`: A `Vector` of `Symbol`s such that `var_names[i]` is the `i`th dimension of `fields`.
- `vars`: A `Dict` mapping variable names to ranges they take.
- `fields`: A `Dict` mapping field names to their arrays.
- `vars_units`: A `Dict` mapping variable names to their units. (Optional)
- `fields_units`: A `Dict` mapping field names to their units.
- `time_start`: For time-dependent fields, this is the `DateTime` of the first entry of the time variable.
- `ref`: An [`EquirectangularReference`](@ref) providing the translation between spherical and equirectangular coordinates of the `vars`.
"""
struct GriddedField{R<:AbstractRange, A<:AbstractArray, Q<:Union{EquirectangularReference{<:Real}, Nothing}}
    var_names::Vector{Symbol}
    vars::Dict{Symbol, R}
    fields::Dict{Symbol, A}
    vars_units::Union{Dict{Symbol, String}, Nothing}
    fields_units::Union{Dict{Symbol, String}, Nothing}
    time_start::Union{DateTime, Nothing}
    ref::Q
end

"""
    GriddedField(infile, var_names, field_names; time_index, time2datetime, NaN_replacement, var_units, field_units, ref)  

Construct a [`GriddedField`](@ref) from the data in `infile`. The data in the file must be gridded beforehand. In addition, the `GriddedField` is 
constructed with the same variable and field names as the names in the input file.

For time-dependent fields, the time grid is rescaled to be the number of steps since the initial time. Refer to the optional arguments.

### Arguments

- `infile`: A `String` with the name of the input file. 
- `var_names`: A `Vector` of `String`s such that `var_names[i]` is the name of the `i`th dimension of the fields and `infile` contains a key with the same name.
- `field_names`: A `Vector` of `String`s such that `infile` contains keys with the same names.

### Optional Arguments 

All of the optional arguments default to `nothing`.

- `time_index`: An `Integer` `i` such that `var_names[i]` identifies the time variable.
- `time2datetime`: A `Function` which maps `var_names[time_index]` to `DateTime`s. Used to compute the initial time.
- `NaN_replacement`: Any `NaN`s encountered in the fields will be replaced by this value.
- `var_units`: A `Vector` of `String`s such that `var_units`[i] gives the units of `var_names`[i].
- `field_units`: A `Vector` of `String`s such that `field_names`[i] gives the units of `field_names`[i].
- `ref`: An [`EquirectangularReference`](@ref) for translation between spherical and equirectangular coordinates.   
"""
function GriddedField(
    infile::String,
    var_names::Vector{String},
    field_names::Vector{String};
    time_index::Union{Integer, Nothing} = nothing,
    time2datetime::Union{Function, Nothing} = nothing, 
    NaN_replacement::Any = nothing,
    var_units::Union{Vector{String}, Nothing} = nothing,
    field_units::Union{Vector{String}, Nothing} = nothing,
    ref::Union{EquirectangularReference{<:Real}, Nothing} = nothing
    )

    # ensure that the provided file is a mat file
    extension = infile[findlast(==('.'), infile)+1:end]
    @assert extension in ["mat"] "Require a .mat file."

    data = matopen(infile)

    # ensure that all the relevant variables are contained in the given file
    data_keys = collect(keys(data))

    for var_name in var_names
        @assert var_name in data_keys "`var_name` $(var_name) not in $(infile)"
    end

    for field_name in field_names
        @assert field_name in data_keys "`field_name` $(field_name) not in $(infile)"
    end   

    vars = read(data, var_names...) .|> vec # want vectors (not N x 1 matrices) so map each var to `vec`
    fields = read(data, field_names...)
    close(data)

    # ensure dimensions are consistent
    for i = 1:length(field_names)
        @assert length(var_names) == ndims(fields[i]) "Field $(field_names[i]) has size $(ndims(fields[i])) but have $(length(var_names)) variables."
    end    

    # clean variable grids to ranges
    vars = vars .|> reduce_vector_to_range

    # ensure that the ranges are increasing, otherwise reverse
    for i = 1:length(vars)
        if step(vars[i]) < 0
            reverse!.(fields, dims = i)
            reverse!(vars[i])
        end
    end

    # NaN replacement
    if NaN_replacement !== nothing
        for field in fields
            field[isnan.(field)] .= NaN_replacement
        end
    end

    # rescale time
    if time_index !== nothing
        time_start = first(vars[time_index])
        vars[time_index] = vars[time_index] .- time_start

        if time2datetime !== nothing
            time_start = time2datetime(time_start)
        else
            time_start = nothing
        end
    else 
        time_start = nothing
    end

    # ensure data are of matching types
    vars = promote(vars...) |> collect
    fields = promote(fields...) |> collect

    # symbolizing
    var_names_symb = Symbol.(var_names)
    vars_dict = Dict(Symbol.(var_names) .=> vars)
    fields_dict = Dict(Symbol.(field_names) .=> fields)
    vars_units_dict = Dict(Symbol.(var_names) .=> var_units)
    fields_units_dict = Dict(Symbol.(field_names) .=> field_units)
 
    return GriddedField(var_names_symb, vars_dict, fields_dict, vars_units_dict, fields_units_dict, time_start, ref)
end

