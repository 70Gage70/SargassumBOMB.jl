using MAT
using Interpolations

include("interpolant-helpers.jl")
include(joinpath(@__DIR__, "..", "src", "coordinates.jl"))

############################################################################################

"""
    struct GriddedField{T, A, R}

A container for gridded, possibly time-dependent, field data. 

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

    vars = read(data, var_names...) .|> vec |> collect # want vectors (not N x 1 matrices) so map each var to `vec`
    fields = read(data, field_names...) |> collect
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
            vars[i] = reverse(vars[i])
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

"""
    sph2xy(gridded_field; lon_name = :lon, lat_name = :lat, x_name = :x, y_name = :y, xy_units = "km")

Transform the `lon`, `lat` variables in `gridded_field` to `x`, `y` variables with reference `gridded_field.ref` and
return a new `GriddedField`.

The units of `x` and `y` are the same as `ref.R`.

The `lon` and `lat` variables in `gridded_field` should have names `lat_name` (default `:lat`) and `lon_name` (default `:lon`).

The units of `x` and `y` are updated to `xy_units` (default "km").

The returned `GriddedField` has its variable names updated to `x_name` (default `:x`) and `y_name` (default `:y`).
"""
function sph2xy(
    gridded_field::GriddedField;
    lon_name::Symbol = :lon, 
    lat_name::Symbol = :lat, 
    x_name::Symbol = :x, 
    y_name::Symbol = :y,
    xy_units::String = "km")

    lon = gridded_field.vars[lon_name]
    lat = gridded_field.vars[lat_name]
    x, y = sph2xy(lon, lat, gridded_field.ref)

    new_var_names = Vector{Symbol}()
    new_var_units = typeof(gridded_field.vars_units)()
    new_vars = typeof(gridded_field.vars)()

    for var_name in gridded_field.var_names
        if var_name == lon_name
            push!(new_var_names, x_name)
            new_var_units[x_name] = xy_units
            new_vars[x_name] = x
        elseif var_name == lat_name
            push!(new_var_names, y_name)
            new_var_units[y_name] = xy_units
            new_vars[y_name] = y
        else
            push!(new_var_names, var_name)
            new_vars[var_name] = gridded_field.vars[var_name]
            new_var_units[var_name] = gridded_field.vars_units[var_name]
        end
    end

    return GriddedField(new_var_names, new_vars, gridded_field.fields, new_var_units, gridded_field.fields_units, gridded_field.time_start, gridded_field.ref)
end

"""
    struct InterpolatedField{T, A, R}

Indentical to [`GriddedField`](@ref) except `fields` is a map from field names to interpolating functions 
rather than arrays.

### Constructors 

Refer to `interpolate(gridded_field::GriddedField)`.
"""
struct InterpolatedField{R<:AbstractRange, I<:AbstractInterpolation, Q<:Union{EquirectangularReference{<:Real}, Nothing}}
    var_names::Vector{Symbol}
    vars::Dict{Symbol, R}
    fields::Dict{Symbol, I}
    vars_units::Union{Dict{Symbol, String}, Nothing}
    fields_units::Union{Dict{Symbol, String}, Nothing}
    time_start::Union{DateTime, Nothing}
    ref::Q
end

"""
    interpolate(gridded_field; interpolant_type, extrapolate_value)

Create an `InterpolatedField` from `gridded_field::GriddedField` using interpolation `interpolant_type` and with 
a constant extrapolat `extrapolate_value`.

### Optional Arguments

- `interpolant_type`: Two convenience flags are provided, `"cubic"` and `"nearest` which refer to cubic BSpline and nearest-neighbor interpolation, respectively. Alternatively, any `Interpolations.InterpolationType` can be provided. Default `"cubic"`.
- `extrapolate_value`: A constant extrapolation is performed with this value. Default `"0.0"`.
"""
function interpolate(
    gridded_field::GriddedField; 
    interpolant_type::Union{String, Interpolations.InterpolationType} = "cubic",
    extrapolate_value::Real = 0.0)

    if interpolant_type isa String
        @assert interpolant_type in ["cubic", "nearest"] "kwarg `interpolant_type` should be either $("cubic"), $("nearest") or an `Interpolations.InterpolationType`."

        if interpolant_type == "cubic"
            spline = BSpline(Cubic(Interpolations.Line(OnGrid())))
        elseif interpolant_type == "nearest"
            spline = BSpline(Constant)
        end
    else
        spline = interpolant_type
    end

    vars = [gridded_field.vars[name] for name in gridded_field.var_names]
    field_names = collect(keys(gridded_field.fields))
    interps = [
        extrapolate(
            scale(
                Interpolations.interpolate(gridded_field.fields[name], spline), 
            vars...), 
        extrapolate_value) 
        for name in field_names]
    fields_dict = Dict(field_names .=> interps)

    return InterpolatedField(
        gridded_field.var_names, 
        gridded_field.vars, 
        fields_dict, 
        gridded_field.vars_units, 
        gridded_field.fields_units, 
        gridded_field.time_start, 
        gridded_field.ref)
end

function Base.show(io::IO, x::Union{GriddedField, InterpolatedField})
    if x isa GriddedField
        println(io, "GriddedField")
    else
        println(io, "InterpolatedField")
    end

    if x.vars_units !== nothing
        println(io, " Variables: ", [(name, x.vars_units[name]) for name in x.var_names])
    else
        println(io, " Variables: ", x.var_names)
    end

    println(io, " Dimensions: ", [length(x.vars[name]) for name in x.var_names])

    if x.fields_units !== nothing
        println(io, " Fields: ", [(name, x.fields_units[name]) for name in keys(x.fields)])
    else
        println(io, " Fields: ", collect(keys(x.fields)))
    end

    if x.time_start !== nothing
        println(io, " Time Start: ", x.time_start)
    else
        println(io, " Time Start: ", "Time-independent")
    end

    if x.ref !== nothing
        println(io, " Ref: ", "lon0 = ", x.ref.lon0, ", lat0 = ", x.ref.lat0)
    end  

end