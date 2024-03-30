"""
    rtr2mat(rtr, outfile; force)

Write the [`RaftTrajectory`](@ref) in `rtr` to `outfile` which must be a `.mat` file.

This writes the raw, unbinned trajectory data.

### Optional Arguments

- `force`: If `true`, delete `outfile` if it already exists. Default `false`.
"""
function rtr2mat(rtr::RaftTrajectory, outfile::String; force::Bool = false)

    @argcheck endswith(outfile, ".mat") "Must output a .mat file."
    !force && @argcheck !isfile(outfile) "This file already exists. Pass `force = true` to overwrite it."

    matdict_t = rtr.t
    matdict_n_clumps = rtr.n_clumps
    matdict_com = [rtr.com.xy ;; rtr.com.t]
    matdict_info = "The trajectory data is a matrix with N rows and 4 columns. 

    Each row is a data point of the form (i, x, y, t) where i is the index of the clump.

    Time t is measured in $(UNITS["time"]) since $(T_REF.x)." 

    matdict_traj = Vector{Vector{Float64}}()
    for clump_id in sort(collect(keys(rtr.trajectories)))
        traj = rtr.trajectories[clump_id]
        for i = 1:length(traj.t)
            push!(matdict_traj, [clump_id, traj.xy[i,1], traj.xy[i,2], traj.t[i]])
        end
    end

    matdict_traj = stack(matdict_traj, dims = 1)

    mat_dict = Dict(
        "trajectories" => matdict_traj,
        "times" => matdict_t,
        "n_clumps" => matdict_n_clumps,
        "com" => matdict_com,
        "info" => matdict_info)

    isfile(outfile) && rm(outfile)
    matwrite(outfile, mat_dict)

    return nothing
end

"""
    rtr2nc(rtr, outfile, lon_bins, lat_bins; force)

Write the [`RaftTrajectory`](@ref) in `rtr` to `outfile` which must be a `.nc` file.

The data are binned by passing `lon_bins` and `lat_bins` to [`bins`](@ref).

It is required that `rtr.t` is linearly spaced.

### Optional Arguments

- `force`: If `true`, delete `outfile` if it already exists. Default `false`.
"""
function rtr2nc(rtr::RaftTrajectory, outfile::String, lon_bins::StepRangeLen, lat_bins::StepRangeLen; force::Bool = false)
    try 
        vec2range(rtr.t)
    catch
        error("`rtr.t` must be linearly spaced.")
    end

    @argcheck endswith(outfile, ".nc") "Must output a .nc file."
    !force && @argcheck !isfile(outfile) "This file already exists. Pass `force = true` to overwrite it."

    times = rtr.t
    data = zeros(UInt16, length(lon_bins)-1, length(lat_bins)-1, length(times))

    for i = 1:length(times)
        data[:,:,i] .= Integer.(bins(time_slice(rtr, (times[i], times[i])), lon_bins, lat_bins))
    end
    
    lon_bins_c = [(lon_bins[i + 1] - lon_bins[i])/2 for i = 1:length(lon_bins) - 1]
    lat_bins_c = [(lat_bins[i + 1] - lat_bins[i])/2 for i = 1:length(lat_bins) - 1]
    
    # attributes
    lonatts = Dict(
        "longname" => "Longitude",
        "units"    => "degrees east",
        "min"      => minimum(lon_bins_c),
        "max"      => maximum(lon_bins_c))
    latatts = Dict(
        "longname" => "Latitude",
        "units"    => "degrees north",
        "min"      => minimum(lat_bins_c),
        "max"      => maximum(lat_bins_c))
    timeatts = Dict(
        "longname" => "Time",
        "units"    => "$(UNITS["time"]) since $(T_REF.x)",
        "example"  => "time = 6677 is $(time2datetime(6677))",
        "min"      => minimum(times),
        "max"      => maximum(times),
        "min_datetime"      => string(minimum(time2datetime.(times))),
        "max_datetime"      => string(maximum(time2datetime.(times))))
    
    data_atts = Dict(
        "longname" => "Trajectory counts",
        "units"    => "number")
    
    # writing to file
    isfile(outfile) && rm(outfile)
    
    nccreate(outfile, 
        "data",
        "lon", lon_bins_c, lonatts,
        "lat", lat_bins_c, latatts, 
        "time", times, timeatts, 
        atts = data_atts,
        gatts = Dict(
            "info" =>       "This file was generated by the Julia package SargassumBOMB.jl. It contains trajectory counts of Sargassum clumps.",
            "github" =>     "https://github.com/70Gage70/SargassumBOMB.jl"),
        t = UInt16)
    
    ncwrite(data, outfile, "data")

    return nothing
end