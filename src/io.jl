"""
    rtr2mat(rtr, outfile; force)

Write the [`RaftTrajectory`](@ref) in `rtr` to `outfile` which must be a `.mat` file.

This writes the raw, unbinned trajectory data.

### Optional Arguments

- `force`: If `true`, delete `outfile` if it already exists. Default `false`.
"""
function rtr2mat(rtr::RaftTrajectory, outfile::String; force::Bool = false)

    @assert endswith(outfile, ".mat") "Must output a .mat file."
    !force && @assert !isfile(outfile) "This file already exists. Pass `force = true` to overwrite it."

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