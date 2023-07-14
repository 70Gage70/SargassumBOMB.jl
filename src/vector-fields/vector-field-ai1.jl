#######################

function MaterialDerivativeX(water_vf::VectorField2DInterpolantEQR; interpolant_type = BSpline(Cubic(Line(OnGrid()))))
    x, y, time = (water_vf.x, water_vf.y, water_vf.time)
    data = [MaterialDerivativeX(water_vf, x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

function MaterialDerivativeY(water_vf::VectorField2DInterpolantEQR; interpolant_type = BSpline(Cubic(Line(OnGrid()))))
    x, y, time = (water_vf.x, water_vf.y, water_vf.time)
    data = [MaterialDerivativeY(water_vf, x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

function Vorticity(water_vf::VectorField2DInterpolantEQR; interpolant_type = BSpline(Cubic(Line(OnGrid()))))
    x, y, time = (water_vf.x, water_vf.y, water_vf.time)
    data = [Vorticity(water_vf, x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

function WindWaterAlphaX(
    wind_vf::VectorField2DInterpolantEQR, 
    water_vf::VectorField2DInterpolantEQR,
    clump_parameters::ClumpParameters; 
    interpolant_type = BSpline(Cubic(Line(OnGrid()))))

    α = clump_parameters.α

    x, y, time = (wind_vf.x, wind_vf.y, wind_vf.time)
    data = [(1 - α) * water_vf.u(x, y, t) + α * wind_vf.u(x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

function WindWaterAlphaY(
    wind_vf::VectorField2DInterpolantEQR, 
    water_vf::VectorField2DInterpolantEQR,
    clump_parameters::ClumpParameters; 
    interpolant_type = BSpline(Cubic(Line(OnGrid()))))

    α = clump_parameters.α

    x, y, time = (wind_vf.x, wind_vf.y, wind_vf.time)
    data = [(1 - α) * water_vf.v(x, y, t) + α * wind_vf.v(x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

function DDtWindWaterAlphaX(
    wind_vf::VectorField2DInterpolantEQR, 
    water_vf::VectorField2DInterpolantEQR,
    clump_parameters::ClumpParameters; 
    interpolant_type = BSpline(Cubic(Line(OnGrid()))))

    α = clump_parameters.α

    x, y, time = (wind_vf.x, wind_vf.y, wind_vf.time)
    data = [(1 - α) * MaterialDerivativeX(water_vf, x, y, t) + α * MaterialDerivativeX(wind_vf, x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

function DDtWindWaterAlphaY(
    wind_vf::VectorField2DInterpolantEQR, 
    water_vf::VectorField2DInterpolantEQR,
    clump_parameters::ClumpParameters; 
    interpolant_type = BSpline(Cubic(Line(OnGrid()))))

    α = clump_parameters.α

    x, y, time = (wind_vf.x, wind_vf.y, wind_vf.time)
    data = [(1 - α) * MaterialDerivativeY(water_vf, x, y, t) + α * MaterialDerivativeY(wind_vf, x, y, t) for x in x, y in y, t in time]
    return interpolate_field(x, y, time, data, interpolant_type = interpolant_type)
end

##############

function construct_itp_EQR(
    water_infile::String = water_file_default, 
    wind_infile::String = wind_file_default, 
    ref::EquirectangularReference = ref_default, 
    clump_parameters::ClumpParameters = ClumpParameters(ref); 
    outfile::String = "itp.jld2")

    @info "Constructing wind interpolant."
    wind_itp = VectorField2DGridSPH(wind_infile, lon_alias = "Lon", lat_alias = "Lat", lon_lat_time_order = [2, 1, 3])
    wind_itp = VectorField2DInterpolantEQR(wind_itp, ref)
    
    @info "Constructing water interpolant."
    water_itp = VectorField2DGridSPH(water_infile, lon_lat_time_order = [2, 1, 3])
    water_itp = VectorField2DInterpolantEQR(water_itp, ref)
    
    @info "Constructing DDt."
    MDX = MaterialDerivativeX(water_itp)
    MDY = MaterialDerivativeY(water_itp)
    V = Vorticity(water_itp)

    @info "Constructing U."
    UX = WindWaterAlphaX(wind_itp, water_itp, clump_parameters)    
    UY = WindWaterAlphaY(wind_itp, water_itp, clump_parameters) 
    DUX = DDtWindWaterAlphaX(wind_itp, water_itp, clump_parameters)
    DUY = DDtWindWaterAlphaY(wind_itp, water_itp, clump_parameters)

    jldsave(outfile;
        v_x = water_itp.u,
        v_y =  water_itp.v,
        Dv_xDt = MDX,
        Dv_yDt = MDY,
        u_x = UX,
        u_y = UY,
        Du_xDt = DUX, 
        Du_yDt = DUY,
        ω = V, 
        ref = ref
    )

    @info "Interpolants written to $(outfile)."

    return nothing

end