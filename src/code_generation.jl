function generate_clumps(i)
    s = 
    """
    eqsC[$i] = [
        ddt(x[$i]) ~ u_x(x[$i], y[$i], t, α[$i]) + τ[$i] * (
            R[$i]*Dv_xDt(x[$i], y[$i], t) - R[$i]*(f[$i] + ω(x[$i], y[$i], t)/3)*v_y(x[$i], y[$i], t) - Du_xDt(x[$i], y[$i], t, α[$i]) + (f[$i] + R[$i]*ω(x[$i], y[$i], t)/3)*u_y(x[$i], y[$i], t, α[$i]) + Fx[$i]
        ),
        ddt(y[$i]) ~ u_y(x[$i], y[$i], t, α[$i]) + τ[$i] * (
            R[$i]*Dv_yDt(x[$i], y[$i], t) + R[$i]*(f[$i] + ω(x[$i], y[$i], t)/3)*v_x(x[$i], y[$i], t) - Du_yDt(x[$i], y[$i], t, α[$i]) - (f[$i] + R[$i]*ω(x[$i], y[$i], t)/3)*u_x(x[$i], y[$i], t, α[$i]) + Fy[$i]
        )
        ]
    """

    return Meta.parse(s)
end

function generate_forces(i)
    sx = "Fx[$i] ~ "
    for j = 1:25
        if j != i
            if j == 1 || (j == 2 && i == 1)
                sx = *(sx, "spring_force_x(x[$i], x[$j], y[$i], y[$j], springs[$i, $j])")
            else
                sx = *(sx, " + spring_force_x(x[$i], x[$j], y[$i], y[$j], springs[$i, $j])")
            end
        end
    end

    sy = "Fy[$i] ~ "
    for j = 1:25
        if j != i
            if j == 1 || (j == 2 && i == 1)
                sy = *(sy, "spring_force_y(x[$i], x[$j], y[$i], y[$j], springs[$i, $j])")
            else
                sy = *(sy, " + spring_force_y(x[$i], x[$j], y[$i], y[$j], springs[$i, $j])")
            end
        end
    end

    s = "eqsF[$i] = [$sx ,\n $sy]"

    return Meta.parse(s)
end

function generate_eqs()
    s = "[eqsC[1][1],eqsF[1][1],eqsC[1][2],eqsF[1][2]"
    for i = 2:25
        for j = 1:2
            s = *(s, ",eqsC[$i][$j]")
            s = *(s, ",eqsF[$i][$j]")
        end
    end

    s = *(s, "]")

    return Meta.parse(s)
end

