
using PhysicsEducation
using GLMakie
using Makie.Colors
function michelsoninterferometer_figure()
    A, λ = 1.0, 600.0e-9 #unit m
    d = 0.08e-3 #m
    t0 = 0.0
    p = [A, λ, d]
    xy = Vector{Float64}[collect(-0.2:1e-3:0.2), collect(-0.2:1e-3:0.2) ]
 
    doubleslit = WaveFunction(michelsoninterferometer_intensity, t0, xy, p)
    dso = WaveObservable([doubleslit]; current_step=0, Δt=0.01)

    ps = [Dict(2=>4e-7:1e-8:8e-6, 3=>1e-5:1e-6:0.0001)]
    pnames = [Dict(2=>"波长λ", 3=>"移动距离d")]
    formats=[Dict(2=>"{:.2e}m", 3=>"{:.2e}m")]
    colors = [:grays]

    fig = Figure(;resolution = (1000, 800))

    intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        formats=formats,
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true),
        slider=:y,
        sliderkwargs=(;linestyle=:dash, color=:white)
    )
    colsize!(fig.layout, 1, Aspect(1,1.0))
    Label(fig[0, :], "迈克尔逊干涉", fontsize = 30)
    return fig, dso
end
fig, = michelsoninterferometer_figure()
fig