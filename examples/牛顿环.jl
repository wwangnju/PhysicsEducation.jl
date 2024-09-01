using PhysicsEducation
using GLMakie
using Makie.Colors
function newton_ring_figure()
    λ = 600.0e-9 #unit m
    R = 1.0
    t0 = 0.0
    p = [λ, R]
    xy = Vector{Float64}[collect(-1e-3:1e-6:1e-3), collect(-1e-3:1e-6:1e-3) ]
 
    doubleslit = WaveFunction(newton_ring_intensity, t0, xy, p)
    dso = WaveObservable([doubleslit])

    ps = [Dict(1=>4e-7:1e-8:8e-6, 2=>0.5:0.1:4.0)]
    pnames = [Dict(1=>"波长λ", 2=>"曲率半径R")]
    formats=[Dict(1=>"{:.2e}m", 2=>"{:.1f}m")]
    colors = [:grays]

    fig = Figure(;size = (1000, 1000), fontsize=15)

    intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        formats=formats,
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true),
        scale=(x->x),
        slider=:y,
        sliderkwargs=(;linestyle=:dash, color=:white)
    )
    # Label(fig[0, :], "Newton's Ring", fontsize = 30)
    Label(fig[0, :], "牛顿环", fontsize = 30)
    colsize!(fig.layout, 1, Aspect(1,1.0) )
    return fig, dso
end
fig, = newton_ring_figure()
fig