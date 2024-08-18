using PhysicsEducation
using GLMakie
using Makie.Colors
function circular_fraunhofer_figure()
    λ = 600.0e-9 #unit m
    W, L = 7.5e-4, 1.0
    t0 = 0.0
    p = [λ, W, L]
    xy = Vector{Float64}[collect(-0.01:1e-5:0.010), collect(-0.01:1e-4:0.01) ]
 
    slit = WaveFunction(circular_fraunhofer, t0, xy, p)
    dso = WaveObservable([slit])
    ps = [Dict(1=>4e-7:1e-8:8e-6, 2=>1e-5:1e-5:1e-3, 3=>0.5:0.01:2.0)]
    pnames = [Dict(1=>"波长λ", 2=>"孔宽W", 3=>"焦距L")]
    colors = [:grays]

    fig = Figure(;resolution = (1000, 800))

    ax, ax1 = intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true, colorrange=(0,0.005)),
        scale=x->(x),
        slider=:y,
        sliderkwargs=(;linestyle=:dash, color=:white)
    )
    # Label(fig[0, :], "Circular Aperture", fontsize = 30)
    Label(fig[0, :], "圆孔夫琅和费衍射", fontsize = 30)
    vlines!(ax1, [1.22*λ/W], color=:red)
    return fig, dso
end
fig, =circular_fraunhofer_figure()
fig