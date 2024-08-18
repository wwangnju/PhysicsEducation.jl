using PhysicsEducation
using GLMakie
using Makie.Colors
function slit_grating_figure()
    λ = 600.0e-9 #unit m
    W, L = 1e-6, 1.0
    d = 2e-6
    N = 10
    t0 = 0.0
    p = [λ, d, W, L, N]
    xy = Vector{Float64}[collect(-1.0:1e-3:1.0), collect(-1:0.2:1) ]
 
    slit = WaveFunction(slit_grating, t0, xy, p)
    dso = WaveObservable([slit])
    ps = [Dict(1=>4e-7:1e-8:8e-6, 2=>1e-7:1e-7:1e-5, 3=>0.0:1e-7:5e-6, 5=>1:1:100, 4=>0.01:0.01:3.0)]
    pnames = [Dict(1=>"波长λ", 2=>"光栅常数d", 3=>"缝宽W", 5=>"缝数N", 4=>"缝屏距L")]
    colors = [:grays]

    fig = Figure(;size = (1000, 1000), fontsize=15)

    ax, ax1 = intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true),
        scale=x->(x),
        slider=:y,
        sliderkwargs=(;linestyle=:dash, color=:white)
    )
    # Label(fig[0, :], "Slit Grating", fontsize = 30)
    Label(fig[0, :], "光栅夫琅和费衍射", fontsize = 30)
    return fig, dso
end
fig, = slit_grating_figure()
fig