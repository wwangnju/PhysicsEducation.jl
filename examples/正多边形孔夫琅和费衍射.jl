using PhysicsEducation
using GLMakie
using Makie.Colors
function regular_polygon_fraunhofer_figure()
    λ = 600.0e-9 #unit m
    R, L = 4e-3, 10.0
    n = 6
    t0 = 0.0
    p = [λ, R, L, n ]
    xy = Vector{Float64}[collect(-0.01:1e-5:0.010), collect(-0.01:1e-4:0.01) ]
 
    slit = WaveFunction(regular_polygon, t0, xy, p)
    dso = WaveObservable([slit])
    ps = [Dict(1=>4e-7:1e-8:8e-6, 2=>1e-3:1e-4:1e-2, 3=>1:0.1:20.0, 4=>3:1:100)]
    pnames = [Dict(1=>"波长λ", 2=>"孔半宽R", 3=>"屏孔间距L", 4=>"正多边形n")]
    formats=[Dict(1=>"{:.2e}m", 2=>"{:.2e}m", 3=>"{:.2f}m", 4=>"{:}边")]
    colors = [:grays]

    fig = Figure(;resolution = (1200, 900),fontsize=30)

    ax, ax1 = intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        formats = formats,
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true,spinewidth=3.0), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true, colorrange=(0,0.05)),
        scale=x->(x),
        slider=:y,
        sliderkwargs=(;linestyle=:dash, color=:white)
    )

    Label(fig[0, :], "正多边形孔夫琅和费衍射", fontsize = 30)
    return fig, dso
end

fig, = regular_polygon_fraunhofer_figure()
fig