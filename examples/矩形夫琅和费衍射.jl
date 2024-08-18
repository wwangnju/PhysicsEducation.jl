using PhysicsEducation
using GLMakie
using Makie.Colors
function rectangular_fraunhofer_figure()
    λ = 600.0e-9 #unit m
    W, H, L = 4e-6, 4e-6,1.0
    t0 = 0.0
    p = [λ, W, H, L]
    xy = Vector{Float64}[collect(-1.0:1e-2:1.0), collect(-1:0.01:1) ]
 
    slit = WaveFunction(rectangular_fraunhofer, t0, xy, p)
    dso = WaveObservable([slit])
    ps = [Dict(1=>4e-7:1e-8:8e-6, 2=>1e-6:1e-6:1e-5, 3=>1e-6:1e-6:1e-5, 4=>0.5:0.1:3.0)]
    pnames = [Dict(1=>"波长λ", 2=>"矩宽W", 3=>"矩高H", 4=>"矩屏间距L")]
    colors = [:grays]

    fig = Figure(;resolution = (1000, 800))

    intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true, colorrange=(0,0.005)),
        scale=x->(x),
        slider=:y,
        sliderkwargs=(;linestyle=:dash, color=:white)
    )
    # Label(fig[0, :], "Rectangular Aperture", fontsize = 30)
    Label(fig[0, :], "矩形光圈衍射", fontsize = 30)
    return fig, dso
end
fig, = rectangular_fraunhofer_figure()
fig