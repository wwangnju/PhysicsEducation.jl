using PhysicsEducation
using Test
using GLMakie
# using Makie.Colors
@testset "PhysicsEducation.jl" begin

function lissajous_figure()
    Δt, tail = 0.01, 5000

    A₁,  θ₁ = 2., 0.0
    A₂, θ₂ = 2., 0.0
    ω₁, ω₂ = 0.2, 0.2

    t0, x0 = 0.0, 5
    
    p = [ω₁, A₁, θ₁, ω₂, A₂, θ₂]
    
    lsa = DynamicalFunction(lissajous, t0, p)
    
    fig = Figure(;resolution = (1000, 1000))
    
    dso = TrajectoryObservable([lsa]; tail=tail, current_step=0, Δt=Δt)
    
    ps = [Dict(1=>0.1:0.1:1, 2=>1.0:1:4.0, 3=> 0.0:0.05:2pi, 4=>0.1:0.1:1, 5=>1.0:1:4.0, 6=> 0.0:0.05:2pi)]
    pnames = [Dict(1=>(L"f_x"), 2=>(L"A_x"), 3=>L"θ_x", 4=>L"f_y", 5=>L"A_y", 6=>L"θ_y")]
    colors = [:magenta, ]
    trajectory_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(;aspect=1.0, limits=((-x0, x0), (-x0, x0)), xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true)
    )
    timeserieslayout = fig[:, 2] = GridLayout()
    fs = [1, 2]

    timeseries_plots!(
        timeserieslayout, dso, fs; 
        colors=colors, 
        timeseries_names=["x/m", "y/m"], 
        linekwargs=(linewidth=3.0, label=nothing),
        timeseries_ylims=[(-x0, x0), (-x0, x0)], 
        axis=(;xminorgridvisible=true, yminorgridvisible=true),
        timeunit = 1,
        timelabel="时间/s"
    )
    Label(fig[0, :], "李萨如图形", fontsize = 30)
    return fig, dso
end
fig, =lissajous_figure()
save("lissajous.png", fig)

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
    colors = [:viridis]

    fig = Figure(;size = (1000, 1000), fontsize=15)
    intensity_plot!(fig, dso; 
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
fig1, = slit_grating_figure()
save("slit_grating.png", fig1)
end
