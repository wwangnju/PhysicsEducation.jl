module PhysicsEducation

include("SimpleHarmonicMotion.jl")
using .SimpleHarmonicMotion
export DynamicalFunction, current_state, current_parameters, TrajectoryObservable, step!, set_state!, reinit!
export set_parameter!, harmonicmotion_plot!, simpleharmonic, timeseries_plots!



end
