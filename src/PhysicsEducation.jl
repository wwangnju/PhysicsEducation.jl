module PhysicsEducation

export set_parameter!, trajectory_plot!, current_state, current_parameters, step!, set_state!, reinit!

function set_parameter! end
function trajectory_plot! end
function current_state end
function current_parameters end
function step! end
function set_state! end
function reinit! end

# abstract type AbstractStateObservable end


include("Vibration.jl")
using .Vibration
export DynamicalFunction, TrajectoryObservable
export simpleharmonic, timeseries_plots!, superposition_plot!, lissajous, dampedvibration, forcevibration

include("Wave.jl")
using .Wave
export harmonicwave1D, harmonicwave2D
export WaveObservable, WaveFunction, intensity_plot!, doubleslit_intensity, michelsoninterferometer_intensity
export doubleslit_fraunhofer, singleslit_fraunhofer, rectangular_fraunhofer, circular_fraunhofer, slit_grating, newton_ring_intensity

end
