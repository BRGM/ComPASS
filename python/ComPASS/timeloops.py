import ComPASS
from ComPASS.utils.units import day, year
mpi = ComPASS.mpi

def check_well_pressure():
    p_min, p_max = float('Inf'), -float('Inf')
    for state in [ComPASS.node_states(), ComPASS.fracture_states(), ComPASS.cell_states()]:
        if state.p.shape[0]>0:
            p_min = min(p_min, state.p.min())
            p_max = max(p_max, state.p.max())
    # whp = well head pressure
    ComPASS.production_whp()[:] = p_min - 1.
    ComPASS.injection_whp()[:] = p_max + 1.

def standard_loop(final_time, initial_timestep=1., output_frequency = None, nb_output = 10, nitermax = None, tstart=0):
    if output_frequency is None:
        nb_output = max(2, nb_output)
        output_frequency = (max(tstart, final_time) - tstart) / (nb_output - 1)
    assert output_frequency is not None and output_frequency>0
    # this is necessary for well operating on pressures
    check_well_pressure()
    t = tstart
    timestep = initial_timestep
    n = 0
    t_output = 0
    @mpi.on_master_proc
    def print_iteration_info():
        print()
        print('Time Step (iteration):', n)
        print('Current time: %.1f y' % (t / year), ' -> final time:', final_time / year, 'y')
        print('Timestep: %.3g s = %.3f d = %.3f y' % (timestep, timestep / day, timestep / year))
    while t <= final_time and (nitermax is None or n < nitermax):
        if t >= t_output:
            ComPASS.output_visualization_files(n)
            # WARNING / CHECKME we may loose some outputs
            while (t_output < t):
                t_output = t_output + output_frequency
        n+= 1
        print_iteration_info()
        ComPASS.make_timestep(timestep)
        t = ComPASS.get_current_time()
        timestep = ComPASS.get_timestep()
        ComPASS.timestep_summary()
    # Output final time
    ComPASS.output_visualization_files(n)
    mpi.synchronize()