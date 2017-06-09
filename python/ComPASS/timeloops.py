import ComPASS
from ComPASS.utils.units import day, year
mpi = ComPASS.mpi

def standard_loop(final_time, output_frequency = None, nb_output = 10, nitermax = None, tstart=0):
    if output_frequency is None:
        nb_output = max(2, nb_output)
        output_frequency = (max(tstart, final_time) - tstart) / (nb_output - 1)
    assert output_frequency is not None and output_frequency>0
    t = tstart
    n = 0
    t_output = 0
    @mpi.on_master_proc
    def print_iteration_info():
        print()
        print('Time Step (iteration):', n)
        print('Current time: %.1f years' % (t / year), ' -> final time:', final_time / year)
        dt = ComPASS.get_timestep()
        print('Timestep: %.3f days = %.3f years' % (dt / day, dt / year))
    while t <= final_time and (nitermax is None or n < nitermax):
        if t >= t_output:
            ComPASS.output_visualization_files(n)
            # WARNING / CHECKME we may loose some outputs
            while (t_output < t):
                t_output = t_output + output_frequency
        n+= 1
        print_iteration_info()
        ComPASS.make_timestep()
        t = ComPASS.get_current_time()
        ComPASS.timestep_summary()
    # Output final time
    ComPASS.output_visualization_files(n)
    mpi.synchronize()