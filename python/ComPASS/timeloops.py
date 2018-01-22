#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
from ComPASS.utils.units import day, year
mpi = ComPASS.mpi
Dumper = ComPASS.dumps.Dumper

def check_well_pressure():
    p_min, p_max = float('Inf'), -float('Inf')
    for state in [ComPASS.node_states(), ComPASS.fracture_states(), ComPASS.cell_states()]:
        if state.p.shape[0]>0:
            p_min = min(p_min, state.p.min())
            p_max = max(p_max, state.p.max())
    # whp = well head pressure
    ComPASS.production_whp()[:] = p_min - 1.
    ComPASS.injection_whp()[:] = p_max + 1.

class Snapshooter:

    def __init__(self, dumper=None):
        if dumper is None:
            dumper = Dumper()
        dumper.start_simulation()
        self.dumper = dumper
        self.nb_snaphots = 0
        self.latest_snapshot_time = None
        self.filename = self.dumper.to_output_directory('snapshots')
        if mpi.is_on_master_proc:
            with open(self.filename, 'w') as f:
                pass

    def new_tag(self):
        tag = '%05d' % self.nb_snaphots
        self.nb_snaphots+= 1
        return tag

    def shoot(self, t):
        tag = self.new_tag()
        assert self.latest_snapshot_time is None or t > self.latest_snapshot_time
        self.latest_snapshot_time = t
        if mpi.is_on_master_proc:
            with open(self.filename, 'a') as f:
                print(tag, '%.12g'%t, file=f)
        self.dumper.dump_states(tag)


def standard_loop(final_time, initial_timestep=1.,
                  output_period = None, output_every = None,
                  nb_output = None, nitermax = None, tstart=0, dumper=None,
                  iteration_callbacks = None, output_callbacks = None,
                 ):
    if output_period is None:
        if nb_output is None:
            output_period = final_time    
        else:
            nb_output = max(2, nb_output)
            output_period = (max(tstart, final_time) - tstart) / (nb_output - 1)
    assert not(output_period is None or output_period<=0)
    if iteration_callbacks is None:
        iteration_callbacks = tuple()
    if output_callbacks is None:
        output_callbacks = tuple()
    # this is necessary for well operating on pressures
    check_well_pressure()
    t = tstart
    timestep = initial_timestep
    n = 0
    t_output = 0
    shooter = Snapshooter(dumper)
    @mpi.on_master_proc
    def print_iteration_info():
        print()
        print('Time Step (iteration):', n)
        print('Current time: %.1f y' % (t / year), ' -> final time:', final_time / year, 'y')
        print('Timestep: %.3g s = %.3f d = %.3f y' % (timestep, timestep / day, timestep / year))
    while t < final_time and (nitermax is None or n < nitermax):
        if t >= t_output or not(output_every is None or n%output_every>0):
            for callback in output_callbacks:
                callback(n, t)
            shooter.shoot(t)
            # WARNING / CHECKME we may loose some outputs
            while (t_output < t):
                t_output = t_output + output_period
        n+= 1
        print_iteration_info()
        ComPASS.make_timestep(timestep)
        t = ComPASS.get_current_time()
        timestep = ComPASS.get_timestep()
        ComPASS.timestep_summary()
        for callback in iteration_callbacks:
            callback(n, t)
    # Output final time
    if shooter.latest_snapshot_time < t:
        for callback in output_callbacks:
            callback(n, t)
        shooter.shoot(t)
    mpi.synchronize()
