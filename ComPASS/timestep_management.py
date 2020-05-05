
class FixedTimeStep:

    def __init__(self, step):
        self.fixed_step = step
        self.step = step
    
    @property
    def current_step(self):
        return self.step
    
    def steps(self, upper_bound=None):
        self.step = self.fixed_step
        if upper_bound is not None:
            self.step = min(self.step, upper_bound)
        return (self.step,)


class TimeStepManager:

    def __init__(self,
            initial_timestep, maximum_timestep=None,
            increase_factor=1.2, decrease_factor=0.8,
            minimum_timestep=1E-14, # CHECKME: machine precision, should not happen when solving 
    ):
        self.previous = None
        self.step = initial_timestep
        assert self.current_step > 0
        self.minimum = minimum_timestep
        self.maximum = maximum_timestep
        assert increase_factor > 1
        self.increase_factor = increase_factor
        assert 0 < decrease_factor < 1
        self.decrease_factor = decrease_factor

    @property
    def current_step(self):
        return self.step

    @current_step.setter
    def current_step(self, step):
        self.previous = None
        self.step = step

    def steps(self, upper_bound=None):
        if self.previous:
            self.previous*= self.increase_factor
            self.step = self.previous
            if self.maximum is not None:
                self.step = min(self.step, self.maximum)
        self.previous = self.step
        if upper_bound is not None:
            self.step = min(self.step, upper_bound)
        def attempts():
            while self.current_step > self.minimum:
                yield self.step
                self.step*= self.decrease_factor
                self.previous = self.step
        return attempts()
        

if __name__=='__main__':
    print(list(FixedTimeStep(1).steps()))
    tsm = TimeStepManager(1, decrease_factor=0.1)
    print(list(tsm.steps()))
    tsm.current_step = 2E-6
    print(list(tsm.steps()))
