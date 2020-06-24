import time

from .. import messages
from .. import mpi


class Timer:
    def __init__(
        self,
        collect_events=True,
        output=None,
        output_on_master=True,
        banner=None,
        fieldwidth=20,
    ):
        """
        Timer instance constructor.

        :param collect_events: boolean that tells if timed event should be collected
        :param output: a stream where timed events can be printed defaults to None (no messages)
        :param output_on_master: boolean that tells if messages are output only on master proc (if False
         the proc id is appended to the banner
        :param banner: a string that can be prepended to the output messages, defaults to None
        :param fieldwidth: the width used to print event name
        """
        if not (collect_events or output):
            message.warning(
                """Your timer is probably useless
            as it neither outputs timed event nor collect them!"""
            )
        self.start = time.process_time_ns()
        self.events = [] if collect_events else None
        self.fieldwidth = fieldwidth
        self.output = output
        self.dump = output is not None and (
            not output_on_master or mpi.is_on_master_proc
        )
        if banner is None:
            banner = ""
        else:
            banner = banner + " > "
        self.banner = banner + (
            "" if output_on_master else "proc % 5d > " % mpi.proc_rank
        )

    def _format_event(self, name, t):
        banner = self.banner
        w = self.fieldwidth
        return "%s%s: %.1fs" % (banner, name.ljust(w)[:w], t)

    def __call__(self, name="no name"):
        self.event(name)

    def event(self, name="no name"):
        now = 1e-9 * (time.process_time_ns() - self.start)
        if self.events is not None:
            self.events.append((name, now))
        if self.dump:
            assert self.output is not None
            print(self._format_event(name, now), file=self.output)

    def __repr__(self):
        if self.events is not None:
            return "\n".join([self._format_event(name, t) for name, t in self.events])


if __name__ == "__main__":
    import sys

    def do_some_work(n=5000):
        for i in range(n):
            for j in range(n):
                pass

    # default timer will ouput on master proc without banner
    t1 = Timer()
    # a timer can have a banner that is prepended to the event id
    t2 = Timer(output=sys.stdout, output_on_master=False, banner="T2")
    t1("an event")
    do_some_work()
    t1("another event")
    t2("another event")
    if mpi.is_on_master_proc:
        do_some_work()
    t1("yet another event")
    do_some_work()
    print()
    print("All events recorded by t1 (on proc %d):" % mpi.proc_rank)
    print(t1)
