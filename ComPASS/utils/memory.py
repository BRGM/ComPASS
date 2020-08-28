import os
import resource
import psutil

# adapted from: https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
def sizeof_fmt(num, suffix="b"):
    for unit in ["", "k", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, "Yi", suffix)


class MemStatus:
    def __init__(self, context=None, skip_first=False):
        self.process = psutil.Process(os.getpid())
        self.context = context
        self.latest_rss = self.rss
        self.latest_peak = self.peak
        self.rss_increase = None if skip_first else 0
        self.peak_increase = None if skip_first else 0

    @property
    def rss(self):
        return self.process.memory_info().rss

    @property
    def rss_str(self):
        return sizeof_fmt(self.rss)

    @property
    def peak(self):
        return (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
        )  # FIXME: Linux only

    @property
    def peak_str(self):
        return sizeof_fmt(self.peak)

    @property
    def context_str(self):
        if self.context is None:
            return "Memory usage:"
        return f"Memory usage ({self.context}):"

    def reset(self, verbose=False):
        rss = self.rss
        if verbose and rss != self.latest_rss:
            print(
                f"{self.context_str} reset rss with delta:",
                f"{sizeof_fmt(rss - self.latest_rss)}",
            )
        self.latest_rss = rss
        peak = self.peak
        if verbose and peak != self.latest_peak:
            print(
                f"{self.context_str} reset peak with delta:",
                f"{sizeof_fmt(peak - self.latest_peak)}",
            )
        self.latest_peak = peak

    def update(self, verbose=False):
        rss = self.rss
        if self.rss_increase is None:
            self.rss_increase = 0
        else:
            self.rss_increase += rss - self.latest_rss
        if verbose:
            print(
                f"{self.context_str} current rss: {sizeof_fmt(rss)}",
                f"delta: {sizeof_fmt(rss - self.latest_rss)}",
                f"(total increase: {sizeof_fmt(self.rss_increase)})",
            )
        self.latest_rss = rss
        peak = self.peak
        if self.peak_increase is None:
            self.peak_increase = 0
        else:
            assert peak >= self.latest_peak
            self.peak_increase += peak - self.latest_peak
        if verbose:
            print(
                f"{self.context_str} current peak: {sizeof_fmt(peak)}",
                f"delta: {sizeof_fmt(peak - self.latest_peak)}",
                f"(total increase: {sizeof_fmt(self.peak_increase)})",
            )
        self.latest_peak = peak

    def __str__(self):
        return (
            f"{self.context_str} rss: {self.rss_str} "
            + f"rss increase: {sizeof_fmt(self.rss_increase)}\n"
            + f"{self.context_str} peak: {self.peak_str} "
            + f"rss increase: {sizeof_fmt(self.peak_increase)}"
        )
