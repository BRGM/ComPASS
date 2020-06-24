import os
import numpy as np


class RasterHeader:
    def __init__(self, filename):
        with open(filename) as f:
            for _ in range(6):
                l = f.readline().split()
                name, value = l[0].strip().lower(), l[1].strip()
                if name in ["ncols", "nrows"]:
                    value = int(value)
                elif name in [
                    "cellsize",
                    "xllcorner",
                    "yllcorner",
                    "xllcenter",
                    "yllcenter",
                ]:
                    value = value.replace(",", ".")
                    value = float(value)
                elif name != "nodata_value":
                    print("Unknown name:", name)
                    assert False
                setattr(self, name, value)


class RasterInfo:
    def __init__(self, filename):
        self.source = filename
        self.header = None
        if os.path.exists(filename):
            self.header = RasterHeader(filename)

    @property
    def llcorner(self):
        header = self.header
        if hasattr(header, "xllcorner"):
            return header.xllcorner, header.yllcorner
        assert hasattr(header, "xllcenter")
        offset = 0.5 * header.cellsize
        return header.xllcenter - offset, header.yllcenter - offset

    @property
    def urcorner(self):
        header = self.header
        Lx = header.cellsize * header.ncols
        Ly = header.cellsize * header.nrows
        if hasattr(header, "xllcorner"):
            return header.xllcorner + Lx, header.yllcorner + Ly
        assert hasattr(header, "xllcenter")
        offset = 0.5 * header.cellsize
        return header.xllcenter - offset + Lx, header.yllcenter - offset + Ly

    @property
    def llcenter(self):
        header = self.header
        if hasattr(header, "xllcenter"):
            return header.xllcenter, header.yllcenter
        assert hasattr(header, "xllcorner")
        offset = 0.5 * header.cellsize
        return header.xllcorner + offset, header.yllcorner + offset

    @property
    def urcenter(self):
        header = self.header
        Lx = header.cellsize * (header.ncols - 1)
        Ly = header.cellsize * (header.nrows - 1)
        if hasattr(header, "xllcenter"):
            return header.xllcenter + Lx, header.yllcenter + Ly
        assert hasattr(header, "xllcorner")
        offset = 0.5 * header.cellsize
        return header.xllcorner + offset + Lx, header.yllcorner + offset + Ly

    @property
    def nodata(self):
        return self.header.nodata_value

    def centers(self):
        header = self.header
        nx = header.ncols
        ny = header.nrows
        dx = header.cellsize
        Lx = dx * nx
        Ly = dx * ny
        O = self.llcenter
        return np.meshgrid(
            np.arange(O[0], O[0] + Lx - 0.5 * dx, dx),
            np.arange(O[1], O[1] + Ly - 0.5 * dx, dx)[::-1],
        )

    def intersect(self, box):
        bllc, burc = box
        xmin, ymin = self.llcorner
        xmax, ymax = self.urcorner
        if xmax < bllc[0]:
            return False
        if ymax < bllc[1]:
            return False
        if xmin > burc[0]:
            return False
        if ymin > burc[1]:
            return False
        return True
