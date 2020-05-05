def grid_center(grid):
    return tuple(grid.origin[i] + 0.5 * grid.extent[i] for i in range(3))


def on_xvalue(x):
    return lambda pts: pts[:, 0] == x


def on_yvalue(y):
    return lambda pts: pts[:, 1] == y


def on_zvalue(z):
    return lambda pts: pts[:, 2] == z


def on_xmin(grid):
    return on_xvalue(grid.origin[0])


def on_xmax(grid):
    return on_xvalue(grid.origin[0] + grid.extent[0])


def on_ymin(grid):
    return on_yvalue(grid.origin[1])


def on_ymax(grid):
    return on_yvalue(grid.origin[1] + grid.extent[1])


def on_zmin(grid):
    return on_zvalue(grid.origin[2])


def on_zmax(grid):
    return on_zvalue(grid.origin[2] + grid.extent[2])


def on_vertical_boundaries(grid):
    return (
        lambda pts: on_xmin(grid)(pts)
        | on_xmax(grid)(pts)
        | on_ymin(grid)(pts)
        | on_ymax(grid)(pts)
    )


def on_any_boundary(grid):
    return (
        lambda pts: on_vertical_boundaries(grid)(pts)
        | on_zmin(grid)(pts)
        | on_zmax(grid)(pts)
    )


def vertical_boundaries(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_vertical_boundaries(grid)(vertices)

    return select


def xmin_boundary(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_xmin(grid)(vertices)

    return select


def xmax_boundary(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_xmax(grid)(vertices)

    return select


def ymin_boundary(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_ymin(grid)(vertices)

    return select


def ymax_boundary(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_ymax(grid)(vertices)

    return select


def bottom_boundary(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_zmin(grid)(vertices)

    return select


def top_boundary(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_zmax(grid)(vertices)

    return select


def all_boundaries(simulation, grid):
    def select():
        vertices = simulation.global_vertices()
        return on_any_boundary(grid)(vertices)

    return select
