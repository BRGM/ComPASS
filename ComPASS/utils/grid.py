def grid_center(grid):
    return tuple(grid.origin[i] + 0.5 * grid.extent[i] for i in range(3))


def below(axis, lim):
    return lambda pts: pts[:, axis] <= lim


def above(axis, lim):
    return lambda pts: pts[:, axis] >= lim


def on_xvalue(x):
    return lambda pts: pts[:, 0] == x


def on_yvalue(y):
    return lambda pts: pts[:, 1] == y


def on_zvalue(z):
    return lambda pts: pts[:, 2] == z


def below_axis_limit(grid, axis, epsilon=0):
    return below(axis, grid.origin[axis] + epsilon)


def above_axis_limit(grid, axis, epsilon=0):
    return above(axis, grid.origin[axis] + grid.extent[axis] - epsilon)


def on_xmin(grid, epsilon=0):
    return below_axis_limit(grid, 0, epsilon)


def on_xmax(grid, epsilon=0):
    return above_axis_limit(grid, 0, epsilon)


def on_ymin(grid, epsilon=0):
    return below_axis_limit(grid, 1, epsilon)


def on_ymax(grid, epsilon=0):
    return above_axis_limit(grid, 1, epsilon)


def on_zmin(grid, epsilon=0):
    return below_axis_limit(grid, 2, epsilon)


def on_zmax(grid, epsilon=0):
    return above_axis_limit(grid, 2, epsilon)


def on_vertical_boundaries(grid, epsilon=0):
    return (
        lambda pts: on_xmin(grid, epsilon)(pts)
        | on_xmax(grid, epsilon)(pts)
        | on_ymin(grid, epsilon)(pts)
        | on_ymax(grid, epsilon)(pts)
    )


def on_any_boundary(grid, epsilon=0):
    return (
        lambda pts: on_vertical_boundaries(grid, epsilon)(pts)
        | on_zmin(grid, epsilon)(pts)
        | on_zmax(grid, epsilon)(pts)
    )


def vertical_boundaries(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_vertical_boundaries(grid, epsilon)(vertices)

    return select


def xmin_boundary(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_xmin(grid, epsilon)(vertices)

    return select


def xmax_boundary(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_xmax(grid, epsilon)(vertices)

    return select


def ymin_boundary(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_ymin(grid, epsilon)(vertices)

    return select


def ymax_boundary(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_ymax(grid, epsilon)(vertices)

    return select


def bottom_boundary(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_zmin(grid, epsilon)(vertices)

    return select


def top_boundary(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_zmax(grid, epsilon)(vertices)

    return select


def all_boundaries(simulation, grid, epsilon=0):
    def select():
        vertices = simulation.global_vertices()
        return on_any_boundary(grid, epsilon)(vertices)

    return select
