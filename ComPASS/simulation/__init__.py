# _not_available = object()
# def __getattr__(name):
#     return getattr(Simulation(), name)
#     value = getattr(state, name, _not_available)
#     if value is not _not_available:
#         return value
#     value = getattr(simulation_wrapper, name, _not_available)
#     if value is not _not_available:
#         return value
#     raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
