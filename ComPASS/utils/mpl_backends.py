# The following implements the trick provided by:
# https://stackoverflow.com/questions/3285193/how-to-switch-backends-in-matplotlib-python
# to access pyplot with an available backend

preferred_backends = ["TkAgg", "Qt4Agg", "Agg"]


def import_pyplot(raise_exception=True):
    for backend in preferred_backends:
        try:
            import matplotlib

            matplotlib.use(backend, warn=False, force=True)
            from matplotlib import pyplot

            break
        except:
            continue
    else:
        if raise_exception:
            raise Exception("No suitable matplotlib backend found!")
        pyplot = None
    return pyplot
