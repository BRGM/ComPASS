from collections import namedtuple
import re


def parametric_functor(*params, **documented_params):
    assert len(params) == 0 or len(documented_params) == 0
    documentation = []
    defaults = []
    for name, info in documented_params.items():
        if type(info) is str:
            assert len(defaults) == 0
            documentation.append((name, info))
        elif type(info) is tuple:
            doc, default = info
            documentation.append((name, doc))
            defaults.append(default)
        else:
            documentation.append((name, ""))
            defaults.append(info)
    assert len(documentation) == len(documentation)
    if len(params) == 0:
        assert len(documentation) > 1
        params = [name for name, _ in documentation]
    else:
        assert len(params) == 1
        params = params[0]
    lines = []
    offset = len(documentation) - len(defaults)
    for k, (name, doc) in enumerate(documentation):
        l = f":param {name}: {doc}"
        if k >= offset:
            l = f"{l}, defaults to {str(defaults[k - offset])}"
        lines.append(l)
    documentation = "\n".join(lines)

    def append_documentation(functor_documentation):
        if functor_documentation is None:
            if len(documentation) == 0:
                return None
            else:
                return documentation
        if len(documentation) == 0:
            return functor_documentation
        return "\n".join([functor_documentation, documentation])

    def decorator(functor):
        class decorated_functor(namedtuple("params", params, defaults=defaults)):
            def __call__(self, *args, **kwargs):
                return functor.__call__(self, *args, **kwargs)

            def __repr__(self):
                parameter_values = re.match(
                    f"{self.__class__.__name__}(.*)", super().__repr__()
                ).group(1)
                return f"{functor.__name__}{parameter_values}"

        decorated_functor.__doc__ = append_documentation(functor.__doc__)
        call_documentation = f"Calls the underlying {functor.__name__} functor."
        if functor.__call__.__doc__ is not None:
            call_documentation = "\n".join(
                [call_documentation, functor.__call__.__doc__]
            )
        decorated_functor.__call__.__doc__ = call_documentation
        return decorated_functor

    return decorator


if __name__ == "__main__":

    @parametric_functor(
        a="slope", b=("intersection with ordinates axis", 0),
    )
    class Line2D:
        """Line equation"""

        def __call__(self, x):
            """The line equation evaluation."""
            return self.slope * x + self.b
