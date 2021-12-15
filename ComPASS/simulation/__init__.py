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

from enum import Enum, auto


class AlignmentMethod(Enum):
    """
    Used in Jacobian.F90
    The idea is to have postive diagonal using linear combinations
    (alternative is to used inverse of block = LC of)
    good for LU O (pas bonne pour amg)
    not used if not preconditionner (but avoid pivoting)
    For the manual alignment method
    it is necessary to give a three-dimension matrix: aligmat in DefModel
    aligmat(:,:,ic) is the alignment matrix for context ic
    the index order of aligmat(:,:,ic) is (col,row), it allows us
    to define aligmat(:,:,ic) without considering that the matrix
    in Fortran is column-major
    """

    manual = auto()
    inverse_diagonal = auto()
