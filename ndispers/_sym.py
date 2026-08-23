"""
Lazy access to sympy.

The dispersion formulas are written as sympy expressions, but the functions
the package actually evaluates are generated from them ahead of time
(``ndispers/_compiled``, built by ``tools/compile_media.py``). sympy is thus
only needed to look at an expression symbolically, to differentiate one, or
to evaluate a medium that has no compiled module (a user subclass, say). It
is imported the first time any of that happens, and ``pip install
ndispers[sym]`` provides it.
"""


def _sympy():
    try:
        import sympy
    except ImportError as e:          # pragma: no cover
        raise ImportError(
            "sympy is needed for symbolic access to the dispersion formulas "
            "(n_expr and friends, or a medium without a compiled module): "
            "pip install ndispers[sym]") from e
    return sympy


class _LazyModule:
    """Stands in for the sympy module; imports it on first attribute access."""
    __slots__ = ()

    def __getattr__(self, name):
        return getattr(_sympy(), name)


sympy = _LazyModule()


class _LazySymbol:
    """A sympy Symbol that is only created when it is used in an expression."""
    __slots__ = ("name",)

    def __init__(self, name):
        self.name = name

    def _sympy_(self):                # sympify() hook: lets sympy functions accept it
        return _sympy().Symbol(self.name)

    def __repr__(self):
        return self.name

    def __pos__(self):           return self._sympy_()
    def __neg__(self):           return -self._sympy_()
    def __add__(self, o):        return self._sympy_() + o
    def __radd__(self, o):       return o + self._sympy_()
    def __sub__(self, o):        return self._sympy_() - o
    def __rsub__(self, o):       return o - self._sympy_()
    def __mul__(self, o):        return self._sympy_() * o
    def __rmul__(self, o):       return o * self._sympy_()
    def __truediv__(self, o):    return self._sympy_() / o
    def __rtruediv__(self, o):   return o / self._sympy_()
    def __pow__(self, o):        return self._sympy_() ** o
    def __rpow__(self, o):       return o ** self._sympy_()


wl = _LazySymbol('lambda')
phi = _LazySymbol('phi')
theta = _LazySymbol('theta')
T = _LazySymbol('T')
