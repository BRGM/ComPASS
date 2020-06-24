# -*- coding: utf-8 -*-
"""
Useful addition for teaching algorithmics.
"""
from __future__ import unicode_literals

from keyword import iskeyword
from numpy import array, empty
from re import compile as re_compile
from sys import version_info

_PUBLIC_NAME_RE = re_compile(r"(?![0-9_])\w+$")


def _is_public_name(name):
    """Return True is name is a public attribute name"""
    return _PUBLIC_NAME_RE.match(name) and not iskeyword(name)


class Struct(object):
    """A super-class for struct-like classes

    """

    def __setattr__(self, name, val):
        """Prevent new attributes"""
        dtype = self.__class__.__dict__.get(name)
        if _is_public_name(name):
            if dtype is None:
                raise AttributeError("Can not add attribute {} to Struct".format(name))
            elif not isinstance(val, dtype):
                if dtype == float and isinstance(val, int):
                    val = float(val)
                else:
                    raise TypeError(
                        "Attribute {} of Struct {} must be {}".format(
                            name, self.__class__.__name__, dtype.__name__
                        )
                    )
        super(Struct, self).__setattr__(name, val)

    def __init__(self, **kw):
        for key, val in kw.items():
            setattr(self, key, val)
        missing = set()
        for keycls in self.__class__.__dict__:
            if _is_public_name(keycls) and keycls not in self.__dict__:
                missing.add(keycls)
        if missing:
            missing = ", ".join(missing)
            raise ValueError(
                "Missing attribute(s) {} in Struct {}".format(
                    missing, self.__class__.__name__
                )
            )

    def __repr__(self):
        fields = []
        for key, val in self.__dict__.items():
            if _is_public_name(key):
                fields.append("%s=%r" % (key, val))
        fields = ", ".join(fields)
        return "%s(%s)" % (self.__class__.__name__, fields)
        self.__class__.__name__


######## TESTS


def test_is_public_name():
    """Testing _is_public_name"""
    assert _is_public_name("foo")
    assert _is_public_name("foo_123")
    assert not _is_public_name("_private")
    assert not _is_public_name("foo bar")
    assert not _is_public_name("foo-bar")
    assert not _is_public_name("1foobar")
    assert not _is_public_name("while")

    def only_from_py3(test):
        """Test should pass only from Python 3.x"""
        if version_info[0] >= 3:
            return test
        else:
            return not test

    assert only_from_py3(_is_public_name("hé"))
    assert only_from_py3(_is_public_name("é"))


class Point(Struct):
    x = float
    y = float


def test_struct_empty():
    s1 = Struct()
    s1._private = 42
    try:
        s1.foo = 42
        assert False, "AttributeError expected"
    except AttributeError:
        pass


def test_corrrect_init():
    p = Point(x=0.0, y=0.0)


def test_missing_init_attribute():
    class Point(Struct):
        x = float
        y = float

    try:
        p = Point(x=0.0)
        assert False, "ValueError expected"
    except ValueError:
        pass


def test_spurious_init_attribute():
    try:
        p = Point(x=0.0, y=0.0, z=0.0)
        assert False, "AttributeError expected"
    except AttributeError:
        pass


def test_corrrect_attributes():
    p = Point(x=0.0, y=0.0)
    p._private = "FOO"
    p.x = 1.2


def test_incorrect_attribute():
    p = Point(x=0.0, y=0.0)
    try:
        p.z = 0.0
        assert False, "AttributeError expected"
    except AttributeError:
        pass


def test_incorrect_type():
    p = Point(x=0.0, y=0.0)
    try:
        p.x = "foo"
        assert False, "TypeError expected"
    except TypeError:
        pass


def test_int_cast_to_float():
    p = Point(x=0.0, y=0.0)
    p.x = 1
