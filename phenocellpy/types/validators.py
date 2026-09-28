"""
Reusable single-field checks for the config dataclasses.

Each check has the signature ``check(name, value)`` and raises on an invalid value. They are attached to dataclass
fields through ``field(metadata={"checks": (...)})`` and run by :class:`phenocellpy.types.base.Validated`.

All checks except :func:`not_none` accept ``None``, as ``None`` means "use the model default" for most optional
parameters. Stack :func:`not_none` in front of a check to make the field required.
"""

from numbers import Real


def not_none(name, value):
    if value is None:
        raise ValueError(f"'{name}' must be defined. Got None.")


def _number(name, value):
    # bool is a subclass of int, but a flag passed as a number is a bug
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"'{name}' must be a number. Got {type(value).__name__}: {value!r}.")


def positive(name, value):
    if value is None:
        return
    _number(name, value)
    if value <= 0:
        raise ValueError(f"'{name}' must be greater than 0. Got {value}.")


def non_negative(name, value):
    if value is None:
        return
    _number(name, value)
    if value < 0:
        raise ValueError(f"'{name}' must be >= 0. Got {value}.")


def fraction(name, value):
    if value is None:
        return
    _number(name, value)
    if not 0 <= value <= 1:
        raise ValueError(f"'{name}' must be in range [0, 1]. Got {value}.")


def boolean(name, value):
    if value is None:
        return
    if not isinstance(value, bool):
        raise TypeError(f"'{name}' must be a bool. Got {type(value).__name__}: {value!r}.")


def integer(name, value):
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"'{name}' must be an int. Got {type(value).__name__}: {value!r}.")


def non_empty_str(name, value):
    if value is None:
        return
    if not isinstance(value, str) or not value:
        raise ValueError(f"'{name}' must be a non-empty string. Got {value!r}.")
