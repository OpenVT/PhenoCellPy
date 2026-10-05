from dataclasses import fields, replace
from typing import Dict

class Validated:
    """
    Mixin for the config dataclasses. Every assignment to a field, including the ones done by the generated
    ``__init__``, runs the single-field checks declared in the field's ``metadata["checks"]``. After construction,
    :meth:`validate` runs the checks that need more than one field.

    Subclasses override :meth:`validate`, never ``__post_init__``, so the cross-field checks always run.
    """

    def __setattr__(self, name, value):
        field = self.__dataclass_fields__.get(name)
        if field is not None:
            for check in field.metadata.get("checks", ()):
                check(name, value)
        super().__setattr__(name, value)

    def __post_init__(self):
        self.validate()

    def validate(self):
        """Cross-field checks. Override in subclasses."""
        pass

    def fill_none_from(self, defaults):
        """
        Returns a copy of this config where the fields that are None take the value of the same field in `defaults`.
        Nested configs are filled field by field.

        :param defaults: Config of the same class to take the missing values from
        """
        changes = {}
        for f in fields(self):
            value, default = getattr(self, f.name), getattr(defaults, f.name)
            if value is None:
                changes[f.name] = default
            elif isinstance(value, Validated) and isinstance(default, Validated):
                changes[f.name] = value.fill_none_from(default)
        return replace(self, **changes)

    @classmethod
    def _reject_unknown_keys(cls, data: Dict, known=None):
        """
        Raises a readable error for keys in `data` that are not fields of `cls` (e.g., a typo in a JSON file),
        instead of the ``TypeError`` from the generated ``__init__``.
        """
        if not isinstance(data, dict):
            raise TypeError(f"{cls.__name__} expects a dict. Got {type(data).__name__}: {data!r}.")
        if known is None:
            known = {f.name for f in fields(cls)}
        unknown = set(data) - set(known)
        if unknown:
            raise ValueError(f"Unknown key(s) for {cls.__name__}: {sorted(unknown)}. Valid keys: {sorted(known)}.")
