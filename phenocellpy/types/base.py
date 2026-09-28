from dataclasses import fields


class Validated:
    """
    Mixin for the config dataclasses. On construction it runs the single-field checks declared in each field's
    ``metadata["checks"]``, then calls :meth:`validate` for checks that need more than one field.

    Subclasses override :meth:`validate`, never ``__post_init__``, so the field checks always run.
    """

    def __post_init__(self):
        for f in fields(self):
            for check in f.metadata.get("checks", ()):
                check(f.name, getattr(self, f.name))
        self.validate()

    def validate(self):
        """Cross-field checks. Override in subclasses."""
        pass

    @classmethod
    def _reject_unknown_keys(cls, data: dict, known=None):
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
