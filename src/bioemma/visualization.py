from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Mapping


@dataclass(frozen=True)
class VisualizationOptions:
    """Layout options for generated Escher visualizations."""

    scaling_factor: float = 4.0
    axis_epsilon: float = 2.0
    markers_dist: float = 10.0
    metabolite_label_shift: tuple[float, float] = (10.0, 10.0)
    reaction_label_shift: tuple[float, float] = (10.0, 10.0)
    canvas_margin_x: float = 160.0
    canvas_margin_y: float = 160.0
    multimarker_distance_fraction: float = 0.3
    use_constant_multimarker_distance: bool = False
    constant_multimarker_distance: float = 300.0
    axis_offset: float = 20.0


def resolve_visualization_options(
    options: VisualizationOptions | Mapping[str, Any] | None = None,
    **overrides: Any,
) -> VisualizationOptions:
    if options is None:
        resolved = VisualizationOptions()
    elif isinstance(options, VisualizationOptions):
        resolved = options
    elif isinstance(options, Mapping):
        resolved = VisualizationOptions(**dict(options))
    else:
        raise TypeError("visualization_options must be a VisualizationOptions instance or mapping.")

    normalized_overrides = {
        key: _normalize_option_value(key, value)
        for key, value in overrides.items()
        if value is not None
    }
    if not normalized_overrides:
        return resolved
    return replace(resolved, **normalized_overrides)


def _normalize_option_value(key: str, value: Any) -> Any:
    if key in {"metabolite_label_shift", "reaction_label_shift"}:
        return _coerce_pair(value, key)
    if key == "use_constant_multimarker_distance":
        return bool(value)
    return float(value) if isinstance(value, int) else value


def _coerce_pair(value: Any, name: str) -> tuple[float, float]:
    if isinstance(value, str):
        parts = value.split(",")
    else:
        parts = list(value)
    if len(parts) != 2:
        raise ValueError(f"{name} must contain exactly two numbers.")
    return (float(parts[0]), float(parts[1]))
