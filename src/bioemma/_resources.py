from importlib.resources import files
from pathlib import Path


def resource_path(name: str) -> Path:
    return files("bioemma").joinpath("resources", name)
