"""
rimopy version handler module

Handles versioning for the `rimopy` package.
"""
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("rimopy")
except PackageNotFoundError:
    __version__ = "0.0.0"
