try:
    from .version import version as __version__

    __all__ = ["__version__"]
except ImportError:
    __version__ = ''
