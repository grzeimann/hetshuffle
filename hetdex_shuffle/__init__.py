from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from importlib.metadata import version

__version__ = version('hetdex-shuffle')

__svn_revision__ = '$Revision: 228 $'

__full_version__ = '''Shuffle version: {}.
Svn revision: {}'''.format(__version__, __svn_revision__)
