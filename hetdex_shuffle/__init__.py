from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pkg_resources
__version__ = pkg_resources.get_distribution('hetdex-shuffle').version

__svn_revision__ = '$Revision: 228 $'

__full_version__ = '''Shuffle version: {}.
Svn revision: {}'''.format(__version__, __svn_revision__)
