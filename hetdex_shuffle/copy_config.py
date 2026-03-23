'''Entry point to copy the configuration files in the directory from where is
run'''

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse as ap
import io
import os
import textwrap as tw

import pkg_resources
import six

from . import __full_version__

try:  # python 3 doesn't have raw_input
    input = raw_input
except NameError:
    pass

COPY_FILES = ['fplane.txt', 'shuffle.cfg']


def main(argv=None):
    '''Entry point

    Parameters
    ----------
    argv : list of strings
        command line
    '''
    args = parse(argv=argv)
    copy(args)


def parse(argv=None):
    """Create the parser and parse the command line arguments

    Parameters
    ----------
    argv : list
        command line, if None taken from ``sys.argv``

    Returns
    -------
    Namespace
        parsed command line
    """
    # shared options

    # main parser
    parser = ap.ArgumentParser(description="""Copy the shuffle configuration
                               files in the current directory""",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--verbose", '-v', action="store_true",
                        help="Increase verbosity")
    parser.add_argument('--version', '-V', action='version',
                        version=__full_version__)

    parser.add_argument('-t', '--to-dir', default='.',
                        help="Directory where to copy the data")
    parser.add_argument('-f', '--force', action='store_true',
                        help="""Force overwriting existing configuration files.
                        If this option is used, all local modifications will be
                        lost""")

    args = parser.parse_args(args=argv)

    return args


def _str_to_unicode(str_):
    '''In python 2 convert the input string to unicode, otherwise returns the
    input'''
    if six.PY2 and isinstance(str_, six.string_types):
        return unicode(str_)
    else:
        return str_


def copy(args):
    non_written_files, written_files = [], []
    for cf in COPY_FILES:
        ofile = os.path.join(args.to_dir, cf)
        if os.path.exists(ofile) and not args.force:
            overwrite = False
            msg = "Do you want to overwrite file '{}' (y/[n])? ".format(cf)
            while True:
                try:
                    answer = input(msg)
                    if not answer:
                        answer = 'n'
                    if answer.lower() == 'y':
                        overwrite = True
                        break
                    elif answer.lower() == 'n':
                        break
                    else:
                        continue
                except EOFError:
                    break
            if not overwrite:
                non_written_files.append(cf)
                continue
        ifile = pkg_resources.resource_string(__name__,
                                              os.path.join('configs', cf))
        with io.open(ofile, 'w', newline=None) as of:
            of.write(_str_to_unicode(ifile))
        written_files.append(cf)

    if written_files:
        msg = "Copied files:  "
        msg = tw.TextWrapper(initial_indent=msg,
                             subsequent_indent=" "*len(msg))
        print(msg.fill(", ".join(written_files)))
        if non_written_files:
            msg = "Skipped files: "
            msg = tw.TextWrapper(initial_indent=msg,
                                 subsequent_indent=" "*len(msg))
            print(msg.fill(", ".join(non_written_files)))
    else:
        msg = "No file copied to directory '{}'"
        print(msg.format(args.to_dir))


if __name__ == '__main__':
    main()
