HETDEX shot shuffling code by

Maximilian Fabricius (mxhf@mpe.mpg.de) Aug. 2011
Yi-Hao Chen (yihao@mpe.mpg.de) Nov. 2012
Francesco Montesano (montefra@mpe.mpg.de) 2014-2015
Jan Snigula (snigula@ex-mac05.mpe.mpg.de) 2015
Sungryong Hong (shong@astro.as.utexas.edu) 2015
Greg Zeimann (gregz@astro.as.utexas.edu) 2015
Karl Gebhardt (gebhardt@astro.as.utexas.edu) 2015
Daniel Farrow (dfarrow@mpe.mpg.de) 2015

Disclaimer
---------------------------------------
This code is very much in an alpha state. Feel free to change it, mess with it
or update it. I tried testing it and verifying its results.
However, there is no guarantees that the results are ok or that it want crash.

0) TO UPDATE EXISTING PACKAGE AT MOUNTAIN

- Login as hetdex (su - hetdex)
- go to /opt/het/hetdex/src/Virus/hetdexshuffle
- svn up
- pip install -U --no-deps --install-option="--prefix=/opt/het/hetdex" .

0a) TO COMMIT AN UPDATE AT THE MOUNTAIN
- Login as hetdex (su - htedex)
- go to /opt/het/hetdex/src/Virus/hetdexshuffle
- type "./touch_commit.sh"
- Write a message above the ---- line and save
- The change should be committed
- To commit just a single file: svn commit FILEPATH/FILENAME -m "WRITE MESSAGE HERE"

1) Installation
---------------------------------------

``hetdexshuffle`` can be installed using ``pip``, either from the checkout svn
repository or directly from the svn server.

For the first way, first get the source code with the following command:

    svn checkout svn://luna.mpe.mpg.de/hetdexshuffle/trunk hetdexshuffle

then cd into the newly created ``hetdexshuffle`` and run:

    pip install .

If you might need to use the
``--extra-index-url https://gate.mpe.mpg.de/pypi/simple/`` option to find some
of the dependences and the ``--user`` option if you are not running the above
command in a virtual/conda environment or as a sudo-er (be aware that running
``sudo pip`` is discouraged and my break some of your OS/distro packages if not
used with caution).
To update an existing shuffle installation, add the ``--upgrade`` options.

Alternatively you can install shuffle directly from the svn server with

    pip install svn+svn://luna.mpe.mpg.de/hetdexshuffle/trunk#egg=shuffle

If you are developing shuffle or don't want to reinstall ``hetdexshuffle`` every
time you update the svn repository that you have checked out, it is possible to
install it in "developer" mode:

    pip install -e .

This command should be rerun if the setup.py file changes.

2) Dependences
---------------------------------------

``hetdexshuffle`` depends on
    
    numpy       # numerics
    scipy
    astropy     # fits files
    matplotlib  # plotting
    pyds9       # interaction with ds9
    pyhetdex    # parse the fplane file
    astroquery  # Query Skyview
    six         # python 2/3 compatibility

If you want to be able to connect also with the curewise database, you need also
the following packages

    cx_Oracle   # connection with the database

``cx_Oracle`` can be installed asking pip to install the ``oracle`` optional
dependences:

    pip install .[oracle]

Note: to be able to install cx_Oracle you need Oracle itself. See e.g. here:
https://gist.github.com/thom-nic/6011715

2) Configuration
---------------------------------------

With the package we distribute a set of files, that can be retrieved with the
``shuffle_config`` command.

``shuffle.cfg`` is the configuration file that drives the execution of the
code. Some of the entries are the names of the output files or instructions
whether to create plots or open DS9 for visualisation, e.g.:

* [General]outfile: contains the new shot positions.
* [General]outinfo: contains extra information such as the coordinates of the
  actual selected guide and wavefront sensor stars.
  Those files will not be overwritten but appended to if they are already in place.
* [General]visualize: plot the shuffle result on SDSS or DSS image.

If you want to skip the shuffling and just search for guide and wavefront
sensor stars, set the ``radius`` command line argument to ``0``

Magnitude limits are set in section "MagLimits".

Note: Currently there are only magnitude limits, no colors.

``shuffle_config`` also copies the ``fplane.txt``. This contains the position
of the VIRUS IFUs and the mapping between the various ids for each one. It also
contain an entry for the ACAM (ifuslot id 000) and for the two LRS2 units (ifuslot
ids 056 and 066).

3) Running the code
---------------------------------------

The main executable installed by ``hetdexshuffle`` is ``do_shuffle``.

See ``do_shuffle -h`` for a list of options.

If the configuration or the fplane file are not found, ``do_shuffle`` will
warn the user and use the version shipped with the package.

Alternative ``shuffle.cfg`` can be provided with the ``-c/--config`` option.
Alternative fplane files can be provided to shuffle via the [General]fplane_file
configuration option.

For each shot, whose position, size are provided from the command line, shuffle
will retrieve guide stars, wavefront sensor stars and IFU stars candidates from
an online archive.

Currently is is set up to work with DR9. If outside the SDSS coverage, USNO is
ued. Note that USNO_A2 catalog has only magnitudes in 2 bands(B,R), thus the
criteria for finding stars are different. It is also possible to use a local
catalogue.

It is also possible to query the curewise database, if your system configure to
point to it and the oracle option is installed.

When running shuffle might try to query some online database, so it might fail
if the internet connection is absent.

With ``hetdexshuffle`` also comes the ``shuffle_distances`` executable.

ADD DOCUMENTATION HERE

4) Basic structure
---------------------------------------
# shuffle/shuffle.py
Main functions for shuffling.

# shuffle/findStars.py
Other functions to download the star catalog from SDSS DR7 (or USNO_A2) server,
to find star candidates, to translate position...etc.
The SQL query tool "sqlcl.py" is used to communicate with SDSS server.

# shuffle/visualize.py
Functions for visualization. It downloads the SDSS DR9 image and overplots the shuffle result.
If the coordinate is outside SDSS coverage, DSS image is then used.
It requires "matplotlib" library.

# shuffle/parang.py
Function for calculating the parallactic angle.

# shuffle/handle_local_catalogue.py
Handles reading in local catalogues of stars (instead of online databases)

# shuffle/do_shuffle_target.py
entry point for the do_shuffle executable

# shuffle/obj_dist_ifu_vis.py
entry point for the shuffle_distances executable

# shuffle/copy_config.py
entry point for the shuffle_config executable. The list of files to copy is
hard-coded in this file and must be updated if the number and/or name of
configuration files change

The main routine is "do" in "shuffle.py". It loops over the tracks for the given
target and shuffles them.

shuffle(...) is called by do(...) and moves the focal plane such to maximize the
number of stars that fall into the IFUs. It uses a very simple grid search
algorithm. Grid size and step width are set in shuffle.cfg.  shuffle(...) loops
over all grid points and finds the number of stars that do fall into IFUs. At
the end it calls pickShuffle(...). Currently pickShuffle sorts all tested grid
positions by number of IFU stars and shuffling distance and will pick a) and
largest number of stars and b) on distance. We may want to change this.

Finally pickShuffle(...) calls findGuideWFSSol(...) to find guide and wavefront
sensor stars. If it does not find a solution the shot will be considered as failed.

The code should be written following the indication from the official python
style guide (PEP8, https://www.python.org/dev/peps/pep-0008/)

5) A note on versioning and the svn repository
----------------------------------------------

The shuffle software is a moving target, so together with a version number
following the semantic versioning scheme (http://semver.org/) we use the svn
revision number.

To make sure that the svn revision number is correctly incorporated in the code
and accessible, we have added the ``$Revision$`` keyword to the
``shuffle/__init__.py`` file
(http://svnbook.red-bean.com/en/1.7/svn.advanced.props.special.keywords.html)
and, on the server side, we use a pre-commit hook to refuse commits that do not
update the aforementioned file. If a commit is refused a message similar to the
following is shown:

    Transmitting file data ..svn: E165001: Commit failed (details follow):
    svn: E165001: Commit blocked by pre-commit hook (exit code 1) with output:

    The file 'shuffle/__init__.py' must be modified.
    Run `touch_setting.sh` to touch the file and recommit.
    svn: E165001: Your commit message was left in a temporary file:
    svn: E165001:    'hetdexshuffle/svn-commit.tmp'

Running the ``touch_setting.sh`` executable will recommit the changes also
touching the file.
The executable can be used as a substitute for ``svn commit``: it touches the
file ``shuffle/__init__.py`` and run ``svn commit``. For more information run
``touch_setting.sh -h``

The revision number is available as in the example:

    import shuffle
    shuffle.__version__  # standard version
    shuffle.__svn_revision  # svn revision number
    shuffle.__full_version__  # combination of the above

6) Using a local catalogue
--------------------------

A local catalogue can be used as the source of stars for shuffle, this avoids the 
requirement of having a connection to an external database. For the finding charts
an external database is still required.

The location of the local catalogue is set by the localcat parameter in the config
file (set to None for an external catalogue). Currently only FITS catalogues, from SDSS
with the columns: 'objid', 'ra', 'dec', 'u', 'g', 'r', 'i', 'z' can be used. 

Note: To add non-SDSS catalogues src/handle_local_catalogue.py needs to be extended, and 
src/shuffle.py:do has to be modified to set the 'type' variable to something other than 
the default SDSS. 
