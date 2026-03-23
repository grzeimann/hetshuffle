# if setup tools is not installed, bootstrap it
try:
    import setuptools
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages

install_requires = ['numpy', 'scipy', 'astropy', 'astroquery', 'matplotlib',
                    'pyds9', 'pyhetdex', 'six']

extras_require = {'oracle': ['cx_Oracle', ]}

entry_points = {'console_scripts': ['do_shuffle = '
                                    'hetdex_shuffle.do_shuffle_target:main',
                                    'shuffle_config = '
                                    'hetdex_shuffle.copy_config:main',
                                    'shuffle_distances ='
                                    'hetdex_shuffle.obj_dist_ifu_vis:main',
                                    ]}

setup(
    # package description and version
    name="hetdex-shuffle",
    version='0.3.5',
    author="HETDEX collaboration",
    author_email="",
    description="The shuffle code",
    long_description=open("README.txt").read(),

    # list of packages and data
    packages=find_packages(),
    # get from the MANIFEST.in file which extra file to include
    include_package_data=True,
    # don't zip when installing
    zip_safe=False,

    # entry points: creates vhc script upon installation
    entry_points=entry_points,
    # dependences
    install_requires=install_requires,
    extras_require=extras_require,

    classifiers=["Development Status :: 3 - Alpha",
                 "Environment :: Console",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License (GPL)",
                 "Operating System :: Unix",
                 "Programming Language :: Python :: 2.7",
                 # "Programming Language :: Python :: 3.4",
                 # "Programming Language :: Python :: 3.5",
                 "Topic :: Scientific/Engineering :: Astronomy",
                 ]
)
