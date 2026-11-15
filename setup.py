# if setup tools is not installed, bootstrap it
try:
    import setuptools
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages

install_requires = ['numpy<2', 'scipy', 'astropy', 'astroquery', 'matplotlib',
                    'reproject', 'psycopg2-binary', 'six']

extras_require = {
    'oracle': ['cx_Oracle', ],
    'ds9': ['pyds9', ],
}

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

    python_requires=">=3.12",
    classifiers=["Development Status :: 3 - Alpha",
                 "Environment :: Console",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: BSD License",
                 "Operating System :: Unix",
                 "Programming Language :: Python :: 3",
                 "Programming Language :: Python :: 3.12",
                 "Topic :: Scientific/Engineering :: Astronomy",
                 ]
)
