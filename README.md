# HET Shuffle (hetdex-shuffle)

The HET shuffle tool selects viable guide and wavefront-sensor stars and produces verified target placements for HET observations. This README is written for on-mountain usage by astronomers/TOs setting up HET tracks.

If you are looking for developer-oriented details, see README.txt and the configs under hetdex_shuffle/configs/.


- Console entry points installed:
  - do_shuffle – main command to compute a shuffle for a target
  - shuffle_config – copies example configuration files to your working directory
  - shuffle_distances – helper for IFU/probe distance visualization


Quick start
- Ensure you have the required configuration files (shuffle.cfg and a valid fplane file).
- Install the package (see Installation below).
- Run do_shuffle with your target coordinates, track, IFU slot, optional offsets, and desired options.


Installation
You can install either with pip (PEP 517) or create a dedicated Conda environment.

Option A: pip (system/virtualenv)
- Python 3.8+
- From a checkout of this repository:
  - pip install .
  - Optional Oracle DB support: pip install .[oracle]
- Notes about dependencies
  - pyhetdex is listed as a normal dependency. If your site requires a custom index, install it with:
    - pip install --extra-index-url https://gate.mpe.mpg.de/pypi pyhetdex
  - pyds9 and ds9 require X11 if you intend to use DS9 visualization.

Option B: Conda environment
- conda env create -f environment.yml
- conda activate hetshuffle
- If pyhetdex is not found, install via pip with the extra index as needed:
  - pip install --extra-index-url https://gate.mpe.mpg.de/pypi pyhetdex


Configuration files you will need
- shuffle.cfg: Main runtime configuration. A ready-to-edit default lives at hetdex_shuffle/configs/shuffle.cfg. You can copy it to your working directory with:
  - shuffle_config
- fplane file: e.g., fplane.txt or fplaneall.txt from your site’s distribution. Path is passed at runtime via --fplane_file.


Command-line usage
Basic syntax (positional arguments then options):
- do_shuffle RA DEC RADIUS TRACK IFUSLOT [X_OFFSET] [Y_OFFSET] [TARGET] [options]
  - RA: Right Ascension of the target, in hours (hh:mm:ss.ss)
  - DEC: Declination of the target, in degrees (±dd:mm:ss.ss)
  - RADIUS: Search radius in degrees for suitable guide/WFS stars (typical 0.0–1.0). 0.0 means use exact target area and rely on probes.
  - TRACK: 0 = East, 1 = West
  - IFUSLOT: Three-digit IFU slot ID (e.g., 058). Boresight is 000.
  - X_OFFSET, Y_OFFSET: Optional small ACAM-plane offsets from the IFU center, in arcseconds (default 0.0, 0.0)
  - TARGET: Optional name for the observation (default "Test")

Common options
- --config PATH_TO/shuffle.cfg
- --fplane_file PATH_TO/fplaneall.txt
- --catalog NAME (e.g., GAIA, PANSTARRS)
- -v or --verbose (can be repeated)
- --visualize PATH_TO/output.png to save sky overlays
- --visualize_ACAM PATH_TO/acam.png to save ACAM-frame overlay
- --visualize_probestars PATH_TO/probes.png for probe-star view
- --epoch, --pmRA, --pmDEC for proper-motion handling

Example (generalized)
- do_shuffle 08:00:16.08 "+40:29:56.00" 0 0 058 25.0 25.0 CARLA_J0800+4029_3 \
  --catalog PANSTARRS -v \
  --fplane_file ~/work/Backup/cure/Remedy/fplaneall.txt \
  --config ~/work/code/hetshuffle/hetdex_shuffle/configs/shuffle.cfg

Notes on the example
- RA and DEC are quoted if they contain + or spaces.
- RADIUS=0, TRACK=0 (East), IFUSLOT=058, offsets of +25.0" in X and Y, target name CARLA_J0800+4029_3.
- Using the PANSTARRS catalog for star selection.
- Using a specific fplane file and a known shuffle.cfg.


Expected output
What you should expect when a shuffle completes successfully:
- Terminal log
  - A summary of parsed inputs, selected guide and WFS stars, and computed probe angles.
  - Warnings if targets fall outside patrol regions or if magnitude/sep cuts are tight.
- Output files as directed by shuffle.cfg [General] section
  - outfile: A text file listing the new shot positions; entries are appended if the file already exists.
  - outinfo: A companion text file that includes extra information such as the coordinates and IDs of the selected guide and wavefront-sensor stars.
  - Optional visualization products if enabled or passed via CLI:
    - --visualize: Sky view PNG (DSS/SDSS/HiPS background) with IFU/probe overlays.
    - --visualize_ACAM: ACAM-frame PNG. The tool also writes a CSV next to this image with columns:
      RA, Dec, acam_x, acam_y, img_x, img_y, g, r, i
      Only stars that fall within the ACAM footprint are listed.
    - --visualize_probestars: PNG highlighting probe-star geometry and selections.
- DS9
  - If configured, a DS9 window may be launched and annotated according to your shuffle.cfg settings.

Return codes
- 0 on success. Non-zero if configuration or geometric constraints fail, or if catalogs are unavailable.

Tips and troubleshooting
- If catalogs time out, retry or switch to a different catalog via --catalog.
- If no viable guide/WFS stars are found, try increasing RADIUS slightly or adjusting magnitude limits in shuffle.cfg.
- Ensure the IFUSLOT is three digits (e.g., 058). The code normalizes it internally, but it’s good practice to supply the padded form.
- For proper motion targets, supply --epoch, --pmRA, and --pmDEC to propagate to the observation epoch.
- For ACAM orientation/offsets, relevant parameters are in shuffle.cfg under [General] and [offsets].

Where things live
- Example configs: hetdex_shuffle/configs/
- Main executable module: hetdex_shuffle/do_shuffle_target.py
- Visualization code: hetdex_shuffle/visualize.py

Acknowledgements and license
- HETDEX collaboration; see LICENSE for terms.
- Portions rely on Astropy, Astroquery, Matplotlib, NumPy, and pyhetdex.