# Muon Beam Dump Calculations
### v0.0.2

Code by Rikab Gambhir
In collaboration with Samuel Alipour-fard and Cari Cesarotti 

# Usage

For the full analysis pipeline, see `muon_beam_dump.ipynb`. This pipeline automates cross section calculations, event yield calculations, and contour plot making.

For the specific models considered in `XXXX.XXXX`, and to reproduce all result plots within the paper, see the notebook `models.ipynb`.

# For Cari and Sam:

This is all my code for everything related to muon beam dumps

Take a look at `muon_beam_dump.ipynb` for the primary analysis pipeline. This does everything from calculating cross sections using (I)WW to calculating cross sections to calculating event yields to plotmaking. This is also where you get to specify all experimental parameters. Hopefully the user interface is straightforward, and everything should work out of the box. All data (cross sections, event yields, plots, and parameter files) is saved in the `experiments` directory. I have some data for a Water experiment already in there you can look at. 

The files `cross_sections.py` and `event_yields.py` are where most of the mathematics is actually done. You might find those interesting.

The rest of the files in this github are either old junk or just used to make auxillary plots, and can probably be ignored.



# Dependencies

(TODO: Write this section properly.)

Only standard python dependices (numpy, matplotlib, etc) are used. You may need to `pip install tqdm` which is a less-common python package for loading bars (I didnt get around to making this optional but nothing in the analysis depends on it). Finally, the plot colors will be uglier by default unless you have `rikabtools` installed, which is not possible at the moment.


# Changelog
v0.0.2: 13 June 2023. Exclusion plots.
v0.0.1: 30 May 2023. Uploaded to private repo.
