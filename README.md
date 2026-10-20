# GWOOM

Set up the environment:

For conda
```bash
conda env create -f environment.yml
conda activate gwboom
```

For python venv
```bash
pip install -e .
```


Tile a skymap:

```bash
python gwoom/tile_skymap.py path/to/skymap -r 5 --level 0.9 --plot
```

Tile a skymap and search every tile in ZTF + LSST with the filters in
`filters/filters.py`, keeping objects first detected no earlier than 1 h before
the trigger (read from the skymap header, or `--t0`). ZTF/LSST matches are
merged and written to a text file:

```bash
python gwoom/multisurvey_query.py ./skymaps/GRB261013a.fits -r 5 --level 0.9 \
    --max-pre-hours 1 --dt-days 2 -o transients.txt --tiles-out tiles.csv
```

To-Do

1. Generalize the cone search
2. Display the candidates in a PDF
