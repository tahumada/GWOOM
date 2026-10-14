# GWOOM

Set up the environment:

```bash
conda env create -f environment.yml
conda activate gwboom
```

Tile a skymap:

```bash
python tile_skymap.py ./skymaps/GRB261013a.fit -r 5 --level 0.9 --plot tiles.png
```

To run one cone search using filter `ztf_basic.json`:

```bash
python boom_filter_cone_search.py --filter ztf_basic.json --ra 135 --dec 14 --radius 5 --t0 "2026-10-14 00:00:00" --dt-days 1 --limit 20 --save-json test.json
```
