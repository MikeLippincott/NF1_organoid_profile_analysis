# NF1_organoid_profile_analysis
This repo contains analysis code of profiles generated in multiple image-based profiling repos.

- 3D image-based profiles of NF1 organoids were generated from: [NF1 3D Organoid profiling pipeline](https://github.com/WayScience/NF1_3D_organoid_profiling_pipeline)
- 2D image-based profiles of NF1 organoids were generated from: [NF1 2D organoid profiling pipeline](https://github.com/WayScience/NF1_2D_organoid_profiling_pipeline)
- Cell types will be predicted using: [Cell type prediction in NF1 organoids](https://github.com/WayScience/NF1_organoid_cell_segmentation)

## Repo modules
- `0.download_data`: Download image-based profiles from the internet (not yet available)
- `1.EDA`: Exploratory data analysis of image-based profiles

## Computational environment
Notebooks in this repo are split between Python and R, each managed by a separate environment.

### Python (uv) and R (uvr) setup
Python notebooks use a `uv`-managed virtual environment defined in `pyproject.toml`/`uv.lock`.

```bash
source uv_setup.sh
```

This creates `.venv` and registers a `python3` Jupyter kernel.
Select this kernel when running the Python notebooks.

This also creates a `.uvr` directory for the R environment.

For more about uv and uvr, see:
- [uv documentation](https://github.com/astral-sh/uv)
- [uvr documentation](https://github.com/nbafrank/uvr)
