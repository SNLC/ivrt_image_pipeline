# In Vitro Rabies Tracing Quantification Pipeline

Callaway lab in vitro rabies tracing quantification pipeline. 

Step 1: Colocalization test across channels. Groups XY coordinates within user-provided radius (default: 4 pixels) into cell objects with one-hot vector corresponding to features based on which channel(s) these coordinates originate in.  

Step 2: Grouping of starter cells into regions. Each starter cell (cell object positive for helper virus and rabies) is queried against every other starter cell in the image. If a starter cell falls within the region radius of another starter cell, those starter cells are grouped into a region.  

Step 3: Assignment of non-starter cells to regions. Each non-starter cell object is assigned a "Closest Region ID" corresponding to the "Region ID" of the nearest starter cell. If the distance between a non-starter cell object and the nearest starter cell is less than the user provided presynaptic radius (default: 300 pixels ~= 500 uM) than that non-starter cell object will be assigned a "Region ID." Cell objects with a distance greater than the user provided presynaptic radius will be assigned no region ID. Data output as 'region.csv'

Step 4: Statistics for the presynaptic target cell type (target of histology or in situ hybridization) are calculated for each region and output as a summary file ('summary.csv')

# Getting Started

### Prerequisites

```
Python 3 (only)
Conda or Pip
numpy
pandas
networkx
scikit-learn
```

### Installing
Installation inside virtual environment is recommended. See Python's virtualenv, pipenv, python
environment wrapper (PEW), or conda.

```
cd /path/to/github/repo/ivrt_image_pipeline
conda create --name ivrt-env --file conda_requirements.txt 
source activate ivrt-env

OR: 
pip install -r pip_requirements.txt
```

### Command line arguments

```
ivrt.py arguments:
  -h, --help            show this help message and exit
  --path                path to folder containing CSV [optional][default: current working directory]
  --alignment_radius    max XY distance in pixels within which XY point annotations across channels are 'aligned' to form same cell object [optional][default: 4 pixels]
  --region_radius       max XY distance in pixels within which multiple starter cells are grouped together into the same region (modified agglomerative clustering) [optional][default: 2 x presynaptic radius]
  --presynaptic_radius  max XY distance in pixels within which a presynaptic cell is assigned to an adjacent starter cell's region [optional][default: 300 pixels]
```

## Authors

* **Michael Lingelbach** - *Initial work* - [mjlbach](https://github.com/mjlbach)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
