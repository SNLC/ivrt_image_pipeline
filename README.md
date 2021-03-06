# In Vitro Rabies Tracing Quantification Pipeline

Callaway lab in vitro rabies tracing quantification pipeline. Takes CSV files corresponding to XY point annotations for multichannel imaging experiments. 

Step 1: Colocalization test across channels. Groups XY coordinates within user-provided radius (default: 4 pixels) and assigns each XY coordinate features based on which channel these coordinates originate in.  

Step 2: Colocalization test across starter cells. Each starter cell is queried against every other starter cell on the dish. If the starter cell falls within a regional radius, it is added to the cluster.  

Step 3: Grouping of non-starter cells by region. Each starter cell is assigned a Closest Region ID based on the nearest adjacent starter cell. Each starter cell is given  


# Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Python 3
Conda or PIP
numpy
pandas
networkx
scikit-learn
```

### Installing

```
conda install numpy pandas scikit-learn networkx
```

### Command line arguments

```
optional arguments:
  -h, --help            show this help message and exit
  --path                path to folder containing CSV
  --alignment_radius    max xy distance in pixels for two points to be 'aligned' to same cell
  --presynaptic_radius  Max XY distance which a presynaptic cell be from a starter cell and  
                        still be classified in that starter cell's region
  --region_radius       Max XY distance in which two starter cells are grouped together into  
                        the same region
```

## Authors

* **Michael Lingelbach** - *Initial work* - [mjlbach](https://github.com/mjlbach)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


