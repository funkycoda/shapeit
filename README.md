# shapeit

Creates shapefiles for either points or polygons with a 20KM buffer

## How to use:

### Using virtualenv (recommended)
Best way to use this is to use [virtualenv](http://virtualenv.readthedocs.io/en/stable/installation/)

Clone the repo
```
git clone https://github.com/funkycoda/shapeit.git
```

Navigate to the folder and enable it for ```virtualenv``` 
```
cd shapeit
virtualenv .
source bin/activate
```

Install the necessary dependencies
```
pip install -r requirements.txt
```

Generate polygons
```
python shape_it.py -o output examplePoints.csv
python shape_it.py -g MasterTax_SciName -lat lat_dec -lon lon_dec -o testing exampleForRanges.csv
```

When done, deactivate from ```virtualenv``` 
```
deactivate
```

### Without using virtualenv
Clone the repo
```
git clone https://github.com/funkycoda/shapeit.git
```

Navigate to the folder
```
cd shapeit
```

Install the necessary dependencies
```
pip install -r requirements.txt
```

Generate polygons
```
python shape_it.py -o output examplePoints.csv
python shape_it.py -g MasterTax_SciName -lat lat_dec -lon lon_dec -o testing exampleForRanges.csv
```




