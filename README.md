# Extended Metal Oxide Interface (xMOI)
Extended Metal Oxide Interface (xMOI) is a Python library designed to model, calculate and visualize the properties of metal (oxide) interfaces, particularly in the context of electrocatalysis. As our expectation, the library will provide tools for constructing surface/interface models, performing different calculations using various computational chemistry packages, and analyzing the results.

## Installation
You can install xMOI using pip:
```bash
git clone https://github.com/yourusername/xMOI.git
cd xMOI
pip install .
```


## Usage

```python
from xMOI import SurfaceModel

# Create a surface model
model = SurfaceModel(bulk_structure="path/to/bulk/structure")

# Generate surface models
model.generate_surface_models()

# Analyze the results
model.analyze_results()
```