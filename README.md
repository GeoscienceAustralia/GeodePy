# ![GeodePy](https://github.com/GeoscienceAustralia/GeodePy/blob/master/docs/geodepy-logo.png)

GeodePy is a python package for precise geodetic and survey computations.

## Documentation

See [here](https://geodepy.readthedocs.io/) for documentation around downloading and using GeodePy.

## Features

GeodePy includes a variety of features for geodesy and geospatial data manipulation, including:

* Converting between coordinate types
* Transforming between datums
* Calculating geodetic distances and bearings
* Working with geoid models
* Surveying calculations
* Various classes for angles, coordinates, and datums
* Statistics
* And more!

## Installation

GeodePy is available on PyPi:

```console
$ pip install geodepy
```

## Dependencies

This package requires the following PyPI Packages installed:

```
NumPy
SciPy
```

Additionally, the geodepy.height module requires the GDAL library (tested using GDAL 3.0.4). For more information, see [here](https://gdal.org/index.html) for information about GDAL, [here](https://anaconda.org/conda-forge/gdal) for Anaconda support for GDAL and [here](http://www.gisinternals.com/release.php) for GDAL Binaries for Windows.

## Testing

Run: `python -m unittest discover geodepy/tests/ --verbose`

## Authors

* **Craig Harrison** - *Project Management* - [harry093](https://github.com/harry093)
* **Josh Batchelor** - *Initial Work, Geodesy and Surveying* - [BatchelorJ](https://github.com/BatchelorJ)
* **Jonathan Mettes** - *Testing, Integration and Deployment* - [jmettes](https://github.com/jmettes)
* **Jack McCubbine** - *Height Module* - [JackMcCubbineGA](https://github.com/JackMcCubbineGA)
* **Kyran Cook** - *Documentation and Uplift* - [Kyran-Cook](https://github.com/Kyran-Cook)

See also the list of [contributors](https://github.com/GeoscienceAustralia/geodepy/graphs/contributors) who participated in this project.

## License

Copyright 2018-2020 Geoscience Australia

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

