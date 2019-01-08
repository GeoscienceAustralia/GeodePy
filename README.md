# ![GeodePy](https://github.com/GeoscienceAustralia/GeodePy/blob/master/docs/geodepy-logo.png)

[![Travis](https://img.shields.io/travis/GeoscienceAustralia/GeodePy/master.svg?label=Travis%20CI)](https://travis-ci.org/GeoscienceAustralia/GeodePy) [![Coverage Status](https://coveralls.io/repos/github/GeoscienceAustralia/GeodePy/badge.svg)](https://coveralls.io/github/GeoscienceAustralia/GeodePy)

This is a package of tools for manipulating geospatial datasets using Python and tested in Python 3.6.4.

### Dependencies

This package requires the following non-standard Python modules installed:

```
NumPy
```

### Testing

Run: `python -m unittest discover geodepy/tests/`

### Tutorials

See [here](https://github.com/GeoscienceAustralia/GeodePy/tree/master/docs/tutorials) for worked examples of common GeodePy functions and routines.

## API

```
cd api/
virtualenv env
source env/bin/activate
pip install -r requirements.txt
zappa deploy dev
```

For subsequent updating run: `zappa update dev`

### Authors

* **Craig Harrison** - *Project Management* - [harry093](https://github.com/harry093)
* **Josh Batchelor** - *Initial Work, Geodesy and Surveying* - [BatchelorJ](https://github.com/BatchelorJ)
* **Jonathan Mettes** - *Testing, Integration and Deployment* - [jmettes](https://github.com/jmettes)

See also the list of [contributors](https://github.com/GeoscienceAustralia/geodepy/graphs/contributors) who participated in this project.

### License

Copyright 2018 Geoscience Australia

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

