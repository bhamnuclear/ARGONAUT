ARGONAUT
========

Argonaut is an R-matrix level profile calculator tool written in python and designed for use with [ReverseSisyphus][RS_repo]. It can be used independently for theoretical calculations instead of fitting to data as well, although the use changes slightly (see **example.py**).

How to use ARGONAUT
--------
There are a few ways to download/use ARGONAUT on your local python enviroment. I would recommend using an Anaconda environment as all required packages are usually installed.

The simplest method is to download the **RmTools.py** file from this repo, under the RmTools/ subfolder, and move it to the working folder location of your python code. 

Then when you call your other modules in your python file, you just need to include:

```python
from RmTools import Argonaut
```

See the **example.py** file for an example of how to use Argonaut for theoretical calculations without data.

To use this function with ReverseSisyphus you need to construct a function and then input that function as normal into ReverseSisyphus:

```python
def fitfunc(*fitting_parameters):
    return Argonaut(*fitting_parameters_and_fixed_parameters)

optimised_params = ReverseSisyphus(data,fitfunc,*args)
```

Typically the fitting parameters are defined as the centroid, width and scale in Argonaut.


A more generalisable way to download Argonaut is to clone the repo and use the following command in a terminal in the root directory:
```terminal
pip install .
```

from which then the package can be called as normal:
```python
from RmTools import Argonaut
```


Once installed, as before, call the function normally and utilise as above. 



Background
----------

Please see the following resources for more info on the R-matrix formalism:

'Nuclear Reactions for Astrophysics' by Thompson and Nunes


-------

[RS_repo]: https://github.com/bhamnuclear/ReverseSisyphus/tree/main