
# TED: Transient-Event Detector

Author: Joachim Mortensen  
License: [BSD 3-clause](http://opensource.org/licenses/BSD-3-Clause)

# About

This is the code that I wrote for master's thesis project.

# Requirements

TED runs in Python 2.7.6+, but not in Python 3+.

It depends on the packages in `requirements.txt` which can be installed using `pip` using e.g. `$ [sudo] pip install -r requirements.txt`. the `sudo` part is not needed if the installation is made in a virtual environment.

It also depends on some general-purpose modules that I have written. These are [`mplconf`](https://github.com/xidus/mplconf) and [`setup_and_parse`](https://github.com/xidus/setup_and_parse).

# Installation

Clone repository

    $git clone https://github.com/xidus/ted.git

or manually download the files in the master branch.

## Setup

Add package location to your `$PYTHONPATH` environment variable.

TED itself requires some paths to be set in a config file. 
To do so, copy the file `base.default.yaml` to `base.yaml` and edit the paths inside.

## Running the code

I have written a control script which makes it possible to use TED from the command line.

This script is available at [`ted-control`](#)


---




