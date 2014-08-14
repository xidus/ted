
# Introduction

This set of procedures explains how I set up the work environment on the remote-computing facility that I had the opportunity to use for my master's-thesis work.

It basically contains all successful commands and settings changes that I wrote to get all the necessary tools working.

---

In short, I set up the command-line tool `tsocks` on the `imagediskserver`s in order to communicate with the internet through the gateway server that shields the facility.

My project uses Python, and I explain how I got a Python distribution working on the `image[disk]server`s.

As part of the image-processing procedures, I installed two special tools for handling astrometrical registration of cutouts and sophisticated difference imaging.

In addition to this, I have installed the version-control system Git in order to ease the workflow. I have used this to synchronise code between my own computer and the remote facility which runs the code with the data.

---

To begin with, I created the following directories in my user's home directory `$HOME`

    $ mkdir -p $HOME/downloads/software
    $ mkdir $HOME/.local

All software is downloaded `$HOME/downloads/software`, and the compiled binaries, libraries and header files are put into sub-directories `$HOME/.local`.

For all the special settings that are specific to my project, I have set up a file that supplements the `$HOME/.bashrc` shell-login script. This supplement file is called `$HOME/.bashrc_skyml` and is sourced from `$HOME/.bashrc`.

To begin with, this additional login-script contains the following:

    # For installation purposes
    export PREFIX=$HOME/.local
    export PATH=$PREFIX/bin:$PATH

More lines are added with the procedures below.

---

__Content:__

* [Setting up *tsocks* on the server](#tsocks)
    * [Setting up SSH configuration file](#sshconfig)
* [Setting up Python](#anaconda)
* [Installing HOTPANTS](#HOTPANTS)
* [Installing WCSREMAP](#WCSREMAP)
* [Installing Git on the server](#git)
* [Configure IPython notebook on the server](#ipynb)

---

<a id="tsocks"></a>

## Setting up *tsocks* on the server

This procedure makes it possible to wrap a command-line program with the `tsocks` utility in order to make sockets go through a dynamic-forwarding port that is connected by SSH to a chosen gateway server. It is available on `imagediskserver[1-3]`, but not on `imageserver[1-3]` (IBM Blade servers).

### Installing *tsocks* on *imageserver[1-3]*

*Link on the `tsocks` project site does not work. __Use the one on sourceforge instead.__*

    $ cd $HOME/downloads/software
    $ tsocks wget http://downloads.sourceforge.net/project/tsocks/tsocks/1.8%20beta%205/tsocks-1.8beta5.tar.gz?r=&ts=1392757146&use_mirror=optimate
    $ tar -xzvf tsocks-1.8beta5.tar.gz
    $ cd tsocks-1.8
    $ ./configure --prefix=$PREFIX
    $ make

Now running

    $ make install

gives

    /bin/sh mkinstalldirs  "/home/zpq522/.local/bin"
    /usr/bin/install -c tsocks /home/zpq522/.local/bin
    /bin/sh mkinstalldirs  "/lib"
    /usr/bin/install -c libtsocks.so.1.8 /lib
    /usr/bin/install: cannot create regular file `/lib/libtsocks.so.1.8': Permission denied
    make: *** [installlib] Error 1

which I know no way out off.

*I added the installation directory to `$PATH` instead.*

    PATH=$HOME/downloads/software/tsocks-1.8:$PATH

*__Unfortunately, when using it in my batch script for Office Grid, I get a linking error and nothing can connect to the world wide web anyway :( :( :( !__*

    # Excerpt from my run.sh script
    tsocks /abs/path/to/main.py

Output similar to this:

    ERROR: ld.so: object '/usr/lib/libtsocks.so' from LD_PRELOAD cannot be preloaded: ignored.

I have tried setting the LD_PRELOAD environment variable before the statement like so

    LD_PRELOAD=/home/zpq522/downloads/software/tsocks-1.8:$LD_PRELOAD tsocks wget -O- checkip.dyndns.org

with the following (error) output:

    ERROR: ld.so: object '/home/zpq522/downloads/software/tsocks-1.8' from LD_PRELOAD cannot be preloaded: ignored.
    ERROR: ld.so: object '/home/zpq522/downloads/software/tsocks-1.8' from LD_PRELOAD cannot be preloaded: ignored.
    ERROR: ld.so: object '/usr/lib/libtsocks.so' from LD_PRELOAD cannot be preloaded: ignored.
    ERROR: ld.so: object '/home/zpq522/downloads/software/tsocks-1.8' from LD_PRELOAD cannot be preloaded: ignored.
    --2014-02-26 21:57:06--  http://checkip.dyndns.org/
    Resolving checkip.dyndns.org... 91.198.22.70, 216.146.38.70, 216.146.39.70, ...
    Connecting to checkip.dyndns.org|91.198.22.70|:80... failed: Connection timed out.

This would explain the hang-up errors that I get when I ran my batch script.

"Great".


---


### Using *tsocks*

To utilise `tsocks`, requires a simple configuration file and access to a server that can act as a proxy that forwards connections.

I created the file `$HOME/tsocks.conf` containing

    server = 127.0.0.1
    server_port = 1080

and added the following environment variable to my login script:

    export TSOCKS_CONF_FILE=$HOME/tsocks.conf

In my case, the chosen gateway server is `ask.diku.dk`. I need to make it available for socket connections made by `tsocks`. I do so by establishing a connectiion to the gateway server which allows for dynamic forwarding. This can be done in the following way:

    $ ssh -D 1080 user@server

With this connection open, everything that communicates with port 1080 on localhost is forwarded to the gateway server which decides how traffic moves on from there. In another terminal I can then use `tsocks` to send connections to this port, e.g. like so:

    $ tsocks ipython
    In [1]: import requests
    In [2]: print requests.get('http://checkip.dyndns.org/').content
    <html><head><title>Current IP Check</title></head><body>Current IP Address: 130.225.96.225</body></html>

which wraps an IPython shell (could also be ordinary Python shell), or

    $ tsocks wget http://sdssdp62.fnal.gov/sdsssn/gallery/sn_gallery.200567.x2.oname.jpg

which downloads a gallery of confirmed SNe Ia.

*When logged on to the server, I run GNU screen to have more than one terminal running at the same time.*

---

<a id="sshconfig"></a>

### Using the SSH configuration file

To ease things, I created an asymmetric key-pair on the `imagediskserver`

    $ mkdir ~/.ssh
    $ cd !$
    $ ssh-keygen # create key-pair named ask-ids(.pub)
    ... #
    $ scp ask-ids.pub ask.diku.dk:~/.ssh/.

On `ask.diku.dk`, I added `ask-ids.pub` to `authorized_keys` like so

    $ # On ask.diku.dk
    $ cd ~/.ssh
    $ cat ask-ids.pub >> authorized_keys

Back on `imagediskserver`, I set up the connection in `~/.ssh/config` in the following way

    Host ask
        User zpq522
        HostName ask.diku.dk
        IdentityFile /home/zpq522/.ssh/ask-ids
        DynamicForward 1080

This allows me to connect to `ask.diku.dk` with dynamic forwarding by simply typing

    $ ssh ask

---

__References:__

*   [TSOCKS](http://tsocks.sourceforge.net/about.php)
*   [Web scraping: Reliably and efficiently pull data from pages that don't expect it](https://www.youtube.com/watch?v=52wxGESwQSA) (video, royughly 1h 25m in)
*   `$ man ssh_config`

---

<a id="anaconda"></a>

## Setting up Python

### Software

__Required to obtain a working Python distribution__

* [Anaconda](http://docs.continuum.io/anaconda/index.html)

__Optional software__

* [PIP](http://www.pip-installer.org/en/latest/) : The version that comes with Anaconda does not support SSL
* [virtualenv](http://www.virtualenv.org/en/latest/) : Create isolated Python environments where installed packages do not interfere with the system-wide installation. 
* [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/) : Wraps `virtualenv` to ease using and switching between different virtualenvs.

---


### Required steps to obtain a working Python distribution on the server

The link to the latest installer can be obtained from [this page with all available downloads from Continuum Analytics](http://continuum.io/downloads).

Download and install with the following lines:

    $ cd $HOME/downloads/software
    $ tsocks wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-1.6.1-Linux-x86_64.sh
    $ bash Anaconda-1.6.1-Linux-x86_64.sh
    ...
    >>> yes                 # Agreeing to licence terms
    >>> $PREFIX/anaconda    # Choose a custom install location.,

To make this Python distribution the primary one, pre-pend the path to the binaries to `$PATH` in `~/.bashrc_skyml`

    export PATH=$PREFIX/anaconda/bin:$PATH

__Note: To use the right Python distribution when running the executable Python files I replace the shebang `#!/usr/bin/python` with `#!/usr/bin/env python`, since `env` uses the first `python` command it finds in `$PATH`.__


---

#### Updating existing packages using *conda*

PIP does not work with `tsocks`, and I have not yet found a way to obtain or update packages through this package installer.

Luckily, most of what I need for my project is already included in the Anaconda distribution. It also comes with its own package manager called `conda` and which *does* work with `tsocks`. To update a series of packages, I can simply write something like the following:

    $ tsocks conda update ipython ipython-notebook matplotlib scipy numpy pandas sqlite sqlalchemy yaml astropy

---


#### Installing packages manually

Since not all the packages that I work with are included in the [packages included in Anaconda][anaconda-pkg], and since PIP does not seem viable from behind the firewall, I have had to install these additional packages manually in either of the following ways.

The basic steps are to

1.  Download the tar-ball for the package to the server using `tsocks` and `wget` as described above, perhaps needing the flag `--no-check-certificate` added to the `wget`command[^WGET-ENCRYPT].
2.  Untar, and from the package directory run `$ python setup.py install` (no `sudo` needed, since I am using Anaconda locally).

[anaconda-pkg]: http://docs.continuum.io/anaconda/pkg-docs.html

[^WGET-ENCRYPT]: This is not always needed, e.g. like when I get the tar-ball from the AstroML package at github.com, even though it downloaded using an TLS/SSL-encrypted connection (HTTPS). I do not know why.

##### Example: Install package *mechanize* by downloading it from PyPI

From <https://pypi.python.org/pypi/mechanize> I obtain the path to the latest version of `mechanize`.

On the server, I then type

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/m/mechanize/mechanize-0.2.5.tar.gz#md5=32657f139fc2fb75bcf193b63b8c60b2 --no-check-certificate
    $ tar -xzf mechanize-0.2.5.tar.gz
    $ cd mechanize-0.2.5
    $ python setup.py install

---


##### Example: AstroML (without Git)

Have a look at [the AstroML master branch](https://github.com/astroML/astroML).

    $ cd $HOME/downloads/software
    $ tsocks wget https://github.com/astroML/astroML/archive/master.tar.gz -O astroml.tar.gz
    $ tar -xzf astroml.tar.gz
    $ cd astroML-master
    $ python setup.py install

---


##### Example: AstroML (with Git)

    $ cd $HOME/repositories
    $ tsocks git clone git://github.com/astroML/astroML.git
    $ cd astroML
    $ python setup.py install

---



### Optional steps in order to update *PIP* and install *virtualenv* and *virtualenvwrapper* on the server

__For now, I am using the 'system-wide' local install as my basis for running the programs that I write in this project. However, I did perform the following steps in case I need to test something that may cause the rest of my code to not work.__


#### PIP

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/p/pip/pip-1.4.1.tar.gz --no-check-certificate
    $ tar -xzf pip-1.4.1.tar.gz
    $ cd pip-1.4.1
    $ python setup.py install

#### Virtualenv

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz --no-check-certificate
    $ tar -xzf virtualenv-1.10.1.tar.gz
    $ cd virtualenv-1.10.1
    $ python setup.py install

#### Virtualenvwrapper (and dependencies)

Virtualenvwrapper depend on the following packages which did not come with Anaconda.

* [pbr](https://pypi.python.org/pypi/pbr) (Python Build Reasonableness)
* [virtualenv-clone](https://pypi.python.org/pypi/virtualenv-clone) (Script to clone virtualenvs)
* [stevedore](https://pypi.python.org/pypi/stevedore) (Manage dynamic plugins for Python applications)

The procedure was the same as above

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/p/pbr/pbr-0.5.21.tar.gz --no-check-certificate
    $ tar -xzf pbr-0.5.21.tar.gz
    $ cd pbr-0.5.21.tar.gz
    $ python setup.py install

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/v/virtualenv-clone/virtualenv-clone-0.2.4.tar.gz --no-check-certificate
    $ tar -xzf virtualenv-clone-0.2.4.tar.gz
    $ cd virtualenv-clone-0.2.4
    $ python setup.py install

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/s/stevedore/stevedore-0.11.tar.gz --no-check-certificate
    $ tar -xzf stevedore-0.11.tar.gz
    $ cd stevedore-0.11
    $ python setup.py install

Finally, `virtualenvwrapper` itself is installed in the same way

    $ cd $HOME/downloads/software
    $ tsocks wget https://pypi.python.org/packages/source/v/virtualenvwrapper/virtualenvwrapper-4.1.1.tar.gz --no-check-certificate
    $ tar -xzf virtualenvwrapper-4.1.1.tar.gz
    $ cd virtualenvwrapper-4.1.1
    $ python setup.py install

To make `virtualenvwrapper` work, set the directory that is to contain all the virtualenvs and make the `virtualenvwrapper` scripts available by adding the following to `~/.bashrc_skyml`:

    # Virtualenvwrapper
    export WORKON_HOME=$HOME/.virtualenvs
    source $PREFIX/anaconda/bin/virtualenvwrapper.sh

---


### Local setup on my own machine

I replicated the Anaconda setup and package update except from installing virtualenv and virtualenvwrapper which I already have.

I am using a virtualenv to easily switch between the system-wide python environment and a my work tools, but instead of using the Python executable that I get with the creation of this virtualenv, I create hooks, so that when I activate the environment, the Anaconda distribution is placed first in `$PATH`.


---



<a id="HOTPANTS"></a>

## Installing HOTPANTS on the server

The High Order Transform of PSF ANd Template Subtraction (HOTPANTS).

__References:__

*   [HOTPANTS](http://www.astro.washington.edu/users/becker/v2.0/hotpants.html)
*   [FITSIO home page](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html)

### Dependency: CFITSIO *.h and libcfitsio

__Steps:__

    $ cd $HOME/downloads/software
    $ tsocks wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3350.tar.gz
    $ tar -xzf cfitsio3350.tar.gz
    $ cd cfitsio
    $ ./configure --prefix=$PREFIX
    $ make
    $ make install

---


### HOTPANTS

__Steps:__

    $ cd $HOME/downloads/software
    $ tsocks wget http://www.astro.washington.edu/users/becker/v2.0/software/hotpants_v5.1.10.tar.gz
    $ tar -xzf hotpants_v5.1.10.tar.gz
    $ cd hotpants_v5.1.10b
    $ nano Makefile

In `Makefile` I made the following changes:

    # CFITSIOINCDIR =  ../../include/cfitsio
    # LIBDIR        =  ../../lib/$(ARCH)
    CFITSIOINCDIR = $(PREFIX)/include
    LIBDIR = $(PREFIX)/lib

Exiting and then running

    $ make

I got the following error

    gcc  -funroll-loops -O3 -ansi -pedantic-errors -Wall -I/home/zpq522/.local/include  -c main.c
    In file included from main.c:6:0:
    /home/zpq522/.local/include/fitsio.h:120:18: error: ISO C90 does not support ‘long long’ [-Wlong-long]
    main.c: In function ‘main’:
    main.c:250:4: warning: implicit declaration of function ‘gethostname’ [-Wimplicit-function-declaration]
    make: *** [main.o] Error 1

which is where I left it for now (2013-11-10).

__Update:__

Kim found that

> It was some problems with the gcc parameters in the Makefile - this was done in a non-portable fashion.

The following additional changes to `Makefile` made it compile without error (but not without warnings)

    # COPTS = -funroll-loops -O3 -ansi -pedantic-errors -Wall -I$(CFITSIOINCDIR)
    # LIBS  = -L$(LIBDIR) -lm -lcfitsio
    COPTS = -funroll-loops -O3 -ansi -std=c99 -Wall -I$(CFITSIOINCDIR)
    LIBS  = -L$(LIBDIR) -lcfitsio -lm

There is no `make install`directive, so the compilation output has to be moved manually to the desired directory for executables:

    cp extractkern hotpants maskim $PREFIX/bin/.



---



<a id="WCSREMAP"></a>

## Installing WCSREMAP on the server

> The first step in the process of difference imaging is to astrometrically align two input images. This should ensure that a given astronomical object is centered at the same pixel in both images. The WCS-format astrometric distortions in the input images should model the distortions in the focal plane to as high order as possible. Any inaccuracies in this distortion model will result in mis-aligned images, which in turn will effect the quality of their difference image.
>
> *-- Andrew Becker, author of WCSREMAP*

__References:__

*   [WCSREMAP](http://www.astro.washington.edu/users/becker/v2.0/wcsremap.html)
*   [FITSIO home page](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html)
*   [WCSTools](http://tdc-www.harvard.edu/wcstools/)

### Dependency: CFITSIO *.h and libcfitsio

*Same dependency as for HOTPANTS. [See installation steps.](#HOTPANTS)*

---


### Dependency: WCSTools

__Steps:__

    $ cd $HOME/downloads/software
    $ tsocks wget http://tdc-www.harvard.edu/software/wcstools/wcstools-3.8.7.tar.gz
    $ tar -xzf wcstools-3.8.7.tar.gz
    $ cd wcstools-3.8.7
    $ make

---


### WCSREMAP

__Steps:__

    $ cd $HOME/downloads/software
    $ tsocks wget http://www.astro.washington.edu/users/becker/v2.0/software/wcsremap-1.0.1.tar.gz

Add includes and libraries in `Makefile`:

    # CFLAGS= -DLOWENDIAN -Wall -funroll-loops -O3
    # CLIBS= -lcfitsio -lm -lwcs
    CFLAGS= -I$(PREFIX)/include -DLOWENDIAN -Wall -funroll-loops -O3
    CLIBS= -L$(PREFIX)/lib -L$(HOME)/downloads/software/wcstools-3.8.7/libwcs -lcfitsio -lm -lwcs

Then compile and copy the binary to `$PREFIX/bin`

    $ make
    gcc  -I/home/zpq522/.local/include -DLOWENDIAN -Wall -funroll-loops -O3   -c -o wcsremap.o wcsremap.c
    wcsremap.c:8:16: fatal error: wcs.h: No such file or directory
    compilation terminated.
    make: *** [wcsremap.o] Error 1

For some reason, the libraries are not included on the path... Left it here. __I have it installed on my own machine.__

<del>

    $ cp wcsremap $PREFIX/bin/.

</del>

---



<a id="git"></a>

## Installing Git on the server

Dependencies not already on `imagediskserver[1-3]`[^D]:

* [zlib](http://zlib.net/)
* [GNU's gettext](https://www.gnu.org/software/gettext/)

<!--* either [tcl/tk](http://www.tcl.tk/software/tcltk/download.html) (two distinct packages)-->

<!-- Therefore, whenever any software that I am installing is not finding the necessary files, I just install them manually (locally). -->

[^D]: Found out by letting the compiler tell what is missing and then install each library as needed.

### Dependency: zlib

    $ cd $HOME/downloads/software
    $ tsocks wget http://zlib.net/zlib-1.2.8.tar.gz
    $ tar -xzf zlib-1.2.8.tar.gz
    $ cd zlib-1.2.8
    $ ./configure --prefix=$PREFIX
    $ make
    $ make install

---



### Dependency: gettext

`gettext` *is* installed on `imagediskserver[1-3]`, but it was not used when I ran `make` in the Git directory.

    $ cd $HOME/downloads/software
    $ tsocks wget http://ftp.gnu.org/pub/gnu/gettext/gettext-0.18.3.1.tar.gz
    $ tar -xzf gettext-0.18.3.1.tar.gz
    $ ./configure --prefix=$PREFIX
    $ make
    $ make install

---


### Git

I found out that the following wil suffice in order to make Git compile:

    $ cd $HOME/downloads/software
    $ tsocks wget http://git-core.googlecode.com/files/git-1.8.4.tar.gz
    $ tar -xzf git-1.8.4.tar.gz
    $ cd git-1.8.4/
    $ # using `--with-tcltk`, it wlil only compile with this graphical interface, if it can be found, otherwise it is left out.
    $ ./configure --prefix=$PREFIX --with-zlib=$PREFIX --with-tcltk
    $ make
    $ make install


#### Git config file

I use the following in my `~/.gitconfig` file on the server:

    [user]
    name = [my name] @ SkyML
    email = [email not shown]

    [core]
        editor = nano

    [color]  # Ensure we get color
        diff = auto
        branch = auto
        status = auto

    [http]
        sslVerify = false

    [alias]
        ci = commit
        st = status
        co = checkout
        oneline = log --pretty=oneline
        br = branch
        la = log --pretty=\"format:%ad %h (%an): %s\" --date=short



#### Setting up the remote repository on the server

I create a Git repository for different parts of my project. I am building a Python package of utilities that I call `ted` and which I use in my separate and project-specific codebase which does not have a name. The project-specific codebase contains a mix of many different things (or will do so, when I have finished my project).

I set up a bare repository on another server (at NBI) that I can push to from both the image servers and my own computer.

I show it only for one of my repositories.

On my bare-repository server:

    $ mkdir -p ~/git/ted.git
    $ cd !$
    $ git init --bare

---

On my own machine, I add this repository as my remote origin.

    $ git remote add origin ssh://fys/~/git/ted.git
    $ git push origin master

The last line pushes my master branch to the remote repository, so that is now updated with everything that I had committed so far on my own computer. `fys` is set up as a host in my SSH-config file.

As Git will likely tell you, thelast command above can be reduced to just `$ git push|pull` if something like the following is added to one's local repository configuration. In `.git/config`, I have

    ...
    [remote "origin"]
        url = ssh://fys/~/git/ted.git
        fetch = +refs/heads/*:refs/remotes/origin/*
    [branch "master"]
        remote = origin
        merge = refs/heads/master


---

On the image servers, I once again had to find a slightly more complicated solution. Here is what I did:

    $ cd $HOME/git
    $ tsocks git clone ssh://fys/~/git/ted.git

And to push changes made while on the server, I just use

    $ tsocks git push

---

__Notes:__

*Removing a previously added file `<fname>` from the repository*

    $ git rm --cached <fname>

*Undoing everything up to the last commit*

    $ git reset --hard

*Unstaging (not undoing changes)*

    $ git reset

---

__References:__

*Setting up Git repositories*

* <https://www.digitalocean.com/community/articles/how-to-install-git-on-ubuntu-12-04>
* [this guide](http://tumblr.intranation.com/post/766290565/how-set-up-your-own-private-git-server-linux "Under 'Add your repositories'")

*Clone from behind a firewall*

* <https://julianscorner.com/wiki/programming/git/githubbehindfirewall>
* <http://knowledgefrontier.blogspot.dk/2010/03/how-to-clone-git-repo-behind-firewall.html>

*Untracking files*

* <http://stackoverflow.com/questions/6964297/untrack-files-from-git>


---


<a id="ipynb"></a>

## Setting up a remote IPython Notebook server with a self-signed SSL certificate

__Steps (simplified):__

*Server side:*

* Create an IPython-Notebook profile `nbserver`.
* Create a self-signed SSL certificate for the secure connection between the client and the server.
* Make relevant changes to the notebook-server configuration file.
* Start the notebook server using the `nbserver` profile.
    * to let it run in a seperate process detached from any screen, I start it with an SSH command which executes the start command and falls to the background. Initially, I used something like `ssh -f skymlnb source ~/.bashrc_skyml; ipython notebook --profile=nbserver` from my own machine (where `skymlnb` is my nickname (host) for the hostname `imagediskserver3` or something like it. To control the process I tried to assign the socket to a temporary file to which an exit signal can be sent in order to close the socket. But this does not kill the process on the remote server. To kill it, I therefore log on to the server and kill it by process id.

*Client side (after following the server-side steps):*

* Open a two-way tunnel between the local machine and the remote server running the notebook server.
* Open a browser and point it to the host on the remote server that is running the notebook server
* Login on the notebook server using the password that was created for the SSL certificate.

---


