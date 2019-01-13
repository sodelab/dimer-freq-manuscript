# SI for CO2 dimer anharmonic vibrations

This project produces all of the data from the Anharmonic vibrational structure of the carbon dioxide dimer with a many-body potential energy surface journal article. The project solves the vibrational Schrodinger equation for the CO<sub>2</sub> monomer and dimer using vibration structure theory and a many-body potential energy surface.  

## Getting Started

These instructions will get your copy of the project up and running on your local machine. 

### Prerequisites

The most efficient and straightforward way to run this project is to install [Docker](https://docs.docker.com/) on your machine. This step allows for a docker container to be built with the correct dependencies for the codes. A [docker file](dockerfile) is located in the project folder, and it contains instructions on how to build the project container and install the necessary software. Instructions can be found in the `Installing` section. If you'd rather not use Docker you'll need to install the following packages and repos before running any of the repository scripts. These include:

```
cmake
wget
gcc
gfortran
libopenblas-dev
liblapack-dev
libgsl2
vim
bc
libgsl0-dev
vibrational structure theory programs (MaVi, SINDO, and NITROGEN)
```

Many of these are most likely already installed on your favorite machine. If not however, you'll need to install them along with the all three of the vibrational structure theory programs:

```
[MaVi](https://github.com/keceli/MaVi.git)
[SINDO]()
[NITROGEN](https://www.colorado.edu/nitrogen/)
```

These three programs need to be placed in the [tools](tools) directory of the downloaded repository and installed individually. In the tools directory, there are four other directories each containing code that needs to be individually installed in the following order: `library`, `optimize`, `evaluate` and `anharmonic`. Once all of the programs are installed, proceed to the [step](step) directory to execute the project.

### Installing

This section describes how to install the project using Docker. If you'd prefer to not use Docker, read through the `Prerequisites` section above. 

1. Download the [docker file](dockerfile) build instructions from this repository. 

2. Make sure to download and install the Docker software onto your machine. Installation links can be found at the following [webpage](https://docs.docker.com/). 

2. Once downloaded, you can build an image using the following command:

   `docker build . -t <MYIMAGE>`,
   
   where <MYIMAGE> is your chosen image name. The image may take a few minutes to build, as well as download and install all      of the necessary software, including this [repository](), and the other vibrational structure theory programs. 
   
3. To run the image, enter the following command:

   `docker run -it <MYIMAGE> bash`.
   
   The current project can be found in the `/container` folder. Using the command line, navigate to this folder and then          proceed to the `step` directory. Here you should find a [README](step/README) with instructions on how to execute the          project commands and generate the manuscript data.

Enjoy!

## Authors

* **Olaseni Sode** - *Initial work* - [oosode](https://github.com/sodelab)

## Acknowledgments

* Thanks [Murat Keçeli](https://github.com/keceli)
