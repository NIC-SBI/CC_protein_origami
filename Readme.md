#**PROTEIN ORIGAMI**                                 
##Platform for the design of single chain protein topological polyhedral cages 
---------------------------------------
The computational platform is capable of designing amino-acid sequences and building 3D models for arbitrary polyhedral meshes constructed from a single polypeptide chain. The edges of the polyhedron are realized as coiled-coil dimer building modules. The design strategy consists of several steps:

1. Specify the polyhedral geometry
2. Routing the chain
3. Selection of  the optimal topology and circular permutation 
4. Building 3D model 
5. Refinement/validation of the models via folding simulations 

**PROTEIN ORIGAMI** computational platform performs the first four steps of the design process. Scripts are provided for automatic execution of all the steps, however the design process can be carried out step by step.

###**Dependencies**
---------------------------------------
Scripts require Python 2.7 or Python 3.3+ with numpy and pandas.

Other dependencies:

* Modeller
* Chimera
* Snakemake (only python 3+)
* MdTraj

###**Installation**
---------------------------------------
**PROTEIN ORIGAMI** design software is available as a git repository [ppmod](https://bitbucket.org/l12/poly_modeller). The files can be cloned to any location.

Using the [Anaconda](https://docs.continuum.io/anaconda/install) python distribution is recommended as it simplifies installing further dependencies. Modeller, Snakemake and MdTraj can then easily be installed by running
>$ conda install -c salilab modeller
>
>$ conda install -c bioconda snakemake
>
>$ conda install -c omnia mdtraj 

respectively. [Chimera](https://www.cgl.ucsf.edu/chimera/download.html) has to be installed separately. 


###**Tests**
---------------------------------------
Installation can be tested by executing py.test, which checks if core modules of the software are working appropriately and all dependencies have been installed.


###**Jump start**
---------------------------------------
To run the program user needs to provide an input file containing information on the sequence of the protein origami design. The input consists of four sections:

* **_name_**, specifying the name of the protein origami design 
* **_entire\_sequence_**, where the polypeptide sequence should be inserted 
* **_segments_**, where the sequence is broke down into individual CC segments and linkers, for every CC segment a name should be provided
* **_pairs_**, in this section segment pairing is specified. For every pair the orientation (A for antiparallel and P for parallel) of the CC dimer should be provided, along with the name of CC dimer structure template file and name of the chains in the model structure.

An example of the input file can be found [here](https://bitbucket.org/l12/poly_modeller/src/aaf92e2cd01b7c84ff3a6db5359eefed6fa5d305/examples/APHsh/make_config.py?at=master&fileviewer=file-view-default). Two APH segments are connected by a linker forming a covalently linked CC dimer. 
The user is provided with a [Snakefile](https://bitbucket.org/l12/poly_modeller/src/aaf92e2cd01b7c84ff3a6db5359eefed6fa5d305/examples/APHsh/Snakefile?at=master&fileviewer=file-view-default) which allows for automatic performance of all the design steps through running *snakemake* in the command line. However, the design can then be carried out stepwise by using appropriate scripts. 

###**Terms of use**
---------------------------------------
The code is licensed under GPLv3. Copyright (C) 2016 Ajasja Ljubetič, Igor Drobnak, Jana Aupič. In published work which uses this package please cite ...
