(This is a place holder. Code will be released after the method is published.)

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

* [Modeller](https://salilab.org/modeller/)
* [Chimera](https://www.cgl.ucsf.edu/chimera)
* [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) (only python 3+)
* [MdTraj](http://mdtraj.org)
* Numpy, Scipy, pandas

Testing:
* [py.test](http://docs.pytest.org/en/latest/)

###**Installation**
---------------------------------------
Using the [Anaconda](https://docs.continuum.io/anaconda/install) python distribution is recommended as it simplifies installing further dependencies. Modeller, Snakemake and MdTraj can then easily be installed by running

	
	conda install numpy scipy pandas
	conda install -c salilab modeller
	conda install -c bioconda snakemake
	conda install -c omnia mdtraj 


respectively. [Chimera](https://www.cgl.ucsf.edu/chimera/download.html) has to be installed separately. Chimera must be available on the system path.

**PROTEIN ORIGAMI** design software is available as a git repository [protein origami](https://github.com/NIC-SBI/protein_origami). The files can be cloned to any location. Using command *cd* move to the folder poly_modeller and run 

	python setup.py install

to install the package.


###**Tests**
---------------------------------------
Installation can be tested by executing `py.test`, which checks if core modules of the software are working appropriately and all dependencies have been installed.


###**Jump start**
---------------------------------------
To run the program user needs to provide an input file containing information on the sequence of the protein origami design. The input consists of four sections:

* **_name_**, specifying the name of the protein origami design 
* **_entire\_sequence_**, where the polypeptide sequence should be inserted 
* **_segments_**, where the sequence is broke down into individual CC segments and linkers, for every CC segment a name should be provided
* **_pairs_**, in this section segment pairing is specified. For every pair the orientation (A for antiparallel and P for parallel) of the CC dimer should be provided, along with the name of CC dimer structure template file and name of the chains in the model structure.

An example of the input file can be found under `examples/APHsh/make_config.py`. Two APH segments are connected by a linker forming a covalently linked CC dimer. 
The user is provided with a snake-file which can be run by typing

>snakemake 

in the command line. However, the design can then be carried out stepwise by using appropriate scripts. 
