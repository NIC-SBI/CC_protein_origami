(This is a place holder. Code will be released after publication)

#**CoCoPOD**                                 
## Coiled-Coil Protein Origami Design platform

![Images of designed protein origami](images/for-readme.png)

CoCoPOD is capable of designing amino-acid sequences and building 3D models for arbitrary polyhedral meshes constructed from a single polypeptide chain. The edges of the polyhedron are realized as coiled-coil dimer building modules. The design strategy consists of several steps:

1.	**Specifying the polyhedral geometry**
2.	**Routing the chain**
3.	**Selecting the optimal topology and circular permutation** 
4.	**Selecting the building modules from the CC toolkit**
5.	**Building the 3D model**
6.	**Refining/validating the models via folding simulations** 

**CoCoPOD** performs the first five steps of the design process.  Scripts are provided for automatic execution of all the steps, for maximum flexibility the package functionality is also available from python code.

- [Quick Start](#tut)
- [Video Tutorials](#vtut)
- [Installation](#install)
- [Dependencies](#deps)
- [Tests](#tests)

When using this platform please cite: (TBA).

<a name="tut"></a>
##**Quick Start** 


Two full examples are provided in the [examples](examples/) subfolder. 
### APHsh
[APHsh](examples/APHsh/make_config.py) is centered on building models and also serves as an integration test. An antiparallel APH segment is built  The second 
To run the program user needs to provide an input file (`make_config.py`) containing information on the sequence of the protein origami design. The input consists of four sections:

* **model_name**, specifying the name of the protein origami design 
* **annotaded_sequence**, where the sequence is broke down into individual CC segments and linkers, for every CC segment a name should be provided
* **pairs_info**, in this section segment pairing is specified. For every pair the orientation (A for antiparallel and P for parallel) of the CC dimer should be provided, along with the name of CC dimer structure template file and name of the chains in the model structure.

An example of the input file can be found under [examples/APHsh/make_config.py](examples/APHsh/make_config.py). Two APH segments are connected by a linker forming a covalently linked CC dimer. 
The models can be built and viewed by by typing in the terminal

	cd cocopod/examples/APHsh	
	doit N_fold=1 N_homology=3
	chimera */03-*.pdb

Where `N_fold` is the number of independent folding simulations and `N_homology` the number of independent homology refinements of each folding simulation. The final number of models built thus equals `N_fold`*`N_homology`.

### TET
The [TET example](examples/TET/TET.ipynb) contains a complete tutorial on designing protein origami polyhedral. The tutorial is presented in the from of a python notebook.  The demonstrated steps include loading a geometry (ply) file, enumerating all the typologies and circular permutations, choosing the best topology and constructing and evaluating 3D models. The notebook can be opened by:

	cd cocopod/examples/TET	
	jupyter notebook TET.ipynb
		 

<a name="vtut"></a>
##**Video Tutorials** 

Several video tutorials (screen casts) demonstrate how to effectively use **CoCoPOD**:

[Part 1: installation](https://www.youtube.com/watch?v=lvTj_qRppME)

[![Part 1: installation](http://img.youtube.com/vi/lvTj_qRppME/2.jpg)](http://www.youtube.com/watch?v=lvTj_qRppME)

[Part 2: APHsh model building](https://www.youtube.com/watch?v=1Qa85p165Bk)

[![Part 2: APHsh model building](http://img.youtube.com/vi/1Qa85p165Bk/3.jpg)](http://www.youtube.com/watch?v=1Qa85p165Bk)

[Part 3: _de novo_ tetrahedron design](https://www.youtube.com/watch?v=-aD7mz4-XeY)

[![Part 3: de novo tetrahedron design](http://img.youtube.com/vi/-aD7mz4-XeY/3.jpg)](http://www.youtube.com/watch?v=-aD7mz4-XeY)

 


<a name="install"></a>
##**Installation**

Using the [Anaconda](https://www.continuum.io/downloads) python distribution is recommended as it simplifies installing further dependencies. [Miniconda](http://conda.pydata.org/miniconda.html) also works nicely. [Git for windows](https://git-scm.com/download/win) is also recommended. Dependencies can be installed by running:

	conda install numpy scipy pandas ipython ipython-notebook ipywidgets pyyaml xlrd biopython
	#set the modeller liscence key
	export KEY_MODELLER=XXXX
	conda install -c salilab modeller
	conda install -c omnia mdtraj 
	pip install plyfile doit #not available in conda


respectively. [Chimera](https://www.cgl.ucsf.edu/chimera/download.html) has to be installed separately. Chimera must be available on the system path.
Note: Currently there are some problems installing Modeller via conda on windows. In case of problems use the [standalone installer](https://salilab.org/modeller/download_installation.html).

**CoCoPOD** is available on [github](https://github.com/NIC-SBI/protein_origami). The files can be cloned to any location. To install the package run:

	git clone https://github.com/NIC-SBI/protein_origami.git cocopod	
	cd cocopod
	python setup.py develop

Alternatively a zip file can be download, extracted and installed with `python setup.py develop`.


<a name="deps"></a>
##**Dependencies**
The package requires Python 2.7 or Python 3.3+ with numpy and pandas and works on Windows and Linux. On Windows a bash enviorment is recommended and can be obtained by installing [git for windows](https://git-scm.com/download/win).  

Other dependencies:

* [Modeller](https://salilab.org/modeller/)
* [Chimera](https://www.cgl.ucsf.edu/chimera)
* [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) (only python 3+)
* [MdTraj](http://mdtraj.org)
* [plyfile](https://github.com/dranjan/python-plyfile)
* Numpy, Scipy, pandas, ipython, ipywidgets

Testing:

* [py.test](http://docs.pytest.org/en/latest/)
* [pytest-xdist](https://pypi.python.org/pypi/pytest-xdist) (optional)

<a name="tests"></a>
##**Tests**
Installation can be tested by executing `py.test`, which checks if core modules of the software are working appropriately and all dependencies have been installed.

	conda install pytest
	cd cocopod	
	py.test	



