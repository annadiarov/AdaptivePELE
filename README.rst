============
AdaptivePELE
============


|MIT license| |GitHub release| |PyPI release| |Conda release| |DOI|

AdaptivePELE is a Python module to perform enhancing sampling of molecular
simulation built around the Protein Energy Landscape Exploration method (`PELE <https://pele.bsc.es/pele.wt>`_) developed in the Electronic and Atomic Protein Modelling grop (`EAPM <https://www.bsc.es/discover-bsc/organisation/scientific-structure/electronic-and-atomic-protein-modeling-eapm>`_) at the Barcelona Supercomputing Center (`BSC <https://www.bsc.es>`_).

Updates of this fork respect the original repo
----------------------------------------------
Improvements in MDSimulation workflow, particularly in the equilibration step.

**NVT Equilibration**

- *Previous*:
    - It runs a constrained NVT simulation at the goal temperature of length ``equilibrationLengthNVT``. Temperature is kept constant along all the equilibration steps.
- *Current*:
    - It runs consecutive constrained NVT simulations to warm up the system, so the temperature of the system increases gradually.
    - In each simulation, we increase the temperature by 5K (default). This temperature_step is called ``temperatureStepNVTEquilibration``
    - Each simulation has a length equal to ``equilibrationLengthNVT / n_NVT_temp_increments`` where ``n_NVT_temp_increments = 1 + (temperature - initialTemperatureNVTEquilibration) / temperatureStepNVTEquilibration``

**NPT Equilibration**

- *Previous*:
    - It runs a constrained NPT simulation at the goal temperature of length ``equilibrationLengthNPT``. Constraints are reduced compared to NVT but keep constant along the whole equilibration.
- *Current*:
    - It runs consecutive constrained NPT simulations to reduce the constraints of the system gradually.
    - In each simulation we reduce the constraints 0.5 kcal/(mol*A^2), this temperature_step is called ``constraintStepNPTEquilibration``
    - Each simulation has a length equal to ``equilibrationLengthNPT / n_NPT_constr_reductions`` where ``n_NPT_constr_reductions = 1 + initial_constraints / constraintStepNPTEquilibration``, where initial_constraints is set to ``constraintsNVT`` value.

**Salt Concentration**

- *Previous*:
    - It only neutralized the system adding counter ions.
- *Current*:
    - It adds counter ions to neutralize the system, and computes how many Na+ and Cl- are needed to reach a specific `saltConcentration` (eg. 0.15M).
    - To compute the number of ions needed, we run tleap twice: the first time to solvate the system and extract the volume from the leap.log; the second time to set the correct number of Na+ and Cl- ions.
    - **Warnings:**
        - Adding ions significantly increases tleap running time (ie. to place ~60 ions needs ~1h of computation)
        
**Covalent ligands**

- *Previous*:
    - It only allowed to simulate non-covalent ligands (no attachment to a residue of the protein)
- *Current*:
    - The code has been modified in order to allow loading the frcmod of lib files from possible covalent ligands in the user's system.
    - The force field files of the covalent ligand(s) must be prepared manually (for now) and added to the "constants/MDtemplates" folder from the AdaptivePELE repository.
    - **Preparing protocol:**
    
        - You must first isolate a PDB file with the covalent ligand and execute antechamber like this ``antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -c bcc -pf y -nc 0``
        - The electrostatic potential (ESP) charges of the atoms from the covalent ligand must be calculated; ``$AMBERbin/sqm -O -i sqm.in -o sqm.out``
        - The frcmod file is calculated from the mol2 file; ``parmchk2 -i $1.mol2 -f mol2 -o $1.frcmod``
        - The lib file is calculated from the mol2 file with tleap and the following commands: ``source leaprc.gaff``, ``LIG = loadmol2 LIG.mol2``, ``check LIG``, ``loadamberparams LIG.frcmod``, ``saveoff LIG amber_LIG.lib``
        - Then, the frcmod and lib files must be manually modified to take into account that the edges where the ligand should be bound to the rest of the protein must be included (the files must have a name like the following; "amber\_LIG")
        - The main modifications are: change the indexes of the `!entry.CZ3.unit.connect array int`` and ``!entry.LYP.unit.residueconnect`` sections of the lib file to the atoms at the covalent edges (the N and C atoms of the peptide bond of the backbone), change the ``restype`` attribute from the ``!entry.LYP.unit.residues`` section in the lib file from "?" to "p", adding the bond, angle, and dihedral parameters of the linking atoms of the covalent ligand to the frcmod file (a good and standard practice is to add the bond, angle, and dihedral parameter values of the peptide bonds from the conventional amino acids)

To control these behaviours using the control_file, we included the following simulation params for MDs.

- **NVT Equilibration**
    - ``temperatureStepNVTEquilibration``: Increment of temperature during the warm-up in NVT equilibrartion

        Default: 5 (type ``float``)
    - ``initialTemperatureNVTEquilibration``: Initial temperature at which we start heating the system during NVT equilibrartion

        Default: 5 (type ``float``)

- **NPT Equilibration**
    - ``constraintStepNPTEquilibration``: Reduction of constraints in each in NPT equilibration simulation

        Default: 0.5 (type ``float``)
    - ``finalConstraintValueNPTEquilibration``: Target constraint to be reached in NPT equilibration. If set to 0, all constraints will be released.

        Default: 0 (type ``float``)
    - ``lengthUnconstrainedNPTEquilibration``: Number of steps for to be added to the last NPT simulation, that is the one with constraints ``finalConstraintValueNPTEquilibration``.

        - Default: 500000 (type ``int``) (1ns)
        - Warning: In case you have defined a ligand box, the production step might not be equivalent to this unconstrained NPT step, since ligand constraints are defined in the production step. This should be checked if someone needs it. An easy solution to solve this could be defining the ``finalConstraintValueNPTEquilibration`` to 0.5 or 0.9 instead of 0.
    - DEPRECATED ``constraintsNPT``: This parameter from the original implementation is no longer used, so I removed it.

- **Salt Concentration**
    - `saltConcentration`: Salt concentration to be set up in the topology by adding Na+ and Cl- ions, after neutralizing the system.

        - Default: 0 (type ``float``)

- **Box type**
    - `useCubicBox`: If True, the box will be cubic, if False, the box will be octahedral. *Note*: Octahedral box is still in development.

        - Default: True (type ``bool``)

Usage
-----

AdaptivePELE is called with a control file as input
parameter. The control file is a json document that contains 4 sections:
general parameters, simulation parameters, clustering parameters and spawning
parameters. The first block refers to general parameters of the adaptive run,
while the other three blocks configure the three steps of an adaptive sampling
run, first run a propagation algorithm (simulation), then cluster the
trajectories obtained (clustering) and finally select the best point to start
the next iteration (spawning).

An example of usage::

    python -m AdaptivePELE.adaptiveSampling controlFile.conf

Installation
------------

Install from source, you need to install and compile cython files in the base folder with::

    git clone https://github.com/annadiarov/AdaptivePELE.git
    cd AdaptivePELE
    python setup.py build_ext --inplace

Also, if AdaptivePELE was not installed in a typical library directory, a common option is to add it to your local PYTHONPATH::

    export PYTHONPATH="/location/of/AdaptivePELE:$PYTHONPATH"

**Update from this fork (Recommended)**
You can create a conda environment to run the MDs with cuda using the following command::

    conda env create -f conda_recipe/openmm_77_adaptive_md.yml

The code provided here was tested only with openmm 7.7.0, so it is recommended to use this version to avoid any compatibility issues.
Also, if you have problems building the environment, it's probably because of the cython version. We build the model using 0.29.24 version.

Documentation
-------------

The documentation for AdaptivePELE can be found `here <https://bsc-cns-eapm.github.io/AdaptivePELE/>`_


Contributors
------------
`Daniel Lecina <https://github.com/lecina>`_, `Joan Francesc Gilabert <https://github.com/cescgina>`_, `Oriol Gracia <https://github.com/OriolGraCar>`_, `Daniel Soler <https://github.com/danielSoler93>`_

Mantainer
---------
Joan Francesc Gilabert (cescgina@gmail.com)

Citation 
--------

AdaptivePELE is research software. If you make use of AdaptivePELE in scientific publications, please cite it. The BibTeX reference is::

    @article{Lecina2017,
    author = {Lecina, Daniel and Gilabert, Joan Francesc and Guallar, Victor},
    doi = {10.1038/s41598-017-08445-5},
    issn = {2045-2322},
    journal = {Scientific Reports},
    number = {1},
    pages = {8466},
    pmid = {28814780},
    title = {{Adaptive simulations, towards interactive protein-ligand modeling}},
    url = {http://www.nature.com/articles/s41598-017-08445-5},
    volume = {7},
    year = {2017}
    }


.. |MIT license| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://lbesson.mit-license.org/


.. |GitHub release| image:: https://img.shields.io/github/release/AdaptivePELE/AdaptivePELE.svg
    :target: https://github.com/AdaptivePELE/AdaptivePELE/releases/

.. |PyPI release| image:: https://img.shields.io/pypi/v/AdaptivePELE.svg
    :target: https://pypi.org/project/AdaptivePELE/

.. |DOI| image:: https://zenodo.org/badge/DOI/10.1038/s41598-017-08445-5.svg
  :target: https://doi.org/10.1038/s41598-017-08445-5
  
.. |Conda release| image:: https://anaconda.org/nostrumbiodiscovery/adaptive_pele/badges/version.svg
  :target: https://anaconda.org/NostrumBioDiscovery/adaptive_pele
