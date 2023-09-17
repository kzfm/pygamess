# pygamess

<span class="title-ref">pygamess</span> is a GAMESS wrapper for Python

## Requirements

-   Python 3.7 or later (pygamess &lt;0.5 supports only Python2)
-   RDKit &gt;= 2020.03.5
-   GAMESS > Jun2020R1
-   ruamel.YAML

## Setup

    $ pip install pygamess

set GAMESS\_HOME environment in your .bashrc or .zshrc:

    $ export GAMESS_HOME=/usr/local/gamess

*** 
Windows/Mac users can obtain the pre-compiled binary executables from [GAMESS download site](https://www.msg.chem.iastate.edu/gamess/download.html).
But Linux users need to compile the souce code. 
***

## Test

    $ pytest

## Basic Usage

### Single point calculation

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> r = g.run(m)
    >>> r.total_energy
    -152.127991054

Or use rdkit directly:

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> m = Chem.MolFromSmiles("CCO")
    >>> m = Chem.AddHs(m)
    >>> AllChem.EmbedMolecule(m)
    0
    >>> AllChem.UFFOptimizeMolecule(m,maxIters=200)
    0
    >>> g = Gamess()
    >>> r = g.run(m)
    >>> r.total_energy
    -152.1279910526

The GamessOut object(r) contains the results of the GAMESS calculation
and the RDKit Chem object with the calculation results:

    >>> r.total_energy
    -152.1279910526
    >>> r.HOMO
    -0.3453
    >>> r.nHOMO
    -0.3978
    >>> r.LUMO
    0.5594
    >>> r.nLUMO
    0.6127
    >>> r.dipole_moment
    [0.681619, -0.605188, 1.146253, 1.464497]
    >>> r.orbital_energies
    [-20.2521, -11.0932, -11.0402, -1.286, -0.9614, -0.7909, -0.6269, -0.5716, -0.5347, -0.4976, -0.4705, -0.3978, -0.3453, 0.5594, 0.6127, 0.6639, 0.69, 0.7002, 0.7388, 0.7549, 0.7852]
    >>> r.mulliken_charges
    [-0.171226, 0.024351, -0.298162, 0.049615, 0.055831, 0.061924, 0.042714, 0.061112, 0.173842]
    >>> r.lowdin_charges
    [-0.089262, 0.062464, -0.212689, 0.022024, 0.02752, 0.031644, 0.010235, 0.026611, 0.121453]

The RDKit Chem object has the same information to store these data into
SDF:

    >>> r.mol
    <rdkit.Chem.rdchem.Mol object at 0x7ffd48217f30>
    >>> r.mol.GetProp("total_energy")
    '-152.12799105260001'
    >>> r.mol.GetProp("HOMO")
    '-0.3453'
    >>> r.mol.GetProp("nHOMO")
    '-0.39779999999999999'
    >>> r.mol.GetProp("LUMO")
    '0.55940000000000001'
    >>> r.mol.GetProp("nLUMO")
    '0.61270000000000002'
    >>> r.mol.GetProp("dipole_moment")
    '1.4644969999999999'
    >>> r.mol.GetProp("dx")
    '0.68161899999999997'
    >>> r.mol.GetProp("dy")
    '-0.60518799999999995'
    >>> r.mol.GetProp("dz")
    '1.146253'
    >>> r.mol.GetProp("orbital_energies")
    '-20.2521 -11.0932 -11.0402 -1.286 -0.9614 -0.7909 -0.6269 -0.5716 -0.5347 -0.4976 -0.4705 -0.3978 -0.3453 0.5594 0.6127 0.6639 0.69 0.7002 0.7388 0.7549 0.7852'
    >>> for a in r.mol.GetAtoms():
    ...   print("{}:\t{:.4f}\t{:.4f}".format(a.GetSymbol(), float(a.GetProp("mulliken_charge")), float(a.GetProp("lowdin_charge"))))
    ... 
    C:  -0.1712 -0.0893
    C:  0.0244  0.0625
    O:  -0.2982 -0.2127
    H:  0.0496  0.0220
    H:  0.0558  0.0275
    H:  0.0619  0.0316
    H:  0.0427  0.0102
    H:  0.0611  0.0266
    H:  0.1738  0.1215

### Geometry optimization

Set the run\_type as 'optimize'. This optimization process updates the
coordinates of the molecule:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type('optimize')
    >>> r = g.run(m)
    >>> r.total_energy
    -152.1330661028
    >>> original_conf = m.GetConformer(0)
    >>> optimized_conf = r.mol.GetConformer(0)
    >>> for c in original_conf.GetPositions():
    ...   print(c)
    ... 
    [ 0.91206647 -0.11944851 -0.1294722 ]
    [-0.47153193  0.42043351  0.21118521]
    [-1.44831334 -0.21539324 -0.56715297]
    [ 0.9650486  -1.2050043   0.09903891]
    [1.67955732 0.41189183 0.47186063]
    [ 1.12654515  0.0378618  -1.20780162]
    [-0.67995345  0.27545635  1.29552482]
    [-0.49851483  1.50993204 -0.00177488]
    [-1.58490399 -1.11572948 -0.17140791]
    >>> for c in optimized_conf.GetPositions():
    ...   print(c)
    ... 
    [ 0.91442972 -0.13086468 -0.12174822]
    [-0.48373921  0.42850882  0.23169745]
    [-1.54145595 -0.17763397 -0.52424945]
    [ 0.97874385 -1.18768306  0.12185969]
    [1.67944907 0.39674369 0.43935649]
    [ 1.11452802 -0.00923776 -1.18177409]
    [-0.65534535  0.31918718  1.31002371]
    [-0.51450285  1.49751073  0.00304491]
    [-1.4921073  -1.13653094 -0.2782105 ]

Or pass the options to a constractor::

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess(options={"contrl":{"runtyp":"optimize"}})
    >>> r = g.run(m)
    >>> r.total_energy
    -152.1330661279

### Calculating IR spectra after optimization

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCC(=O)O")
    >>> g = Gamess()
    >>> g.run_type('optimize', hessend=True)
    >>> r = g.run(m)
    >>> r.ir_spectra
    [('1.078', '0.058962'), ('0.361', '0.000420'), ('0.119', '0.000239'), ('0.703', '0.027425'), ('1.694', '0.025676'), ('3.212', '0.064166'), ('80.809', '0.030877'), ('220.540', '0.023447'), ('255.386', '0.510706'), ('362.039', '1.450661'), ('461.248', '0.076788'), ('590.483', '0.244231'), ('743.957', '0.128188'), ('915.147', '0.079289'), ('941.690', '0.122604'), ('1163.260', '0.213250'), ('1262.233', '0.180680'), ('1283.531', '0.070173'), ('1451.777', '0.939433'), ('1514.198', '0.042036'), ('1564.520', '3.659319'), ('1627.441', '0.509608'), ('1736.037', '0.011788'), ('1799.935', '0.021852'), ('1830.887', '0.076255'), ('1834.826', '0.085977'), ('2153.439', '1.875334'), ('3572.386', '0.028347'), ('3609.011', '0.023501'), ('3730.553', '0.009088'), ('3752.383', '0.015918'), ('3758.396', '0.003956'), ('4272.454', '0.229192')]

    ### GAMESS ### 
    #  DFT ANALYTIC HESSIAN PRESENTLY HAS 5 RESTRICTIONS:
    #  $CONTRL: SCFTYP MUST BE EITHER RHF OR UHF
    #  $CONTRL: POINT GROUP SYMMETRY NOT ALLOWED, SET NOSYM=1
    #     $SCF: AO INTEGRAL DIRECT: SET DIRSCF=.TRUE.
    #    $CPHF: AO INTEGRAL DRIVEN: SET CPHF=AO
    #  AND THE FUNCTIONAL MUST NOT BE OF META-GGA TYPE.

### Changing basis sets

Use basis\_sets method:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type = "optimize"
    >>> g.run(m).total_energy
    -152.127991054
    >>> g.basis_sets("3-21G")
    >>> g.run(m).total_energy
    -153.2170653562
    >>> g.basis_sets("6-31G")
    >>> g.run(m).total_energy
    -154.0054866151
    >>> g.basis_sets("6-31G*")
    >>> g.run(m).total_energy
    -154.0702703669
    >>> g.basis_sets("6-31G**")
    >>> g.run(m).total_energy
    -154.0843823698

Or edit the basis attribute directly:

    >>> g.options({'basis':{'gbasis': 'sto', 'ngauss': '3'}})
    >>> g.run(m).total_energy
    -152.127991054

### Changing the electronic state

The ground state (default):

    >>> from pygamess.utils import rdkit_optimize
    >>> from pygamess import Gamess
    >>> g = Gamess()
    >>> g.basis_sets("6-31G*")
    >>> m = rdkit_optimize("CCO")
    >>> g.run_type("optimize")
    >>> r = g.run(m)
    >>> r.total_energy
    -154.0755757352

The cationic state:

    >>> g.scf_type("uhf")
    >>> g.charge(1)
    >>> g.multiplicity(2)
    >>> r = g.run(m)
    >>> r.total_energy
    -153.7367666449

Or:

    >>> g.options({"contrl":{"icharg": 1, "mult": 2, "scftyp": "uhf"}})
    >>> r = g.run(m)
    >>> r.total_energy
    -153.7367666449

The anionic state:

    >>> g.options({"contrl":{"icharg": -1, "mult": 2, "scftyp": "uhf"}})
    >>> r = g.run(m)
    >>> r.total_energy
    -153.9302151707

The triplet state:

    >>> g.options({"contrl":{"icharg": 0, "mult": 3, "scftyp": "uhf"}})
    >>> r = g.run(m)
    >>> r.total_energy
    -153.9581403463

### DFT calculation

B3LYP/6-31G\*:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type("optimize")
    >>> g.basis_sets("6-31G*")
    >>> g.dft_type("B3LYP")
    >>> g.run(m).total_energy
    -154.9387962055

M062X/6-31G\*\*:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.run_type("optimize")
    >>> g.dft_type("M06-2X")
    >>> g.basis_sets("6-31G**")
    >>> g.run(m).total_energy
    -154.9636095207

### PCM calculation

Pygamess currently only supports CPCM, but will support IEFPCM in the
future:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> g.basis_sets("6-31G*")
    >>> g.pcm_type("water")
    >>> g.run_type("optimize")
    >>> r = g.run(m)
    >>> r.total_energy
    -154.0824604616
    >>> r.internal_energy
    -154.0748367584
    >>> r.delta_internal_energy
    0.0
    >>> r.electrostatic_interaction
    -0.0076237032
    >>> r.pierotti_cavitation_energy
    0.0
    >>> r.dispersion_free_energy
    0.0
    >>> r.repulsion_free_energy
    0.0
    >>> r.total_interacion
    -0.0076237032

### TDDFT calculation

The example compound(Methyl yellow) was downloaded from [PubchemQC project](http://pubchemqc.riken.jp/cgi-bin/molecularquery.py?name=methyl+yellow).

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> m = Chem.MolFromMolFile("examples/methyl_yellow.mol", removeHs=False)
    >>> g = Gamess()
    >>> g.dft_type("b3lyp", tddft=True)
    >>> g.basis_sets("6-31G*")
    >>> r = g.run(m)
    >>> r.uv_spectra # (exitation ev, oscillator strength)
    [('2.629', '0.0000'), ('3.217', '0.9349'), ('4.209', '0.0066'), ('4.263', '0.0020'), ('4.424', '0.1041'), ('4.779', '0.1068'), ('4.913', '0.0563'), ('4.940', '0.0001'), ('5.051', '0.0000'), ('5.430', '0.0006')]

- ref: [独習 量子化学計算(Self-study Quantum Chemical Calculations)](https://www.amazon.co.jp/gp/product/B0863C799Z/)

### NMR spectra calculation

Optimizing the compound:

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> g = Gamess()
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("C=CCBr")
    >>> g.run_type("optimize")
    >>> g.dft_type("b3lyp")
    >>> g.basis_sets("6-31G*")
    >>> r = g.run(m)
    >>> with open("examples/C=CCBr.mol", "w") as f:
    ...   f.write(Chem.MolToMolBlock(r.mol))
    ... 

NMR spectra calculation (It takes a long time):

    >>> from pygamess import Gamess
    >>> from rdkit import Chem
    >>> m = Chem.MolFromMolFile("examples/C=CCBr.mol", removeHs=False)
    >>> g = Gamess(num_cores=1) # PARALLEL EXECUTION IS NOT ENABLED.
    >>> g.basis_sets("6-31G*")
    >>> g.run_type("nmr")
    >>> r = g.run(m)
    >>> r.isotropic_shielding
    [79.8218, 68.6661, 157.7233, 2476.7501, 27.0851, 27.2072, 26.4652, 28.7654, 28.7932]

    # NMR MAY BE COMPUTED ONLY FOR SCFTYP=RHF,
    # NO CORRELATION OPTION (DFTTYP, CITYP, CCTYP, MPLEVL) MAY BE CHOSEN
    # NO SEMI-EMPIRICAL OPTION (GBASIS=AM1/PM3/MNDO) MAY BE CHOSEN
    # DIRECT AO INTEGRAL CALCULATION (DIRSCF) IS NOT ENABLED,
    # AND/OR PARALLEL EXECUTION IS NOT ENABLED.

- ref: [独習 量子化学計算(Self-study Quantum Chemical Calculations)](https://www.amazon.co.jp/gp/product/B0863C799Z/)

### Printing GAMESS INPUT

use input method:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> m = rdkit_optimize("CO")
    >>> g = Gamess()
    >>> print(g.input(m))
    $contrl scftyp=rhf runtyp=energy $end
    $basis gbasis=sto ngauss=3 $end
    $system mwords=100 $end
    $DATA
    6324
    C1
    C      6.0     -0.3577002260    0.0075902163   -0.0214817423 
    O      8.0      0.9087355734   -0.5349924519   -0.2611189822 
    H      1.0     -0.5468334701    0.0717914414    1.0721087268 
    H      1.0     -0.4337681128    1.0193437527   -0.4757947304 
    H      1.0     -1.1269974200   -0.6479305528   -0.4789564639 
    H      1.0      1.5565636556    0.0841975943    0.1652431920 
    $END

### Parsing a GAMESS output file

use gparse:

    >>> from pygamess.gamout_parser import gparse
    >>> r = gparse("/somewhere/gamess.out")

### Persistence of calculation results

save sdf file:

    >>> from pygamess import Gamess
    >>> from pygamess.utils import rdkit_optimize
    >>> from rdkit import Chem
    >>> m = rdkit_optimize("CCO")
    >>> g = Gamess()
    >>> r = g.run(m)
    >>> w = Chem.SDWriter("CCO.sdf")
    >>> w.write(r.mol)
    >>> w.close()

load from sdf file:

    >>> from pygamess.utils import sdf2gamout
    >>> r = sdf2gamout("CCO.sdf")
    >>> r = rs[0]
    >>> r.HOMO
    -0.3453
    >>> r.mulliken_charges
    ['-0.17122599999999999', '0.024351000000000001', '-0.29816199999999998', '0.049614999999999999', '0.055830999999999999', '0.061924', '0.042714000000000002', '0.061112', '0.173842']

### Debugging pygamess

set PYGAMESS\_DEBUG environment:

    $ export PYGAMESS_DEBUG=1

This won't remove the all files generated by the GAMESS executable,
including the output files.

set logger level:

    >>> from pygamess import Gamess, logger
    >>> import logging
    >>> logger.setLevel(logging.DEBUG)
    >>> g = Gamess()
    DEBUG:pygamess.gamess:tmpdir: /var/folders/gm/4tcnnyqd09d2jt7p0dtvr28m0000gn/T/tmp889j9c7e

## History

### 0.6.9 (2023-09-18)

- Fix Mulliken/Lowdin population bug

### 0.6.8 (2023-09-10)

- Add Mulliken/Lowdin population

### 0.6.7 (2021-09-12)

- Bug fix (options)
- Support Windows pre-compiled GAMESS (#11)
- OS-dependent termination detection (#12)

### 0.6.6 (2021-08-31)

- Bug fix (SDF -> GamessOut)

### 0.6.5 (2021-08-31)

- Add SDF -> GamessOut object parser for persistence of calculation results

### 0.6.4 (2021-06-29)

- Store UV-spectra NMR-spectra and IR-spectra into Chem object

### 0.6.3 (2021-05-09)

- Support hessend (#8)
- Support TD-DFT (#2)
- Store the energy of each step during structural optimization (#7)
- Add the calculation condition into Chem object (#6)
- Support NMR calculation

### 0.6.2 (2021-05-05)

- Support options
- Add parser description

### 0.6.1 (2021-05-03)

- Fix bug (ModuleNotFoundError: No module named 'pygamess_utils')
- Change README format (rst -> md)

### 0.6.0 (2021-05-02)

-   Support DFT calculation
-   Support PCM calculation (C-PCM only)
-   Improve the parser
-   Support logger levels
-   Change method name from "basis\_set" to "basis\_sets"

### 0.5.0 (2020-09-13)

-   Support Python3

### 0.4.1.1 (2017-09-16)

-   Update Readme

### 0.4.1 (2017-09-16)

-   Bug fix (coordinates problem)

### 0.4.0 (2017-09-13)

-   Change the backend library from openbabel to RDKit

### 0.3.0 (2012-03-31)

-   Use internal rungms (default)
-   Add basis\_set
    method(STO-3G,3-21G,6-31G,6-311G,6-31G\*,6-31G\*\*,AM1,PM3,MNDO)
-   Constructor can accept options
-   Bug fixed (spin multiplicity)

### 0.2.2 (2012-03-30)

-   Add charge settings
-   Change Method name (gamess\_input -&gt; input)

### 0.2.1 (2012-03-23)

-   Bug fix (multiplicity setting for pybel)
-   Bug fix (print error when rungms exec failed)
-   Add document

### 0.2.0 (2012-03-06)

-   Run method accepts OBMol and Pybel-Molecule object

### 0.1.2 (2011-09-23)

-   Add CIS method (and optimization)

### 0.1.1 (2011-08-06)

-   Update document
-   Semiempical method (AM1, PM3, MNDO)
-   Add statpt option
-   Change default error print message (10 lines)

### 0.1 (2011-6-25)

-   First release
