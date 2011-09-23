.. pygamess documentation master file, created by
   sphinx-quickstart on Thu Jun 23 20:39:43 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pygamess's documentation!
====================================

基本的な使いかた
---------------------------------------

モジュールをインポートしてGamessオブジェクトを生成します::

    >>> import pygamess
    >>> g = pygamess.Gamess()

openbabelモジュールを使って分子を読み込みます::

    >>> import openbabel as ob
    >>> obc = ob.OBConversion()
    >>> obc.SetInFormat("mol")
    True
    >>> mol = ob.OBMol()
    >>> obc.ReadFile(mol, "examples/ethane.mol")
    True

gamessの実行にはrunメソッドを使ってOBMolオブジェクトを渡します。Gamess
実行時にエラー終了した場合にはGamessError例外が投げられます::

    >>> try:
    ...     newmol = g.run(mol)
    ... except GamessError, gerr:
    ...     print gerr.value
    ... 

エネルギーはGetEnergyメソッドで取得できます::

    >>> newmol.GetEnergy()
    -78.305307479999996

原子の情報はOBMolのメソッドを利用することでアクセスできます::

    >>> for obatom in ob.OBMolAtomIter(newmol):
    ...     (obatom.GetIdx(), obatom.GetType(), obatom.GetPartialCharge())
    ... 
    (1, 'C3', -0.16967199999999999)
    (2, 'C3', -0.16967199999999999)
    (3, 'HC', 0.056557999999999997)
    (4, 'HC', 0.056559999999999999)
    (5, 'HC', 0.056554)
    (6, 'HC', 0.056559999999999999)
    (7, 'HC', 0.056554)
    (8, 'HC', 0.056557999999999997)

構造最適化計算を行う場合
---------------------------------------------------

Gamessオブジェクトのrun_typeで指定します::

    >>> g.run_type('optimize')
    >>> optimized_mol = g.run(mol)
    >>> optimized_mol.GetEnergy()
    -78.306179642000004

GetEnergiesメソッドで最適化の各構造のエネルギーを取得することができるの
で、振動してないかどうかチェックすることができます。::

    >>> optimized_mol.GetEnergies()
    (-78.305307499999998, -78.306143500000005, -78.306164999999993, -78.306179400000005, -78.306179599999993)

基底関数を変更したい場合
--------------------------------------------------

basis_typeで変更できます::

    >>> g.run_type('energy')
    >>> g.basis_type('sto3g')
    {'gbasis': 'sto', 'ngauss': '3'}
    >>> mol_sto3g = g.run(mol)
    >>> mol_sto3g.GetEnergy()
    -78.305307479999996
    >>> g.basis_type('631g')
    {'gbasis': 'N31', 'ndfunc': '1', 'ngauss': '6'}
    >>> mol_631g = g.run(mol)
    >>> mol_631g.GetEnergy()
    -79.228127109699997
    >>> g.basis_type('631gdp')
    {'gbasis': 'N31', 'ndfunc': '1', 'npfunc': '1', 'ngauss': '6'}
    >>> mol_631gdp = g.run(mol)
    >>> mol_631gdp.GetEnergy()
    -79.237634701499999

又はbasisプロパティを直接変更します::

    >>> g.basis = {'gbasis': 'sto', 'ngauss': '3'}
    >>> mol_sto3g = g.run(mol)
    >>> mol_sto3g.GetEnergy()
    -78.305307479999996

GAMESSのインプットファイルを出力したい場合
----------------------------------------------------

gamess_inputメソッドが使えます

励起状態で構造最適化したい場合
----------------------------------------------------

CISによる構造最適化計算が行えます(予定)。



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

