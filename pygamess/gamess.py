#!/usr/bin/env python
# -*- encoding:utf-8 -*-

import openbabel as ob
import pybel
from tempfile import mkdtemp
from os import removedirs, unlink, system, environ, path, getcwd, chdir, system
from shutil import rmtree
import re
import string
from random import choice


def randstr(n):
    """make a random string"""
    return u''.join(choice('abcdefghijklmnopqrstuvwxyz') for i in xrange(n))


class GamessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Gamess(object):
    """GAMESS WRAPPER"""

    def __init__(self, gamess_path=None):
        self.tempdir = mkdtemp()
        self.debug = environ.get('debug', False)
        self.err_lines = 10

        if self.debug: print self.tempdir
        if gamess_path == None:
            try:
                gamess_path = filter(lambda f: path.isfile(f),
                                     [path.join(d, 'rungms') for d in environ['PATH'].split(':')])[0]
            except IndexError:
                print "rungms not found"
                exit()

        self.gamess = gamess_path
        self.jobname = ''
        self.cwd = getcwd()
        self.contrl = {'scftyp': 'rhf', 'runtyp': 'energy'}
        self.basis = {'gbasis': 'sto', 'ngauss': '3'}
        self.statpt = {'opttol': '0.0001', 'nstep': '20', }
        self.system = {'mwords': '30'}
        self.cis = {'nstate': '1'}
        self.obc = ob.OBConversion()
        self.obc.SetInAndOutFormats("gamout", "gamin")

    def run(self, mol):
        is_pybel = False
        if isinstance(mol,pybel.Molecule):
            mol = mol.OBMol
            is_pybel = True

        self.jobname = randstr(6)

        err_re = re.compile('^( \*\*\*|Error:)')
        eng_re = re.compile('^                       TOTAL ENERGY =')
        if self.basis['gbasis'] in ['am1', 'pm3', 'mndo']:
            eng_re = re.compile(' FINAL .+ ENERGY IS')

        gamin = self.write_file(mol)
        gamout = self.tempdir + "/" + self.jobname + ".out"

        chdir(self.tempdir)
        system("%s %s> %s  2> /dev/null" % (self.gamess, self.jobname, gamout))

        err_message = ""
        total_energy = 0

        with open(gamout, "r") as f:
            err_count = 0
            for l in f:
                if eng_re.match(l):
                    if self.basis['gbasis'] in ['am1', 'pm3', 'mndo']:
                        total_energy = float(l.split()[4])
                    else:
                        total_energy = float(l.split('=')[1])
                if err_re.match(l):
                    err_count = self.err_lines
                if err_count > 0:
                    err_message += l
                    err_count -= 1

        # エラーが出たのでstringを渡したらなおった
        # TypeError: in method 'OBConversion_ReadFile', argument 3 of type 'std::string'
        new_mol = ob.OBMol()
        s = open(gamout).read()
        self.obc.ReadString(new_mol, s)

        # singlepoint だと数値が入んないので対応
        new_mol.SetEnergy(total_energy)

        chdir(self.cwd)
        if not self.debug:
            unlink(gamin)
            unlink(gamout)

        if len(err_message) > 0:
            raise GamessError(err_message)
        else:
            if is_pybel:
                new_mol = pybel.Molecule(new_mol)
            return new_mol

    def print_header(self):
        """ gamess header"""

        header = ""
        header += self.print_section('contrl')
        header += self.print_section('basis')
        header += self.print_section('system')
        if self.contrl['runtyp'] == 'optimize':
            header += self.print_section('statpt')
        if self.contrl.get('citype', None) == 'cis':
            header += self.print_section('cis')
        return header

    def print_section(self, pref):
        d = getattr(self, pref)
        section = " $%s " % pref
        for k, v in d.iteritems():
            section += "%s=%s " % (k, v)
        section += " $end\n"
        return section

    def input(self, mol):
        if isinstance(mol,pybel.Molecule):
            mol = mol.OBMol
        self.contrl['mult'] = mol.GetTotalSpinMultiplicity()
        self.contrl['icharg'] = mol.GetTotalCharge()
        gamin_tmp = self.obc.WriteString(mol)
        h = self.print_header()
        return gamin_tmp.replace(" $CONTRL COORD=CART UNITS=ANGS $END\n", h[:-1])

    def write_file(self, mol):
        gamess_input_str = self.input(mol)
        gamess_input_file = self.tempdir + "/" + self.jobname + ".inp"
        with open(gamess_input_file, "w") as f:
            f.write(gamess_input_str)
        return gamess_input_file

    def __del__(self):
        if not self.debug:
            rmtree(self.tempdir)

    def basis_type(self, basis_type):
        if basis_type in ["sto3g", "STO-3G"]:
            self.basis = {'gbasis': 'sto', 'ngauss': '3'}
        elif basis_type in ["631g", "6-31G(d)", "6-31G"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1'}
        elif basis_type in ["631gdp", "6-31G(d,p)"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1'}
        elif basis_type in ["631+gdp", "6-31G+(d,p)"]:
            self.basis = {'gbasis': 'n31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1', 'diffsp': '.t.', }
        elif basis_type in["am1", "AM1"]:
            self.basis = {'gbasis': 'am1'}
        elif basis_type in["pm3", "PM3"]:
            self.basis = {'gbasis': 'pm3'}
        elif basis_type in["mndo", "MNDO"]:
            self.basis = {'gbasis': 'mndo'}
        else:
            print "basis type not found"
        return self.basis

    def run_type(self, runtype):
        self.contrl['runtyp'] = runtype

    def scf_type(self, scftype):
        self.contrl['scftyp'] = scftype


if __name__ == '__main__':

    g = Gamess()
   #g.basis_type('631gdp')
    #g.run_type('')
    obc = ob.OBConversion()
    obc.SetInFormat("mol")

    mol = ob.OBMol()
    next = obc.ReadFile(mol, "examples/CID_674.sdf")
    print g.input(mol)
    try:
        newmol = g.run(mol)
    except GamessError, gerr:
        print gerr.value

    print newmol.GetEnergy()
    print [(obatom.GetIdx(), obatom.GetType(), obatom.GetPartialCharge()) for obatom in ob.OBMolAtomIter(newmol)]
