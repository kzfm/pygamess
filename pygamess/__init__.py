#!/usr/bin/env python
# -*- encoding:utf-8 -*-

import openbabel as ob
from tempfile import mkstemp, mkdtemp
from os import removedirs, unlink, system, environ, path, getcwd, chdir, system
import re
import string
from random import choice

DEBUG = environ.get('debug', False)


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
        if DEBUG:
            print self.tempdir
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

    def run(self, mol):
        self.jobname = randstr(6)
        obc = ob.OBConversion()
        obc.SetInFormat("gamout")

        err_re = re.compile('^ \*\*\*')
        eng_re = re.compile('^                       TOTAL ENERGY =')

        gamin = self.write_file(mol)
        gamout = self.tempdir + "/" + self.jobname + ".out"

        chdir(self.tempdir)
        system("%s %s> %s  2> /dev/null" % (self.gamess, self.jobname, gamout))

        err_message = ""
        total_energy = 0

        with open(gamout, "r") as f:
            err_flag = False
            for l in f:
                if eng_re.match(l):
                    total_energy = float(l.split('=')[1])
                if err_re.match(l):
                    err_flag = True
                if err_flag == True:
                    err_message += l

        # エラーが出たのでstringを渡したらなおった
        # TypeError: in method 'OBConversion_ReadFile', argument 3 of type 'std::string'
        new_mol = ob.OBMol()
        s = open(gamout).read()
        obc.ReadString(new_mol, s)

        # singlepoint だと数値が入んないので対応
        new_mol.SetEnergy(total_energy)

        chdir(self.cwd)
        if not DEBUG:
            unlink(gamin)
            unlink(gamout)

        if len(err_message) > 0:
            raise GamessError(err_message)
        else:
            return new_mol

    def print_header(self):
        """ gamess header"""

        header = ""
        header += self.print_control_section()
        header += self.print_basis_section()
        header += self.print_system_section()

        return header

    def print_control_section(self):
        control_section = " $contrl "
        for k, v in self.contrl.iteritems():
            control_section += "%s=%s " % (k, v)
        control_section += " $end\n"
        return control_section

    def print_basis_section(self):
        basis_section = " $basis "
        for k, v in self.basis.iteritems():
            basis_section += "%s=%s " % (k, v)
        basis_section += "$end\n"
        return basis_section

    def print_system_section(self):
        system_section = ""
        system_section = " $SYSTEM MWORDS=30 $END\n"
        return system_section

    def print_statpt_section(self):
        statpt_section = ""
        statpt_section = " $STATPT OPTTOL=0.0001 NSTEP=20 $END\n"
        return statpt_section

    def gamess_input(self, mol):
        obc = ob.OBConversion()
        obc.SetOutFormat("gamin")
        self.contrl['mult'] = mol.GetTotalSpinMultiplicity()
        gamin_tmp = obc.WriteString(mol)
        h = self.print_header()
        return gamin_tmp.replace(" $CONTRL COORD=CART UNITS=ANGS $END\n", h[:-1])

    def write_file(self, mol):
        gamess_input_str = self.gamess_input(mol)
        gamess_input_file = self.tempdir + "/" + self.jobname + ".inp"
        with open(gamess_input_file, "w") as f:
            f.write(gamess_input_str)
        return gamess_input_file

    def __del__(self):
        if not DEBUG:
            removedirs(self.tempdir)

    def basis_type(self, basis_type):
        if basis_type in ["sto3g", "STO-3G"]:
            self.basis = {'gbasis': 'sto', 'ngauss': '3'}
        elif basis_type in ["631g", "6-31G(d)", "6-31G"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1'}
        elif basis_type in ["631gdp", "6-31G(d,p)"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1'}
        elif basis_type in["631+gdp", "6-31G+(d,p)"]:
            self.basis = {'gbasis': 'n31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1', 'diffsp': '.true', }
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
    next = obc.ReadFile(mol, "examples/ethane.mol")
    print g.gamess_input(mol)
    try:
        newmol = g.run(mol)
    except GamessError, gerr:
        print gerr.value

    print newmol.GetEnergy()
    print [(obatom.GetIdx(), obatom.GetType(), obatom.GetPartialCharge()) for obatom in ob.OBMolAtomIter(newmol)]
