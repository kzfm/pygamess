#!/usr/bin/env python
# -*- encoding:utf-8 -*-

import openbabel as ob
from tempfile import mkstemp, mkdtemp
from os import removedirs, unlink, system, environ, path, getcwd, chdir, system
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

        if self.debug:
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
        self.statpt = {'opttol': '0.0001', 'nstep': '20', }
        self.system = {'mwords': '30'}
        self.cis = {'nstate': '1'}


    def run(self, mol):
        self.jobname = randstr(6)
        obc = ob.OBConversion()
        obc.SetInFormat("gamout")

        err_re = re.compile('^ \*\*\*')
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
        obc.ReadString(new_mol, s)

        # singlepoint だと数値が入んないので対応
        new_mol.SetEnergy(total_energy)

        chdir(self.cwd)
        if not self.debug:
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
        if self.contrl['runtyp'] == 'optimize':
            header += self.print_statpt_section()
        if self.contrl.get('citype', None) == 'cis':
            header += self.print_cis_section()

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
        system_section = " $system "
        for k, v in self.system.iteritems():
            system_section += "%s=%s " % (k, v)
        system_section += " $end\n"
        return system_section

    def print_statpt_section(self):
        statpt_section = " $statpt "
        for k, v in self.statpt.iteritems():
            statpt_section += "%s=%s " % (k, v)
        statpt_section += " $end\n"
        return statpt_section

    def print_cis_section(self):
        cis_section = " $cis "
        for k, v in self.cis.iteritems():
            cis_section += "%s=%s " % (k, v)
        cis_section += " $end\n"
        return cis_section

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
        if not self.debug:
            removedirs(self.tempdir)

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
    print g.gamess_input(mol)
    try:
        newmol = g.run(mol)
    except GamessError, gerr:
        print gerr.value

    print newmol.GetEnergy()
    print [(obatom.GetIdx(), obatom.GetType(), obatom.GetPartialCharge()) for obatom in ob.OBMolAtomIter(newmol)]
