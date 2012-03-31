#!/usr/bin/env python
# -*- encoding:utf-8 -*-

import openbabel as ob
import pybel
from tempfile import mkdtemp
from shutil import rmtree
import re
import os
import string
import socket
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

    def __init__(self, gamess_path=None, **options):
        self.tempdir = mkdtemp()
        self.debug = os.environ.get('debug', False)
        self.err_lines = 10

        if self.debug:
            print self.tempdir

        # search gamess_path
        # 1. find environ
        # 2. find path which include ddikick.x
        if gamess_path == None:
            gamess_path = os.environ.get('GAMESS_HOME', None)

        if gamess_path == None:
            try:
                gamess_path = filter(lambda f: os.path.isfile(os.path.join(f, 'ddikick.x')),
                                     [d for d in os.environ['PATH'].split(':')])[0]
            except IndexError:
                print "gamess_path not found"
                exit()

        # serch rungms script 
        rungms = None
        try:
            rungms = filter(lambda f: os.path.isfile(f),
                                 [os.path.join(d, 'rungms') for d in os.environ['PATH'].split(':')])[0]
        except IndexError:
            pass

        self.rungms = rungms
        self.gamess_path = gamess_path
        self.jobname = ''
        self.cwd = os.getcwd()
        self.contrl = {'scftyp': 'rhf', 'runtyp': 'energy'}
        self.basis = {'gbasis': 'sto', 'ngauss': '3'}
        self.statpt = {'opttol': '0.0001', 'nstep': '20', }
        self.system = {'mwords': '30'}
        self.cis = {'nstate': '1'}
        self.obc = ob.OBConversion()
        self.obc.SetInAndOutFormats("gamout", "gamin")
        
        # Todo: rewrite this
        self.contrl.update(options.get('contrl',{}))
        self.basis.update(options.get('basis',{}))
        self.statpt.update(options.get('statpt',{}))
        self.system.update(options.get('system',{}))
        self.cis.update(options.get('cis',{}))

    def parse_gamout(self, gamout):
        err_re = re.compile('^( \*\*\*|Error:)')
        eng_re = re.compile('^                       TOTAL ENERGY =')
        if self.basis['gbasis'] in ['am1', 'pm3', 'mndo']:
            eng_re = re.compile(' FINAL .+ ENERGY IS')

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

        # fixed TypeError: in method 'OBConversion_ReadFile', argument 3 of type 'std::string'
        new_mol = ob.OBMol()
        s = open(gamout).read()
        self.obc.ReadString(new_mol, s)

        # set energy when single point calculation
        new_mol.SetEnergy(total_energy)

        if len(err_message) > 0:
            raise GamessError(err_message)
        else:
            return new_mol

    def exec_rungms(self, mol):
        gamin = self.write_file(mol)
        gamout = self.tempdir + "/" + self.jobname + ".out"

        ## exec rungms
        os.chdir(self.tempdir)
        os.system("%s %s> %s  2> /dev/null" % (self.rungms, self.jobname, gamout))

        new_mol = self.parse_gamout(gamout)

        os.chdir(self.cwd)
        if not self.debug:
            os.unlink(gamin)
            os.unlink(gamout)

        return new_mol

    def run(self, mol, use_rungms=False):
        is_pybel = False
        if isinstance(mol, pybel.Molecule):
            mol = mol.OBMol
            is_pybel = True

        self.jobname = randstr(6)

        
        new_mol = self.exec_rungms(mol) if use_rungms else self.py_rungms(mol)

        # openbabel importer bug
        new_mol.SetTotalSpinMultiplicity(mol.GetTotalSpinMultiplicity())

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
        if isinstance(mol, pybel.Molecule):
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

    def basis_set(self, basis_type):
        basis_type = basis_type.upper()
        if basis_type in ["STO3G", "STO-3G"]:
            self.basis = {'gbasis': 'sto', 'ngauss': '3'}
        elif basis_type in ["321G", "3-21G"]:
            self.basis = {'gbasis': 'N21', 'ngauss': '3'}
        elif basis_type in ["631G", "6-31G"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6'}
        elif basis_type in ["6311G", "6-311G"]:
            self.basis = {'gbasis': 'N311', 'ngauss': '6'}
        elif basis_type in ["631G*", "6-31G*", "6-31G(D)", "631G(D)"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1'}
        elif basis_type in ["631G**", "6-31G**", "631GDP", "6-31G(D,P)", "631G(D,P)"]:
            self.basis = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1'}
        elif basis_type in ["631+G**", "6-31+G**", "631+GDP", "6-31+G(D,P)", "631+G(D,P)"]:
            self.basis = {'gbasis': 'n31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1', 'diffsp': '.t.', }
        elif basis_type in["AM1"]:
            self.basis = {'gbasis': 'am1'}
        elif basis_type in["PM3"]:
            self.basis = {'gbasis': 'pm3'}
        elif basis_type in["MNDO"]:
            self.basis = {'gbasis': 'mndo'}
        else:
            print "basis type not found"
        return self.basis

    def run_type(self, runtype):
        self.contrl['runtyp'] = runtype

    def scf_type(self, scftype):
        self.contrl['scftyp'] = scftype

    def py_rungms(self, mol):
        gamin = os.path.join(self.tempdir, self.jobname + ".F05")
        with open(gamin, "w") as f:
            f.write(self.input(mol))

        gamout = os.path.join(self.tempdir, self.jobname + ".out")

        gamess_path = self.gamess_path 

        ddikick = os.path.join(gamess_path, "ddikick.x")
        if not os.path.isfile(ddikick):
            raise IOError("ddikick not found")

        gamesses = [f for f in os.listdir(gamess_path) if f.startswith('gamess') and f.endswith('.x')]
        if len(gamesses) < 1:
            raise IOError("gamess.*.x not found")
        gamess = os.path.join(gamess_path, gamesses[0])
 

        hostname = socket.gethostname()

        setenv_data = [
            (" MAKEFP", "efp"), ("GAMMA", "gamma"), ("TRAJECT", "trj"),
            ("RESTART", "rst"), ("  PUNCH", "dat"), ("  INPUT", "F05"),
            (" AOINTS", "F08"), (" MOINTS", "F09"), ("DICTNRY", "F10"),
            ("DRTFILE", "F11"), ("CIVECTR", "F12"), ("CASINTS", "F13"),
            (" CIINTS", "F14"), (" WORK15", "F15"), (" WORK16", "F16"),
            ("CSFSAVE", "F17"), ("FOCKDER", "F18"), (" WORK19", "F19"),
            (" DASORT", "F20"), ("DFTINTS", "F21"), ("DFTGRID", "F22"),
            (" JKFILE", "F23"), (" ORDINT", "F24"), (" EFPIND", "F25"),
            ("SVPWRK1", "F26"), ("SVPWRK2", "F27"), ("  MLTPL", "F28"),
            (" MLTPLT", "F29"), (" DAFL30", "F30"), (" SOINTX", "F31"),
            (" SOINTY", "F32"), (" SOINTZ", "F33"), (" SORESC", "F34"),
            ("GCILIST", "F37"), ("HESSIAN", "F38"), ("QMMMTEI", "F39"),
            ("SOCCDAT", "F40"), (" AABB41", "F41"), (" BBAA42", "F42"),
            (" BBBB43", "F43"), (" MCQD50", "F50"), (" MCQD51", "F51"),
            (" MCQD52", "F52"), (" MCQD53", "F53"), (" MCQD54", "F54"),
            (" MCQD55", "F55"), (" MCQD56", "F56"), (" MCQD57", "F57"),
            (" MCQD58", "F58"), (" MCQD59", "F59"), (" MCQD60", "F60"),
            ("NMRINT1", "F61"), ("NMRINT2", "F62"), ("NMRINT3", "F63"),
            ("NMRINT4", "F64"), ("NMRINT5", "F65"), ("NMRINT6", "F66"),
            ("ELNUINT", "F67"), ("NUNUINT", "F68"), (" NUMOIN", "F69"),
            (" GMCREF", "F70"), (" GMCO2R", "F71"), (" GMCROC", "F72"),
            (" GMCOOC", "F73"), (" GMCCC0", "F74"), (" GMCHMA", "F75"),
            (" GMCEI1", "F76"), (" GMCEI2", "F77"), (" GMCEOB", "F78"),
            (" GMCEDT", "F79"), (" GMCERF", "F80"), (" GMCHCR", "F81"),
            (" GMCGJK", "F82"), (" GMCGAI", "F83"), (" GMCGEO", "F84"),
            (" GMCTE1", "F85"), (" GMCTE2", "F86"), (" GMCHEF", "F87"),
            (" GMCMOL", "F88"), (" GMCMOS", "F89"), (" GMCWGT", "F90"),
            (" GMCRM2", "F91"), (" GMCRM1", "F92"), (" GMCR00", "F93"),
            (" GMCRP1", "F94"), (" GMCRP2", "F95"), (" GMCVEF", "F96"),
            (" GMCDIN", "F97"), (" GMC2SZ", "F98"), (" GMCCCS", "F99")
            ]

        for e in setenv_data:
            os.environ[e[0].strip()] = os.path.join(self.tempdir, self.jobname + "." + e[1])

        os.environ["ERICFMT"] = os.path.join(gamess_path, "ericfmt.dat")
        os.environ["MCPPATH"] = os.path.join(gamess_path, "mcpdata")
        os.environ["EXTBAS"] = "/dev/null"
        os.environ["NUCBAS"] = "/dev/null"

        exec_string = "%s %s %s -ddi 1 1 %s -scr %s > %s" % \
            (ddikick, gamess, self.jobname, hostname, self.tempdir, gamout)
        os.system(exec_string)

        new_mol = self.parse_gamout(gamout)

        os.chdir(self.cwd)
        if not self.debug:
            os.unlink(gamin)
            os.unlink(gamout)

        return new_mol

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
