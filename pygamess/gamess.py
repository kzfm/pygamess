#!/usr/bin/env python
# -*- encoding:utf-8 -*-

from tempfile import mkdtemp
from shutil import rmtree
from random import choice
from rdkit import Chem
from . import email
import multiprocessing
import subprocess
import datetime
import logging
import socket
import re
import os

logging.basicConfig(level=logging.INFO)  # Configures logging to level INFO if the logging has not been configured
logger = logging.getLogger(__name__)


def randstr(n):
    """make a random string"""
    return ''.join(choice('abcdefghijklmnopqrstuvwxyz') for i in range(n))


class GamessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Gamess:
    """GAMESS WRAPPER"""

    def __init__(self, gamess_path=None, rungms_suffix='', executable_num='00',
        num_cores=None, reset=False, **options):
        self.tempdir = mkdtemp()
        self.debug = os.environ.get('PYGAMESS_DEBUG', False)
        self.executable_num = executable_num
        self.err_lines = 10
        if num_cores is None:
            self.num_cores = multiprocessing.cpu_count()
        else:
            self.num_cores = num_cores

        if self.debug:
            print("tmpdir", self.tempdir)

        # search gamess_path
        # 1. find environ
        # 2. find path which include ddikick.x
        if gamess_path is None:
            gamess_path = os.environ.get('GAMESS_HOME', None)

        if gamess_path is None:
            try:
                gamess_path = list(filter(lambda f: os.path.isfile(os.path.join(f, 'ddikick.x')),
                                     os.environ['PATH'].split(':')))[0]
            except IndexError:
                print("gamess_path not found")
                exit()

        #  search rungms script
        rungms = None
        possible_paths = [os.path.join(d, f'rungms{rungms_suffix}')
                          for d in os.environ['PATH'].split(':')]
        try:
            rungms = list(filter(lambda f: os.path.isfile(f), possible_paths))[0]
        except IndexError:
            pass

        self.rungms = rungms
        self.gamess_path = gamess_path
        self.jobname = ''
        #self.cwd = os.getcwd()

        #Minimal options set. Options specified in the options arguments will be
        # merged into this options set
        self.options = {
            'contrl': {'scftyp': 'rhf', 'runtyp': 'energy'},
            'basis': {'gbasis': 'sto', 'ngauss': '3'},
            'statpt': {'opttol': '0.0001', 'nstep': '20'},
            'system': {'mwords': '30'},
            'cis': {'nstate': '1'}
        }
        """self.contrl = {'scftyp': 'rhf', 'runtyp': 'energy'}
        self.basis = {'gbasis': 'sto', 'ngauss': '3'}
        self.statpt = {'opttol': '0.0001', 'nstep': '20', }
        self.system = {'mwords': '30'}
        self.cis = {'nstate': '1'}"""

        #  Todo: rewrite this
#        for configname in ('contrl', 'basis', 'statpt', 'system', 'cis'):
        for (category_name, categoptions) in options.items():
            if category_name in self.options:
                self.options[category_name].update(categoptions)
            else:
                self.options[category_name]=categoptions
            #dct = getattr(self, configname)
            #dct.update(options.get(configname, {}))

        if reset:
            self.reset()

    def discover_scratch_folders(self):
        scratches = {}
        regexpstring = re.compile(r"^\s*set\s((USER)?SCR)=(.*)")
        command = fr"grep -Pi '^\s*set\s(USER)?SCR=' {self.rungms}"
        with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True) as coprd:
            returncode = coprd.wait()
            for line in coprd.stdout:
                match = regexpstring.match(line.rstrip())
                mgrps = match.groups()
                if mgrps[0] not in scratches:
                    scratches[mgrps[0]] = os.path.expandvars(mgrps[2])
        return scratches

    def reset(self):
        subprocess.run("pkill gamess", shell=True)
        scratches = self.discover_scratch_folders()
        # delcmd = "rm -rf /scr1/lucioric/* /home/lucioric/gamess/gamess/scr/*"
        delcmd = "rm -rf {0}".format(" ".join(f"{pth}/*" for pth in scratches.values()))
        delproc = subprocess.run(delcmd, shell=True)
        delproc.check_returncode()

    def send_report_e_mail(self, num_last_lines, email_configuration, a_priori_exception=None):
        lastlinesrun = subprocess.run(f"tail -n {num_last_lines} {self.gamout}", shell=True,
                                      stdout=subprocess.PIPE, text=True)
        logger.info(f"GAMESS last {num_last_lines} lines: {lastlinesrun.stdout}")
        if a_priori_exception is not None:
            error = a_priori_exception
            which = "error"
            email_body = email_configuration[which]["body"].format(gamess=self,
                error=error, last_lines=lastlinesrun.stdout)
        else:
            error = self.completed_gamess.stderr
            if error:
                which = "fail"
                email_body = email_configuration[which]["body"].format(gamess=self,
                    error_message=error, last_lines=lastlinesrun.stdout)
            else:
                which = "success"
                email_body = email_configuration[which]["body"].format(gamess=self,
                    last_lines=lastlinesrun.stdout)
        email.smtplib_email(email_body, email_configuration[which]["receivers"],
            email_configuration[which]["subject"], email_configuration["smtp"])

    def parse_gamout(self, gamout, mol):
        err_re = re.compile('^( \*\*\*|Error:)')
        eng_re = re.compile('TOTAL ENERGY =(.*)\n')
        coord_re = re.compile('COORDINATES OF ALL ATOMS ARE (.*?)------------\n(.*?)\n\n', re.DOTALL)
        #if self.basis['gbasis'] in ['am1', 'pm3', 'mndo']:
        #    eng_re = re.compile(' FINAL .+ ENERGY IS')

        err_message = ""
        err_count = 0
        total_energy = 0

        nmol = Chem.Mol(mol)
        conf = nmol.GetConformer(0)

        out_str = open(gamout, "r").read()
        for m in eng_re.finditer(out_str):
            nmol.SetDoubleProp("total_energy", float(m.group(1).strip()))

        for m in coord_re.finditer(out_str):
            atom_idx = 0
            for l in m.group(2).split("\n"):
                coord = l.split()
                cds = [float(coord[2]), float(coord[3]), float(coord[4])]
                conf.SetAtomPosition(atom_idx, cds)
                atom_idx += 1

        #         if err_re.match(l):
        #             err_count = self.err_lines
        #             if err_count > 0:
        #                 err_message += l
        #             err_count -= 1

        if len(err_message) > 0:
            raise GamessError(err_message)
        else:
            return nmol

    def run_input(self, gamin, jobname, gamout):
        """"""
        self.jobname = jobname
        self.gamin = gamin
        self.gamout = gamout
        command = "%s %s %s %i > %s" % (self.rungms, gamin, self.executable_num,
            self.num_cores, gamout)
        self.start_time = datetime.datetime.now()
        logger.info(f"Executing {command}")
        self.completed_gamess = subprocess.run(command, shell=True)
        self.end_time = datetime.datetime.now()
        self.elapsed_time = self.end_time - self.start_time
        logger.info(f"Status code: {self.completed_gamess.returncode}")
        if self.completed_gamess.returncode != 0:
            logger.error("Gamess error: {0}".format(self.completed_gamess.stderr))
        #new_mol = self.parse_gamout(gamout)

        #$os.chdir(self.cwd)
        #if not self.debug:
        #    os.unlink(gamin)
        #    os.unlink(gamout)

    def run(self, mol, use_rungms=False):
        self.jobname = randstr(6)
        # the self properties self.gamin and self.gamout will be assigned to by the next line
        new_mol = self.exec_rungms(mol) if use_rungms else self.py_rungms(mol)
        return new_mol

    def print_header(self):
        """ gamess header"""
        header = "{}{}{}".format(self.print_section('contrl'),
                                 self.print_section('basis'),
                                 self.print_section('system'))
        if self.options['contrl']['runtyp'] == 'optimize':
            header += self.print_section('statpt')

        if self.options['contrl'].get('citype', None) == 'cis':
            header += self.print_section('cis')

        return header

    def print_section(self, pref):
        d = self.options[pref]
        section = " ${} ".format(pref)
        for k, v in d.items():
            section += "{}={} ".format(k, v)
        section += "$end\n"
        return section

    def atom_section(self, mol):
        # self.contrl['icharg'] = mol.GetFormalCharge()
        conf = mol.GetConformer(0)
        section = ""
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            section += "{:<3} {:>4}.0   {:> 15.10f} {:> 15.10f} {: 15.10f} \n".format(atom.GetSymbol(), atom.GetAtomicNum(), pos.x, pos.y, pos.z)
        return section

    def input(self, mol):
        return "{0} $DATA\n6324\nC1\n{1} $END\n".format(self.print_header(),
                                                        self.atom_section(mol))

    def write_file(self, mol):
        gamess_input_file = os.path.join(self.tempdir, "{}.inp".format(self.jobname))

        with open(gamess_input_file, "w") as f:
            f.write(self.input(mol))

        return gamess_input_file

    def __del__(self):
        if not self.debug:
            logger.info(f"deleting tempdir {self.tempdir}")
            rmtree(self.tempdir)

    def basis_set(self, basis_type):
        basis_type = basis_type.upper()
        if basis_type in ["STO3G", "STO-3G"]:
            self.options['basis'] = {'gbasis': 'sto', 'ngauss': '3'}
        elif basis_type in ["321G", "3-21G"]:
            self.options['basis'] = {'gbasis': 'N21', 'ngauss': '3'}
        elif basis_type in ["631G", "6-31G"]:
            self.options['basis'] = {'gbasis': 'N31', 'ngauss': '6'}
        elif basis_type in ["6311G", "6-311G"]:
            self.options['basis'] = {'gbasis': 'N311', 'ngauss': '6'}
        elif basis_type in ["631G*", "6-31G*", "6-31G(D)", "631G(D)"]:
            self.options['basis'] = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1'}
        elif basis_type in ["631G**", "6-31G**", "631GDP", "6-31G(D,P)", "631G(D,P)"]:
            self.options['basis'] = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1'}
        elif basis_type in ["631+G**", "6-31+G**", "631+GDP", "6-31+G(D,P)", "631+G(D,P)"]:
            self.options['basis'] = {'gbasis': 'n31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1', 'diffsp': '.t.', }
        elif basis_type in ["AM1"]:
            self.options['basis'] = {'gbasis': 'am1'}
        elif basis_type in ["PM3"]:
            self.options['basis'] = {'gbasis': 'pm3'}
        elif basis_type in ["MNDO"]:
            self.options['basis'] = {'gbasis': 'mndo'}
        else:
            logger.error("basis type not found")
        return self.options['basis']

    def run_type(self, runtype):
        self.options['contrl']['runtyp'] = runtype

    def scf_type(self, scftype):
        self.options['contrl']['scftyp'] = scftype

    def exec_rungms(self, mol):
        self.gamin = self.write_file(mol)
        self.gamout = self.tempdir + "/" + self.jobname + ".out"
        #  exec rungms
        #os.chdir(self.tempdir)
        cmd = "%s %s> %s  2> /dev/null" % (self.rungms, self.jobname, self.gamout)
        logger.info(f"Executeing exec_rungms with command {cmd}")
        self.completed_gamess = subprocess.run(cmd, shell=True, cwd=self.tempdir)
        #os.system("%s %s> %s  2> /dev/null" % (self.rungms, self.jobname, self.gamout))

        new_mol = self.parse_gamout(self.gamout, mol)

        #os.chdir(self.cwd)
        if not self.debug:
            os.unlink(self.gamin)
            os.unlink(self.gamout)

        return new_mol

    def py_rungms(self, mol):
        self.gamin = os.path.join(self.tempdir, self.jobname + ".F05")
        with open(self.gamin, "w") as f:
            f.write(self.input(mol))

        self.gamout = os.path.join(self.tempdir, self.jobname + ".out")
        gamess_path = self.gamess_path

        ddikick = os.path.join(gamess_path, "ddikick.x")

        if not os.path.isfile(ddikick):
            raise IOError("ddikick not found")

        gamesses = [f for f in os.listdir(gamess_path) if f.startswith('gamess') and f.endswith('.x')]
        if len(gamesses) < 1:
            raise IOError("gamess.*.x not found")
        gamess = os.path.join(gamess_path, gamesses[0])

        hostname = socket.gethostname()

        setenv_data = (
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
        )

        for e in setenv_data:
            os.environ[e[0].strip()] = os.path.join(self.tempdir, self.jobname + "." + e[1])

        if os.path.exists(os.path.join(gamess_path, "ericfmt.dat")):
            os.environ["ERICFMT"] = os.path.join(gamess_path, "ericfmt.dat")
        elif os.path.exists(os.path.join(gamess_path, "auxdata", "ericfmt.dat")):
            os.environ["ERICFMT"] = os.path.join(gamess_path, "auxdata", "ericfmt.dat")
        else:
            raise IOError("ericfmt.dat not found")

        os.environ["MCPPATH"] = os.path.join(gamess_path, "mcpdata")
        os.environ["EXTBAS"] = "/dev/null"
        os.environ["NUCBAS"] = "/dev/null"

        exec_string = "%s %s %s -ddi 1 1 %s -scr %s > %s" % \
            (ddikick, gamess, self.jobname, hostname, self.tempdir, self.gamout)

        logger.info(f"Executeing py_rungms with command {exec_string}")
        self.completed_gamess = subprocess.run(exec_string, shell=True)
        #os.system(exec_string)

        new_mol = self.parse_gamout(self.gamout, mol)

        #os.chdir(self.cwd)
        if not self.debug:
            os.unlink(self.gamin)
            os.unlink(self.gamout)

        return new_mol


if __name__ == '__main__':
    g = Gamess()
    mol = Chem.MolFromMolFile("examples/ethane.mol", removeHs=False)
    try:
        newmol = g.run(mol)
    except GamessError as gerr:
        logger.error(f"GamesError: {gerr.value}")

    print(newmol.GetProp("total_energy"))
