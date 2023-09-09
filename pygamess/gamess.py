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
import platform
from .gamout_parser import gparse


logging.basicConfig(level=logging.WARNING)  # Configures logging to level INFO if the logging has not been configured
#logging.basicConfig(level=logging.DEBUG)  # Configures logging to level INFO if the logging has not been configured
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
        num_cores=None, reset=False, options={}):
        self.tempdir = mkdtemp()
        self.debug = os.environ.get('PYGAMESS_DEBUG', False)
        self.executable_num = executable_num
        self.err_lines = 10
        if num_cores is None:
            self.num_cores = multiprocessing.cpu_count()
        else:
            self.num_cores = num_cores
        logger.debug("tmpdir: {0}".format(self.tempdir))
        logger.info("#cores: {0}".format(self.num_cores))

        # search gamess_path
        # 1. find environ
        # 2. find path which include ddikick.x
        if gamess_path is None:
            gamess_path = os.environ.get('GAMESS_HOME', None)

        if gamess_path is None:
            try:
                # NOTE: ddikick is not involved in Windows installation,
                #   so Windows users rely on GAMESS_HOME environment variable.
                gamess_path = list(filter(lambda f: os.path.isfile(os.path.join(f, 'ddikick.x')),
                                     os.environ['PATH'].split(os.path.pathsep)))[0]
            except IndexError:
                print("gamess_path not found")
                exit()

        #  search rungms script
        rungms = None
        # Find rungms script inside GAMESS path
        possible_paths = [os.path.join(gamess_path, f'rungms{rungms_suffix}')]
        # Search entire PATH in case rungms isn't in GAMESS path
        possible_paths.extend([os.path.join(d, f'rungms{rungms_suffix}')
                          for d in os.environ['PATH'].split(os.path.pathsep)])
        # NOTE: Is there really a use case for searching rungms not in GAMESS path?
        #   If not, just set rungms = os.path.join(gamess_path, f'rungms{rungms_suffix}')

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
        self._options = {
            'contrl': {'scftyp': 'rhf', 'runtyp': 'energy', 'maxit': 200},
            'basis': {'gbasis': 'sto', 'ngauss': 3},
            'statpt': {'opttol': '0.0001', 'nstep': 100},
            'system': {'mwords': 300, 'memddi': 0},
            'scf': {'dirscf': '.f.'},
            'cis': {'nstate': 1}
        }

        self._options.update(options)

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
        result = gparse(gamout)
        if not result.success:
            raise GamessError(result.error_message)

        nmol = Chem.Mol(mol)
        conf = nmol.GetConformer(0)
        nmol.SetDoubleProp("total_energy", result.total_energy)
        nmol.SetDoubleProp("HOMO", result.HOMO)
        nmol.SetDoubleProp("LUMO", result.LUMO)
        nmol.SetDoubleProp("nHOMO", result.nHOMO)
        nmol.SetDoubleProp("nLUMO", result.nLUMO)
        nmol.SetDoubleProp("dipole_moment", result.dipole_moment[3])
        nmol.SetDoubleProp("dx", result.dipole_moment[0])
        nmol.SetDoubleProp("dy", result.dipole_moment[1])
        nmol.SetDoubleProp("dz", result.dipole_moment[2])
        nmol.SetProp("orbital_energies", str(result.orbital_energies))
        nmol.SetProp("program", "GAMESS")
        nmol.SetProp("basis", str(self._options['basis']))
        nmol.SetProp("method", str(self._options['contrl']))
        if hasattr(result, "uv_spectra"):
            nmol.SetProp("uv_spectra", str(result.uv_spectra))
        if hasattr(result, "isotropic_shielding"):
            nmol.SetProp("isotropic_shielding", str(result.isotropic_shielding))
        if hasattr(result, "ir_spectra"):
            nmol.SetProp("ir_spectra", str(result.ir_spectra))

        for i, cds in enumerate(result.coordinates):
            conf.SetAtomPosition(i, cds)

        for i, c in enumerate(zip(result.mulliken_charges, result.lowdin_charges, result.mulliken_populations, result.lowdin_populations)):
            atom = nmol.GetAtomWithIdx(i)
            atom.SetDoubleProp("mulliken_charge", c[0])
            atom.SetDoubleProp("lowdin_charge", c[1])
            atom.SetDoubleProp("mulliken_population", c[2])
            atom.SetDoubleProp("lowdin_population", c[3])
        Chem.CreateAtomDoublePropertyList(nmol, "mulliken_charge")
        Chem.CreateAtomDoublePropertyList(nmol, "lowdin_charge")
        Chem.CreateAtomDoublePropertyList(nmol, "mulliken_population")
        Chem.CreateAtomDoublePropertyList(nmol, "lowdin_population")

        result.mol = nmol
        return result

    def run_input(self, gamin, jobname, gamout):
        """"""
        self.jobname = jobname
        self.gamin = gamin
        self.gamout = gamout
        command = "%s %s %s %i > %s" % (self.rungms, self.gamin,
            self.executable_num, self.num_cores, self.gamout)
        self.start_time = datetime.datetime.now()
        logger.info(f"Executing {command}")

        # Run from GAMESS_PATH, lest rungms.gms config file not exist
        self.completed_gamess = subprocess.run(command, shell=True, cwd=self.gamess_path)

        self.end_time = datetime.datetime.now()
        self.elapsed_time = self.end_time - self.start_time
        logger.info(f"Status code: {self.completed_gamess.returncode}")
        if self.completed_gamess.returncode != 0:
            logger.error("Gamess error: {0}".format(self.completed_gamess.stderr))

    def run(self, mol, use_rungms=False):
        # NOTE: Default use_rungms=False incompatible with Windows version
        #   Use use_rungms=True instead.
        self.jobname = randstr(6)
        # the self properties self.gamin and self.gamout will be assigned to by the next line
        new_mol = self.exec_rungms(mol) if use_rungms else self.py_rungms(mol)
        return new_mol

    def print_header(self):
        """ gamess header"""
        # self.contrl['icharg'] = mol.GetFormalCharge()
        # TODO: cope with the charge and the multiplicity of the compound
        header = "{}{}{}{}{}".format(
            self.print_section('contrl'),
            self.print_section('pcm'),
            self.print_section('basis'),
            self.print_section('scf'),
            self.print_section('system'))

        if self._options['contrl']['runtyp'] == 'optimize':
            header += self.print_section('statpt')
            if self._options['statpt']['hssend'] == ".t.":
                header += self.print_section('cphf')

        if self._options['contrl'].get('citype', None) == 'cis':
            header += self.print_section('cis')

        if self._options['contrl'].get('tddft', None) == 'excite':
            header += self.print_section('tddft')

        if self._options['contrl'].get('tddft', None) == 'excite':
            header += self.print_section('tddft')


        return header

    def print_section(self, pref):
        # A line can only be 80 characters long. Anything after that is ignored.
        section = ""
        if pref in self._options:
            ret_flag = True
            d = self._options[pref]
            section = " ${} ".format(pref)
            num_items = len(d)
            for k, v in d.items():
                section += "{}={} ".format(k, v)
                num_items -= 1
                if ret_flag and len(section) > 60 and num_items > 0:
                    section += "\n"
                    ret_flag = False
            section += "$end\n"
        return section

    def atom_section(self, mol):
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

    def basis_sets(self, basis_type):
        basis_type = basis_type.upper()
        if basis_type in ["STO3G", "STO-3G"]:
            self._options['basis'] = {'gbasis': 'sto', 'ngauss': '3'}
        elif basis_type in ["321G", "3-21G"]:
            self._options['basis'] = {'gbasis': 'N21', 'ngauss': '3'}
        elif basis_type in ["631G", "6-31G"]:
            self._options['basis'] = {'gbasis': 'N31', 'ngauss': '6'}
        elif basis_type in ["6311G", "6-311G"]:
            self._options['basis'] = {'gbasis': 'N311', 'ngauss': '6'}
        elif basis_type in ["631G*", "6-31G*", "6-31G(D)", "631G(D)"]:
            self._options['basis'] = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1'}
        elif basis_type in ["631G**", "6-31G**", "631GDP", "6-31G(D,P)", "631G(D,P)"]:
            self._options['basis'] = {'gbasis': 'N31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1'}
        elif basis_type in ["631+G**", "6-31+G**", "631+GDP", "6-31+G(D,P)", "631+G(D,P)"]:
            self._options['basis'] = {'gbasis': 'n31', 'ngauss': '6', 'ndfunc': '1', 'npfunc': '1', 'diffsp': '.true.', }
        elif basis_type in ["AM1"]:
            self._options['basis'] = {'gbasis': 'am1'}
        elif basis_type in ["PM3"]:
            self._options['basis'] = {'gbasis': 'pm3'}
        elif basis_type in ["MNDO"]:
            self._options['basis'] = {'gbasis': 'mndo'}
        else:
            logger.error("basis type not found")
        logger.info(self._options['basis'])

    def run_type(self, runtype, hessend=False):
        self._options['contrl']['runtyp'] = runtype
        if runtype == "optimize":
            self.hessend(hessend)
        logger.debug(self._options['contrl'])

    def dft_type(self, dfttype, tddft=False, nstate=10):
        self._options['contrl']['dfttyp'] = dfttype
        if tddft:
            self._options['contrl']['tddft'] = "excite"
            if "tddft" not in self._options:
                self._options['tddft'] = {}
            self._options['tddft']['nstate'] = nstate
        else:
            self._options['contrl'].pop('tddft', None)
        logger.debug(self._options['contrl'])

    def ci_type(self, citype):
        self._options['contrl']['cityp'] = citype
        logger.debug(self._options['contrl'])

    def scf_type(self, scftype):
        self._options['contrl']['scftyp'] = scftype
        logger.debug(self._options['contrl'])

    def mul_type(self, multype):
        self._options['contrl']['multype'] = multype
        logger.debug(self._options['contrl'])

    def charge(self, charge):
        self._options['contrl']['icharg'] = charge
        logger.debug(self._options['contrl'])

    def multiplicity(self, mult):
        self._options['contrl']['mult'] = mult
        logger.debug(self._options['contrl'])

    def hessend(self, hessend):
        #  DFT ANALYTIC HESSIAN PRESENTLY HAS 5 RESTRICTIONS:
        #  $CONTRL: SCFTYP MUST BE EITHER RHF OR UHF
        #  $CONTRL: POINT GROUP SYMMETRY NOT ALLOWED, SET NOSYM=1
        #     $SCF: AO INTEGRAL DIRECT: SET DIRSCF=.TRUE.
        #    $CPHF: AO INTEGRAL DRIVEN: SET CPHF=AO
        #  AND THE FUNCTIONAL MUST NOT BE OF META-GGA TYPE.
        if hessend:
            self._options['statpt']['hssend'] = ".t."
            self._options['contrl']['nosym'] = 1
            self._options['scf']['dirscf'] = ".t."
            self._options['cphf'] = {"cphf": "AO"}
        else:
            self._options['statpt']['hssend'] = ".f."
            self._options['contrl'].pop('nosym', None)
            self._options['scf']['dirscf'] = ".f."
        logger.debug(self._options['statpt'])

    def pcm_type(self, solvent, ief=-10):
        """C-PCM is normally a better choice than IEF-PCM.  The
        iterative solvers chosen by IEF=-3 or -10 usually reproduce
        the energy of the explicit solvers IEF=3 or 10 to within
        1.0d-8 Hartrees, and will be much faster and use less
        memory for large molecules.  D-PCM should be considered
        obsolete, and choices 1 and 2 are seldom made.
        """
        if solvent == "gas":
            if "pcm" in self._options:
                self._options.pop("pcm")
        else:
            if "pcm" not in self._options:
                self._options['pcm'] = {}
            self._options['pcm']['solvnt'] = solvent
            self._options['pcm']['ief'] = ief
            logger.debug(self._options['pcm'])

    def options(self, options):
        for k, v in options.items():
            if k not in self._options:
                self._options[k] = {}
            self._options[k].update(v)
        logger.debug(self._options)

    def exec_rungms(self, mol):
        self.gamin = self.write_file(mol)
        self.gamout = os.path.join(self.tempdir, f"{self.jobname}.out")
        #  exec rungms
        #os.chdir(self.tempdir)
        cmd = "%s %s %s %i > %s" % (self.rungms, self.gamin,
            self.executable_num, self.num_cores, self.gamout)
            # NOTE: I don't think error re-direction is necessary,
            #   but if it is, cross-platform approach uses os.devnull
        logger.info(f"Executeing exec_rungms with command {cmd}")
        # Run from GAMESS_PATH, lest rungms.gms config file not exist
        self.completed_gamess = subprocess.run(cmd, shell=True, cwd=self.gamess_path)
        #os.system("%s %s> %s  2> /dev/null" % (self.rungms, self.jobname, self.gamout))

        result = self.parse_gamout(self.gamout, mol)

        #os.chdir(self.cwd)
        if not self.debug:
            os.unlink(self.gamin)
            os.unlink(self.gamout)

        return result

    def py_rungms(self, mol):
        # NOTE: Currently incompatible with Windows version, use run(.., use_rungms=True)
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

        AUXDATA = os.path.join(gamess_path, "auxdata")
        os.environ["AUXDATA"] = AUXDATA
        os.environ["DFTBPAR"] = os.path.join(AUXDATA, "DFTB")
        os.environ["ERICFMT"] = os.path.join(AUXDATA, "ericfmt.dat")
        os.environ["MCPPATH"] = os.path.join(AUXDATA, "MCP")
        os.environ["BASPATH"] = os.path.join(AUXDATA, "BASES")
        os.environ["QUANPOL"] = os.path.join(AUXDATA, "QUANPOL")
        os.environ["EXTBAS"] = "/dev/null"
        os.environ["NUCBAS"] = "/dev/null"
        os.environ["EXTCAB"] = "/dev/null"
        os.environ["POSBAS"] = "/dev/null"
        setenv_data = (
                ("MAKEFP", "efp"), ("GAMMA", "gamma"), ("TRAJECT", "trj"), ("RESTART", "rst"),("QMWAVE", "qmw"),("MDDIP", "dip"),
                ("OPTHES1", "hs1"),("OPTHES2", "hs2"), ("XYZ", "xyz"),("INPUT", "F05"),("PUNCH", "dat"),("AOINTS", "F08"),
                ("MOINTS", "F09"),("DICTNRY", "F10"),("DRTFILE", "F11"),("CIVECTR", "F12"),("CASINTS", "F13"),("CIINTS", "F14"),
                ("WORK15", "F15"),("WORK16", "F16"),("CSFSAVE", "F17"),("FOCKDER", "F18"),("WORK19", "F19"),("DASORT", "F20"),
                ("DIABDAT", "F21"),("DFTINTS", "F21"),("DFTGRID", "F22"),("JKFILE", "F23"),("ORDINT", "F24"),("EFPIND", "F25"),
                ("PCMDATA", "F26"),("PCMINTS", "F27"),("SVPWRK1", "F26"),("SVPWRK2", "F27"),("COSCAV ", "F26"),("COSDATA", "cosmo"),
                ("COSPOT ", "pot"),("MLTPL ", "F28"),("MLTPLT", "F29"),("DAFL30", "F30"),("SOINTX", "F31"),("SOINTY", "F32"),
                ("SOINTZ", "F33"),("SORESC", "F34"),("GCILIST", "F37"),("HESSIAN", "F38"),("QMMMTEI", "F39"),("SOCCDAT", "F40"),
                ("AABB41", "F41"),("BBAA42", "F42"),("BBBB43", "F43"),("REMD  ", "F44"),("UNV   ", "F45"),("UNPV  ", "F46"),
                ("MCQD50", "F50"),("MCQD51", "F51"),("MCQD52", "F52"),("MCQD53", "F53"),("MCQD54", "F54"),("MCQD55", "F55"),
                ("MCQD56", "F56"),("MCQD57", "F57"),("MCQD58", "F58"),("MCQD59", "F59"),("MCQD60", "F60"),("MCQD61", "F61"),
                ("MCQD62", "F62"),("MCQD63", "F63"),("MCQD64", "F64"),("NMRINT1", "F61"),("NMRINT2", "F62"),("NMRINT3", "F63"),
                ("NMRINT4", "F64"),("NMRINT5", "F65"),("NMRINT6", "F66"),("DCPHFH2", "F67"),("DCPHF21", "F68"),("ELNUINT", "F67"),
                ("NUNUINT", "F68"),("GVVPT", "F69"),("NUMOIN ", "F69"),("NUMOCAS", "F70"),("NUELMO ", "F71"),("NUELCAS", "F72"),
                ("RIVMAT ", "F51"),("RIT2A  ", "F52"),("RIT3A  ", "F53"),("RIT2B  ", "F54"),("RIT3B  ", "F55"),("DEN2P1", "F70"),
                ("DEN2P2", "F71"),("DEN2P3", "F72"),("DEN2P4", "F73"),("DEN2NM", "F74"),("DEN2OPT", "F75"),("GMCREF", "F70"),
                ("GMCO2R", "F71"),("GMCROC", "F72"),("GMCOOC", "F73"),("GMCCC0", "F74"),("GMCHMA", "F75"),("GMCEI1", "F76"),
                ("GMCEI2", "F77"),("GMCEOB", "F78"),("GMCEDT", "F79"),("GMCERF", "F80"),("GMCHCR", "F81"),("GMCGJK", "F82"),
                ("GMCGAI", "F83"),("GMCGEO", "F84"),("GMCTE1", "F85"),("GMCTE2", "F86"),("GMCHEF", "F87"),("GMCMOL", "F88"),
                ("GMCMOS", "F89"),("GMCWGT", "F90"),("GMCRM2", "F91"),("GMCRM1", "F92"),("GMCR00", "F93"),("GMCRP1", "F94"),
                ("GMCRP2", "F95"),("GMCVEF", "F96"),("GMCDIN", "F97"),("GMC2SZ", "F98"),("GMCCCS", "F99"),("CCREST", "F70"),
                ("CCDIIS", "F71"),("CCINTS", "F72"),("CCT1AMP", "F73"),("CCT2AMP", "F74"),("CCT3AMP", "F75"),("CCVM", "F76"),
                ("CCVE", "F77"),("CCQUADS", "F78"),("QUADSVO", "F79"),("EOMSTAR", "F80"),("EOMVEC1", "F81"),("EOMVEC2", "F82"),
                ("EOMHC1", "F83"),("EOMHC2", "F84"),("EOMHHHH", "F85"),("EOMPPPP", "F86"),("EOMRAMP", "F87"),("EOMRTMP", "F88"),
                ("EOMDG12", "F89"),("MMPP", "F90"),("MMHPP", "F91"),("MMCIVEC", "F92"),("MMCIVC1", "F93"),("MMCIITR", "F94"),
                ("EOMVL1", "F95"),("EOMVL2", "F96"),("EOMLVEC", "F97"),("EOMHL1", "F98"),("EOMHL2", "F99"),("CCVVVV", "F80"),
                ("AMPROCC", "F70"),("ITOPNCC", "F71"),("FOCKMTX", "F72"),("LAMB23", "F73"),("VHHAA", "F74"),("VHHBB", "F75"),
                ("VHHAB", "F76"),("VMAA", "F77"),("VMBB", "F78"),("VMAB", "F79"),("VMBA", "F80"),("VHPRAA", "F81"),
                ("VHPRBB", "F82"),("VHPRAB", "F83"),("VHPLAA", "F84"),("VHPLBB", "F85"),("VHPLAB", "F86"),("VHPLBA", "F87"),
                ("VEAA", "F88"),("VEBB", "F89"),("VEAB", "F90"),("VEBA", "F91"),("VPPPP", "F92"),("INTERM1", "F93"),
                ("INTERM2", "F94"),("INTERM3", "F95"),("ITSPACE", "F96"),("INSTART", "F97"),("ITSPC3", "F98"))

        for e in setenv_data:
            os.environ[e[0].strip()] = os.path.join(self.tempdir, self.jobname + "." + e[1])

        #if os.path.exists(os.path.join(gamess_path, "ericfmt.dat")):
        #    os.environ["ERICFMT"] = os.path.join(gamess_path, "ericfmt.dat")
        #elif os.path.exists(os.path.join(gamess_path, "auxdata", "ericfmt.dat")):
        #    os.environ["ERICFMT"] = os.path.join(gamess_path, "auxdata", "ericfmt.dat")
        #else:
        #    raise IOError("ericfmt.dat not found")

        if "pcm" in self._options:
            os.environ["PCMDATA"] = os.path.join(self.tempdir, "{}.F26".format(self.jobname))
            os.environ["PCMINTS"] = os.path.join(self.tempdir, "{}.F27".format(self.jobname))


        # exec string
        # See rungms examples
        if platform.system() == "Darwin":
            #       4. How to run in a single computer, namely the "localhost", so
            #          this computer needn't have a proper Internet name.
            #          This example also presumes SysV was deliberately *not* chosen
            #          when DDI was compiled, so that host names have to be repeated,
            #          instead of using the simpler localhost:cpus=$NCPU form.
            hostname = "localhost"
            if self.num_cores > 1:
                hostname = " ".join([hostname] * self.num_cores)

            exec_string = "{0} {1} {2} -ddi {3} {3} {4} -scr {5} > {6}".format(ddikick, gamess, self.jobname,
                    self.num_cores, hostname, self.tempdir, self.gamout)
        else:
            if self.num_cores == 1:
                exec_string = "{0} {1} {2} -ddi 1 {3} {4} -scr {5} > {6}".format(ddikick, gamess, self.jobname,
                    self.num_cores, hostname, self.tempdir, self.gamout)
            else:
                exec_string = "{0} {1} {2} -ddi 1 {3} {4}:cpus={3} -scr {5} > {6}".format(ddikick, gamess, self.jobname,
                    self.num_cores, hostname, self.tempdir, self.gamout)

        logger.info(f"Executeing py_rungms with command {exec_string}")
        self.completed_gamess = subprocess.run(exec_string, shell=True)
        #os.system(exec_string)

        result = self.parse_gamout(self.gamout, mol)


        if not self.debug:
            os.unlink(self.gamin)
            os.unlink(self.gamout)

        return result


if __name__ == '__main__':
    g = Gamess()
    mol = Chem.MolFromMolFile("examples/ethane.mol", removeHs=False)
    r = g.run(mol)
    print(r.mol.GetProp("total_energy"))
