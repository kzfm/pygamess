from ruamel.yaml import YAML
import argparse, subprocess, logging
from . import gamess
import os

logger = logging.getLogger(__name__)
parser = argparse.ArgumentParser("pyGamess command line interface")
parser.add_argument('input_file')
parser.add_argument('-n', '--num_cores', default=None, type=int)
parser.add_argument('-x', '--executable_num', default='00', type=str)
parser.add_argument('-s', '--rungms_suffix', default='0', type=str)
parser.add_argument('-o', '--output_file', default=None, type=str)
parser.add_argument('-g', '--pygamess_log_file', default=None, type=str)
parser.add_argument('-e', '--emails_yml', nargs="?", type=str)
parser.add_argument('-l', '--num_last_lines', default=40, type=int)
parser.add_argument('-r', '--reset', action='store_true',
	help="Kill existent gamess processes and deletes their files in scratch folders")
args = parser.parse_args()
inputfilesplutext = os.path.splitext(args.input_file)
if args.output_file is None:
    args.output_file = inputfilesplutext[0] + ".log"
if args.pygamess_log_file is None:
	args.pygamess_log_file = f"pygamess_{inputfilesplutext[0]}.log"
yaml = YAML(typ="safe")
with open(args.emails_yml, "r") as opf:
    configuration = yaml.load(opf)
g = gamess.Gamess(rungms_suffix=args.rungms_suffix, executable_num=args.executable_num,
	num_cores=args.num_cores, reset=args.reset)
rootlogger = logging.getLogger()
logger.info(f"Writing pygamess log to {args.pygamess_log_file}")
rootlogger.addHandler(logging.FileHandler(args.pygamess_log_file, mode="w"))
try:
    g.run_input(args.input_file, inputfilesplutext[0], args.output_file)
    bdy = g.send_report_e_mail(args.num_last_lines, configuration)
except KeyboardInterrupt as exc:
    logger.error(f"Exception when runiing GAMESS: {exc}")
    bdy = g.send_report_e_mail(args.num_last_lines, configuration, a_priori_exception=exc)
