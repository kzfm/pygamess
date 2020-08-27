from ruamel.yaml import YAML
import argparse, subprocess, logging
from . import gamess, email

logger = logging.getLogger(__name__)
parser = argparse.ArgumentParser("pyGamess command line interface")
parser.add_argument('input_file')
parser.add_argument('-n', '--num_cores', default=None, type=int)
parser.add_argument('-x', '--executable_num', default='01', type=str)
parser.add_argument('-s', '--rungms_suffix', default='0', type=str)
parser.add_argument('-o', '--output_file', default=None, type=str)
parser.add_argument('-g', '--pygamess_log_file', default=None, type=str)
parser.add_argument('-e', '--emails_yml', nargs="?", type=str)
parser.add_argument('-r', '--reset', action='store_true', help="Kill existent gamess processes and deletes their files in scratch folders")
args = parser.parse_args()
g = gamess.Gamess(rungms_suffix=args.rungms_suffix, executable_num=args.executable_num, num_cores=args.num_cores)
if args.output_file is None:
    args.output_file = args.input_file + ".log"
if args.pygamess_log_file is None:
	args.pygamess_log_file = f"pygamess_{args.input_file}.log"
logger.info(f"Writting pygamess log to {args.pygamess_log_file}")
rootlogger = logging.getLogger()
rootlogger.addHandler(logging.FileHandler(args.pygamess_log_file, mode="w"))
if args.reset:
	g.reset()
status = g.run_input(args.input_file, args.output_file)
lastlinesrun = subprocess.run(f"tail -n 40 {args.output_file}", shell=True, stdout=subprocess.PIPE)
yaml = YAML(typ="safe")
with open(args.emails_yml, "r") as opf:
    configuration = yaml.load(opf)

if status == 0:
	email.smtplib_email(configuration["success"]["body"].format(args.output_file, lastlinesrun.stdout.decode("UTF-8")), configuration["success"]["receivers"],
		configuration["success"]["subject"], configuration["smtp"])
else:
	email.smtplib_email(configuration["error"]["body"].format(args.output_file, status, lastlinesrun.stdout.decode("UTF-8")), configuration["error"]["receivers"],
		configuration["error"]["subject"], configuration["smtp"])
