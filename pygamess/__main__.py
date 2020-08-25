import argparse, multiprocessing
from . import gamess
parser = argparse.ArgumentParser("pyGamess command line interface")
parser.add_argument('input_file')
parser.add_argument('-n','--num_cores', default=None, type=int)
parser.add_argument('-x', '--executable_num', default='01', type=str)
parser.add_argument('-s', '--rungms_suffix', default='0', type=str)
parser.add_argument('-o', '--output_file', default=None, type=str)
parser.add_argument('-e', '--email', default=None, type=str)
args = parser.parse_args()
g = gamess.Gamess(rungms_suffix=args.rungms_suffix, executable_num=args.executable_num, num_cores=args.num_cores)
if args.output_file is None:
    output_file = args.input_file+".log"
else:
    output_file = args.output_file
status = g.run_input(args.input_file, output_file)