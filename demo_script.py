#!/usr/bin/env python3.4
import argparse
import inspect

def convert_arg_line_to_args(arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        if arg[0] == '#':
            break
        yield arg


def args_to_file(args):
    f = open(args.cfg, "w")
    for arg, value in vars(args).items():
        if not isinstance(value, (list)):
            value = [value]
        val = '{} '*len(value)
        f.write("--{} ".format(arg) + val.format(*value) + "\n")
    f.close


parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars = '@'
    )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument(
    "--POPSIZE", 
    type = int,
    help = "Set the census populations size.",
    default = 1000
    )
parser.add_argument(
    "--nloc", 
    type = int,
    help = "Set the number of unlinked loci.",
    default = 10
    )
parser.add_argument(
    "--outfile", 
    type = str,
    help = "Set the name of the output file.",
    default = "foo"
    )
parser.add_argument(
    "--cfg", 
    type = str,
    help = "Set the name of the configuration file.",
    default = "CONFIG.args"
    )
parser.add_argument(
    "--GENERATIONS", 
    type = int,
    help = "Number of generations to evolve.",
    default = 10001
    )
parser.add_argument(
    "--STEPS", 
    type = int,
    help = "Steps at which to save evolving populations",
    default = 1000
    )
parser.add_argument(
    "--sexrate", 
    type = float,
    nargs = "+",
    help = "Percentage of sexual reproduction",
    default = [0.0, 1.0]
    )
parser.add_argument(
    "--murate", 
    type = float,
    nargs = "+",
    help = "mutation rate per locus",
    default = [1e-05]*10
    )
parser.add_argument(
    "--amin", 
    type = int,
    help = "Minimum number of alleles per locus",
    default = 6
    )
parser.add_argument(
    "--amax", 
    type = int,
    help = "Maximum number of alleles per locus",
    default = 10
    )
parser.add_argument(
    "--rep", 
    type = int,
    help = "Number of replicates per population",
    default = 10
    )

args = parser.parse_args()
args_to_file(args)


# import simuPOP as sim



# def outputstat(pop):
#     'Calculate and output statistics, ignored'
#     return True


# the_args = {'initOps':[
#         sim.InitSex(),
#         sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
#         sim.InitGenotype(freq=[0.5, 0.5]),
#         sim.IdTagger(),
#         sim.PyOutput('Prevalence of disease in each age group:\n'),
# ],
# 'preOps':sim.InfoExec('age += 1'),
# 'matingScheme':sim.HeteroMating([
#     sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
#     sim.RandomMating(ops=[
#         sim.IdTagger(),
#         sim.Recombinator(intensity=1e-4)
#     ], subPops=[(0,1)]),
# ]),
# 'postOps':[
#     sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.3]),
#     sim.PyOperator(func=outputstat)
# ],
# 'gen':100,
# 'numRep':3}
# # describe this evolutionary process
# print(sim.describeEvolProcess(
#     **the_args
# ))     