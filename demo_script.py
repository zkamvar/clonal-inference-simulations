#!/usr/bin/env python3.4
import simuPOP as sim

def outputstat(pop):
    'Calculate and output statistics, ignored'
    return True


the_args = {'initOps':[
        sim.InitSex(),
        sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.IdTagger(),
        sim.PyOutput('Prevalence of disease in each age group:\n'),
],
'preOps':sim.InfoExec('age += 1'),
'matingScheme':sim.HeteroMating([
    sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
    sim.RandomMating(ops=[
        sim.IdTagger(),
        sim.Recombinator(intensity=1e-4)
    ], subPops=[(0,1)]),
]),
'postOps':[
    sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.3]),
    sim.PyOperator(func=outputstat)
],
'gen':100,
'numRep':3}
# describe this evolutionary process
print(sim.describeEvolProcess(
    **the_args
))     