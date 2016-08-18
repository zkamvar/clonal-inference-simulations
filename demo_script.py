#!/usr/bin/env python3.4
import simuPOP as sim
pop = sim.Population([2000], loci=[1]*50)
pop.setVirtualSplitter(sim.RangeSplitter([0, 500]))

NP = "04"
NE = "6"
NG = "03"

head, foot  = r"'", r"\n'"
popsize     = r"Pop Size: %"+NP+"d"
males       = r"Males: %"+NP+"d"
Ne          = r"Ne: %"+NE+".1f (%"+NE+".1f - %"+NE+".1f)"
generations = r"Gen: %"+NG+"d"
reps        = r"Rep: %d"

stats = " | ".join([head, generations, Ne, foot])
stateval  = " % "
stateval += "tuple("
stateval += "[gen] + "
stateval += "Ne_waples89_P1"
stateval += ")"

print(stats + stateval)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.3, 0.7]),
        sim.Stat(effectiveSize=range(50), subPops=[(0,0)],
            vars='Ne_temporal_base'),
    ],
    preOps=[
        sim.Stat(effectiveSize=range(50), subPops=[(0,0)],
            vars=['Ne_waples89_P1', 'Ne_tempoFS_P1'], step=20),
        sim.PyEval(stats + stateval, step=20)
    ],
    matingScheme=sim.RandomMating(),
    gen = 101
)