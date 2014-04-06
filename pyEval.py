#!/usr/bin/env python

import simuPOP as sim
pop = sim.Population(1000, loci=1,
    infoFields=['mother_idx', 'father_idx'])
scheme_of_mating = sim.RandomMating(ops=[
    sim.MendelianGenoTransmitter(),
    sim.ParentsTagger(),
    ])
evals = r"'gen %d, #father %d, #mother %d\n'"
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme = scheme_of_mating,
    postOps=[
        sim.Stat(alleleFreq = 0),
        sim.PyEval(evals + " % (gen, numFather, numMother)",
            stmts="numFather = len(set(pop.indInfo('father_idx')))\n"
                "numMother = len(set(pop.indInfo('mother_idx')))",
            exposePop='pop')
    ],
    gen=3
)