#!/usr/bin/env python

import sys
import argparse
import json
import math
import re
from datetime import datetime

from openmm import *
from openmm.unit import *
from openmm.app import *
from openmm.app.internal.unitcell import computeLengthsAndAngles

class CRDFile(CharmmCrdFile):
    @staticmethod
    def writeFile(topology, crd, positions, file=sys.stdout):
        if isinstance(file, str):
            file = open(file, 'w')
        print("* DATE:   %s" % datetime.now().strftime("%x   %X"), file=file)
        print("*", file=file)
        print("%10d  EXT" % topology.getNumAtoms(), file=file)
        for i, atom in enumerate(topology.atoms()):
            position = [p.value_in_unit(angstroms) for p in positions[i]]
            resid = atom.residue.id + atom.residue.insertionCode
            print("%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8s%20.10f" %
                  (crd.atomno[i], crd.resno[i], crd.resname[i], crd.attype[i], position[0],
                   position[1], position[2], crd.segid[i], resid, 0.0),
                  file=file)

def addBarostat(system, conf):
    if conf['p_type'] == 'ISOTROPIC':
        barostat = MonteCarloBarostat( conf['p_ref']*bar, conf['temperature']*kelvin )
    elif conf['p_type'] == 'MEMBRANE':
        barostat = MonteCarloMembraneBarostat(
                conf['p_ref']*bar, 0.0*bar*nanometers, conf['temperature']*kelvin,
                MonteCarloMembraneBarostat.XYIsotropic,
                MonteCarloMembraneBarostat.ZFree )
    elif conf['p_type'] == 'ANISOTROPIC':
        barostat = MonteCarloAnisotropicBarostat( conf['p_ref']*bar, conf['temperature']*kelvin )
    system.addForce(barostat)

def updateNonbonded(system, top, conf):
    nbforce1 = nbforce2 = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nbforce1 = force
        elif isinstance(force, CustomNonbondedForce) and force.getNumTabulatedFunctions() == 2:
            nbforce2 = force

    exclusions = dict()
    for i in range(nbforce1.getNumExceptions()):
        p1, p2, c, s, e = nbforce1.getExceptionParameters(i)
        exclusions[p1, p2] = i

    atom_list = dict()
    for i, atom in enumerate(top.topology.atoms()):
        resname = atom.residue.name
        if not resname.startswith('HYB'): continue
        m = re.match(r'HYB(\d+)', resname)
        index = int(m[1])-1 if m != None else 0
        if not index in atom_list:
            atom_list[index] = dict()
        atom_list[index][atom.name] = i

    for index, site in enumerate(conf['sites']):
        nsub = len(site)
        for i in range(nsub):
            for j in range(i+1, nsub):
                for iatom in site[i]:
                    for jatom in site[j]:
                        atom_id1 = atom_list.get(index, atom_list[0])[iatom]
                        atom_id2 = atom_list.get(index, atom_list[0])[jatom]
                        if (atom_id1, atom_id2) in exclusions:
                            nbforce1.setExceptionParameters(exclusions[(atom_id1, atom_id2)],
                                                            atom_id1, atom_id2, 0.0, 0.1, 0.0)
                        elif (atom_id2, atom_id1) in exclusions:
                            nbforce1.setExceptionParameters(exclusions[(atom_id2, atom_id1)],
                                                            atom_id2, atom_id1, 0.0, 0.1, 0.0)
                        else:
                            nbforce1.addException(atom_id1, atom_id2, 0.0, 0.1, 0.0)
                            if nbforce2: nbforce2.addExclusion(atom_id1, atom_id2)

def get_box_type(dimensions):
    def farray_equal(array, value=None, tol=1e-3):
        if value is None:
            value = array[0]

        for elem in array:
            if abs(elem - value) > tol:
                return False
        return True

    a, b, c, alpha, beta, gamma = dimensions

    if farray_equal((a, b, c)) and farray_equal((alpha, beta, gamma), 90):
        box_type = 'cubi'
    elif farray_equal((a, b)) and not farray_equal((a, c)) and farray_equal((alpha, beta, gamma), 90):
        box_type = 'tetr'
    elif farray_equal((a, b)) and not farray_equal((a, c)) and farray_equal((alpha, beta), 90) and abs(gamma-120) < 1e-3:
        box_type = 'hexa'
    elif farray_equal((a, b, c)) and farray_equal((alpha, beta, gamma), 109.47):
        box_type = 'octa'
    else:
        raise ValueError('Failed to determine PBC box type.')

    return box_type


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='inpfile', default='input_config.inp', help='Input configuraion file')
    parser.add_argument('-f', dest='sysfile', default='prep/system.xml', help='Input OpenMM system xml file')
    parser.add_argument('-t', dest='topfile', default='prep/system.psf', help='Input topology file')
    parser.add_argument('-c', dest='crdfile', default='prep/system.crd', help='Input coordinate file')
    parser.add_argument('-s', dest='sysinfo', default='prep/sysinfo.str', help='Output system information file')
    parser.add_argument('-o', dest='outfile', default='prep/minimized.crd', help='Output coordinate file')
    parser.add_argument('--platform', help='OpenMM platform (default: CUDA)')
    args = parser.parse_args()

    # Load inputs & parameters
    print("Loading inputs")
    conf = json.load(open(args.inpfile, 'r'))
    top = CharmmPsfFile(args.topfile)
    crd = CharmmCrdFile(args.crdfile)

    system = XmlSerializer.deserialize(open(args.sysfile, 'r').read())
    if conf.get('p_type') != None: addBarostat(system, conf)
    updateNonbonded(system, top, conf)

    integrator = LangevinIntegrator(conf['temperature']*kelvin, 1.0/picosecond, 0.002*picoseconds)

    # Set platform
    DEFAULT_PLATFORMS = ('CUDA', 'OpenCL', 'CPU')
    enabled_platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    platform = args.platform
    if platform in enabled_platforms:
        platform = Platform.getPlatformByName(platform)
    else:
        for platform in DEFAULT_PLATFORMS:
            if platform in enabled_platforms:
                platform = Platform.getPlatformByName(platform)
                break
    if isinstance(platform, str):
        raise RuntimeError("Unable to find any OpenMM platform; exiting")

    print("Using platform:", platform.getName())
    prop = dict(CudaPrecision='single') if platform.getName() == 'CUDA' else dict()

    # Build simulation context
    simulation = Simulation(top.topology, system, integrator, platform, prop)

    # Assign coordinates
    box_vectors = system.getDefaultPeriodicBoxVectors()
    center = [box_vectors[i][i] / 2. for i in range(3)]
    positions = [[p[i]+center[i] for i in range(3)] for p in crd.positions]
    simulation.context.setPositions(positions)

    # Calculate initial system energy
    print("\nInitial system energy")
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    # Energy minimization
    print("\nEnergy minimization: %s steps" % conf['mini_nstep'])
    simulation.minimizeEnergy(maxIterations=conf['mini_nstep'])
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    # Production
    # Use simulated annealing
    print("\nMD run: %s steps" % conf['nstep'])
    simulation.reporters.append(
        StateDataReporter(sys.stdout, 1000, step=True, time=True, potentialEnergy=True,
                          temperature=True, progress=True, remainingTime=True, speed=True,
                          totalSteps=conf['nstep'], separator='\t')
    )
    interval = conf['temperature']/conf['nstep']
    temperature = 0.0
    for i in range(conf['nstep']):
        integrator.setTemperature(temperature*kelvin)
        simulation.step(1)
        temperature += interval

    # Write CRD file
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    unitcell = state.getPeriodicBoxVectors()
    CRDFile.writeFile(top.topology, crd, state.getPositions(), args.outfile)

    dimensions = computeLengthsAndAngles(unitcell)
    dimensions = [d * 10.0 for d in dimensions[:3]] + [d * 180.0 / math.pi for d in dimensions[3:]]

    box_type = get_box_type(dimensions)

    sysinfo = open(args.sysinfo, 'w')
    sysinfo.write('coordinates box %s %.2f %.2f %.2f %.2f %.2f %.2f\n' % (box_type, *dimensions))

