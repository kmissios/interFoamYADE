# Script to simulate the conveying of metal powder for Laser Metal Deposition (LMD) using CFD-DEM
# 2024 (c) by Vasileios Angelidakis <v.angelidakis@qub.ac.uk>
# Based on the script CollapseExample.py of Eduard Puig Montella

#FIXME: Reinstate the OpenFOAM code
#FIXME: Add the functions to implement a Rosin Rammler PSD. For clumps, do I scale particles or just use the clumps? Check their sizes
#FIXME: Add a sphere factory or even better write a function to periodically add a hexaPack at the same location
#DONE: Change contact model to Hertz Mindlin + the corresponding materials + Check with Alex that we want zero adhesion!
#FIXME: Add function to create folder or empty them if already existing. I have this in another script
#FIXME: Discuss with Alex about adding a DomainLimiter. Harmonise all iterPeriod for all functions!!!
#FIXME: Decide the initial position of the particles

import os, shutil
from yadeimport import *
from yade import mpy as mp
from yade import ymport, geom

# If the directory for this particular simulation exists, delete old results
# filename='results'
# filenameVTK=filename+'/spheres/'
# if os.path.exists(filename):
# 	shutil.rmtree(filename)
# os.makedirs(filename)

# #os.makedirs(filename+'/data')
# os.makedirs(filenameVTK)

# path = "spheres"
# isExist = os.path.exists(path)
# if not isExist:
#    os.makedirs(path)
#    print("The new directory"+path+" is created!")

parallelYade=True #mpirun --allow-run-as-root -n 2 python3 scriptMPI.py , if False  python3 scriptMPI.py
numProcOF=2
# numProcOF=1


NSTEPS  = 2000000
saveVTK = 1
radius = 1
# --------------------------------------------------------------------------------------------------
# Create materials
#O.materials.append(ViscElMat(en=0.9, et=0.9, young=2.11e7, poisson=0.3, density=7980, frictionAngle=atan(0.52), label='spheres'))
#O.materials.append(ViscElMat(en=0.9, et=0.9, young=2.11e7, poisson=0.3, density=0   , frictionAngle=atan(0.52), label='walls'    ))
O.materials.append(FrictMat(young=2.1e7, poisson=0.30, density=7980, frictionAngle=atan(0.52), label='spheres'))
O.materials.append(FrictMat(young=1.0e7, poisson=0.23, density=7980, frictionAngle=atan(0.43), label='walls'    ))

# --------------------------------------------------------------------------------------------------
# Generate sphere packing
# mn = Vector3( 0.1, 2, -4) #0.3
# mx = Vector3( 0, 0, 0.0012)

# sp=pack.regularHexa(pack.inCylinder(mn,mx,0.001/6),radius=1*(66.9e-6)/2,gap=5e-5,material='spheres') #FIXME: Check spacing, and D 1mm or 1.25mm
# O.bodies.append(sp)

O.bodies.append(sphere((4,5,-3),radius))
# O.bodies.append(sphere((0,radius*4,0.001),radius))
# O.bodies.append(sphere((0,-radius*4,0.001),radius))
# O.bodies.append(sphere((radius*4,0,0.001),radius))
# O.bodies.append(sphere((-radius*4,0,0.001),radius))
# O.bodies.append(sphere((0,radius*4*2,0.001),radius))
# O.bodies.append(sphere((0,-radius*4*2,0.001),radius))
# O.bodies.append(sphere((radius*4*2,0,0.001),radius))
# O.bodies.append(sphere((-radius*4*2,0,0.001),radius))

# sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere] # FIXME: Amend this to count clumps
print ("Number of particles:", len(O.bodies))

# --------------------------------------------------------------------------------------------------
## Generate facet cylinder # FIXME: Increase its fidelity? What is considered in the CFD simulation?
#H=0.005 # FIXME: Return to reinstate the real length of 0.3m
# H=0.3
# O.bodies.append(geom.facetCylinder(center=(0,0,H/2),radius=1.25e-3/2,height=H,segmentsNumber=32,wallMask=4,material='walls'))

## Import facets
#facets=ymport.stl('pipe.stl',material='walls',fixed=True,wire=True)
#ids_list=O.bodies.append(facets)

# Set up OpenFOAM coupling parameters
fluidCoupling = FoamCoupling()
fluidCoupling.couplingModeParallel = parallelYade
fluidCoupling.isGaussianInterp = True
#use pimpleFoamYade for gaussianInterp (only in serial mode)
sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere]
#wallsLockId = O.bodies.append(box(center= (L_tank, H/2.,W/2.),extents=(0,H/2.,W/2.),fixed=True,wire=False,color = (1.,0.,0.),material='walls'))

'''The yade specific (icoFoamYade, pimpleFoamYade) OpenFOAM solver can be found in $FOAM_USER_APPBIN, (
# full path here, the scond argument, 2 is the number of FoamProcs. '''
# fluidCoupling.SetOpenFoamSolver(os.environ.get('FOAM_USER_APPBIN')+'/icoFoamYade', 2)
# it also work without path after sourcing OFoam's bashrc
fluidCoupling.SetOpenFoamSolver("interFoamYadev2312", numProcOF)


# --------------------------------------------------------------------------------------------------
# Define engines
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()], label='collider'),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		#[Ig2_Sphere_Sphere_ScGeom6D(interactionDetectionFactor=enlFactor),Ig2_Box_Sphere_ScGeom6D(interactionDetectionFactor=enlFactor), Ig2_Wall_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(en=0.9, es=0.9)], #Ip2_ViscElMat_ViscElMat_ViscElPhys
		[Law2_ScGeom_MindlinPhys_Mindlin()], #Law2_ScGeom_ViscElPhys_Basic
		label="InteractionLoop"
	),
#	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.7, timeStepUpdateInterval=100, parallelMode=parallelYade, label="ts"),
	fluidCoupling,  #to be called after timestepper
	NewtonIntegrator(gravity=(0,0,-9.81),damping=0.0,label='newton'), #(-0,-9.81,0)
	VTKRecorder(fileName='results/spheres/3d-vtk-', recorders=['spheres'], parallelMode=parallelYade, iterPeriod=saveVTK)
]

# Example usage
#O.dt=0.20*RayleighWaveTimeStep(); print(O.dt*10e6) # FIXME: Check this
# O.dt=2.13e-6; # FIXME: Check how this compares with 0.2 of the Rayleigh critical timestep
#O.dt=2e-6;

#collider.verletDist = 0.0075
mp.YADE_TIMING = False
mp.FLUID_COUPLING = True
mp.VERBOSE_OUTPUT = False
mp.USE_CPP_INTERS = True
mp.ERASE_REMOTE_MASTER = True
mp.REALLOC_FREQUENCY = 0
mp.fluidBodies = sphereIDs
#mp.commSplit = True
mp.DOMAIN_DECOMPOSITION = True

########### sedimentation ###########
# #
#while 1:
#	mp.mpirun(2000)
#	unb = unbalancedForce()
#	if unb < 0.01 :
#		newton.damping=0.1
#		for b in O.bodies:
#			if b.id ==wallsLockId:
#				mp.bodyErase(b.id)
#				print ("wallsLockId removed!")
#		break

#dataFProfile.dead=0

mp.mpirun(NSTEPS)
mp.mprint("RUN FINISH")
#fluidCoupling.killMPI()
exit()


## --------------------------------------------------------------------------------------------------
## Visualise the scene
#from yade import qt
#v=qt.View()
#v.ortho=True

#v.eyePosition=Vector3(0,-2*H,H/2)
#v.upVector=Vector3(0,0,1)
#v.viewDir=Vector3(0,1,0)

#rndr=qt.Renderer()

#O.saveTmp()
## --------------------------------------------------------------------------------------------------
## Run the script serially
#O.run(10000)
