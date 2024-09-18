#!/usr/bin/python3.8
# encoding: utf-8
# syntax:python

import sys,os,os.path,time
try:
	import _io
except ImportError:
	import io

# Add search path for yade Python-modules
# It allows to use both Yade-versions (packaged and self-compiled one).
# See LP:1254708 for more details
# (old site, fixed bug) https://bugs.launchpad.net/yade/+bug/1254708

sys.path.insert(1,'/home/missios/Downloads/yade/install/lib/x86_64-linux-gnu/yade-Unknown/py')

# get yade path (allow YADE_PREFIX to override)
prefix,suffix='/home/missios/Downloads/yade/install' if 'YADE_PREFIX' not in os.environ else os.environ['YADE_PREFIX'],'-Unknown'
# duplicate some items from yade.config here, so that we can increase verbosity when the c++ part is booting
features,version,debugbuild=' LOGGER USEFUL_ERRORS COMPLEX_MP VTK OPENMP GTS QT5 CGAL PFVFLOW PFVFLOW LINSOLV MPI TWOPHASEFLOW LS_DEM FEMLIKE GL2PS LBMFLOW THERMAL PARTIALSAT POTENTIAL_PARTICLES POTENTIAL_BLOCKS'.split(' '),'Unknown',' '

if (features[0]==''): features=features[1:]

if(True==True):
	import mpmath
	mpmath.mp.dps=(int(15)+1)

libPATH='lib/x86_64-linux-gnu'
if (libPATH[1:] == '{LIBRARY_OUTPUT_PATH}'): libPATH='lib'

## find available builds
libDir=prefix+'/'+libPATH+'/yade'+suffix
if not (os.path.exists(libDir+'/py/yade/__init__.py')):
	raise RuntimeError('Libraries are not found! ('+libDir+'/py/yade/__init__.py, /py/yade/__init__.py)')

# handle command-line options first
try:
	import argparse
except ImportError: # argparse not present, print error message
	raise RuntimeError("\n\nPlease install 'python-argparse' package.\n")
prog = os.path.basename(sys.argv[0])
par=argparse.ArgumentParser(usage='%s [options] [ simulation.xml[.bz2] | script.py [script options]]'%prog,
  prog=prog,description="Yade: open-source platform for dynamic computations. It\
  is an extensible open-source framework for discrete numerical models, focused\
  on Discrete Element Method. The computation parts are written in c++ using\
  flexible object model, allowing independent implementation of new algorithms\
  and interfaces. Python is used for rapid and concise scene construction, \
  simulation control, postprocessing and debugging.\
  Available features: %s.\
  Homepage http://www.yade-dem.org, code hosted at http://www.launchpad.net/yade."%features
  )
par.add_argument('-v','--version',help='Print version and exit.',dest='version',action='store_true')
par.add_argument('-j','--threads',help='Number of OpenMP threads to run; defaults to 1. Equivalent to setting OMP_NUM_THREADS environment variable.',dest='threads',type=int)
par.add_argument('--cores',help='Set number of OpenMP threads (as \-\-threads) and in addition set affinity of threads to the cores given. Please provide a string with comma-separated core-ids.',dest='cores',type=str)
par.add_argument('--update',help='Update deprecated class names in given script(s) using text search & replace. Changed files will be backed up with ~ suffix. Exit when done without running any simulation.',dest='updateScripts',nargs='+')
par.add_argument('--nice',help='Increase nice level (i.e. decrease priority) by given number.',dest='nice',type=int)
par.add_argument('-x',help='Exit when the script finishes',dest='exitAfter',action='store_true')
par.add_argument('-f',help=('Set "Default" filter for logging https://yade-dem.org/doc/prog.html#logging, default is -f3 for all classes.'+
                ('' if 'LOGGER' in features else " (Since this build doesn't use boost::log this option is unavailable, use cmake option -DENABLE_LOGGER=ON to enable)")),dest='verbosity',type=int)
par.add_argument('-n',help="Run without graphical interface (equivalent to unsetting the DISPLAY environment variable)",dest='nogui',action='store_true')
par.add_argument('--test',help="Run regression test suite and exit; the exists status is 0 if all tests pass, 1 if a test fails and 2 for an unspecified exception.",dest="test",action='store_true')
par.add_argument('--check',help='Run a quick series of user-defined check tests as described in /yade/scripts/checks-and-tests/checks/README and https://yade-dem.org/doc/prog.html#regression-tests',dest='check',action='store_true')
par.add_argument('--checkall',help='Run a full series of user-defined check tests as described in /yade/scripts/checks-and-tests/checks/README and https://yade-dem.org/doc/prog.html#regression-tests',dest='checkall',action='store_true')
par.add_argument('--performance',help='Starts a regular test to measure the productivity',dest='performance',action='store_true')
par.add_argument('--quickperformance',help='Starts a quick test to measure the productivity',dest='quickperformance',action='store_true')
par.add_argument('--stdperformance',help='Starts N large tests on 10000 spheres. Do not start this test while other programs are running. In such case it is unlikely that an average with low standard deviation will be found.',dest='stdperformance',type=int)
par.add_argument('script',nargs='?',default='',type=str,help=argparse.SUPPRESS)
par.add_argument('args',nargs=argparse.REMAINDER,help=argparse.SUPPRESS) # see argparse doc, par.disable_interspersed_args() from optargs module
par.add_argument('-L','--libs',help='import libraries at startup before importing yade libs. May be used when the ordering of imports matter (see fixed bug (on old site) https://bugs.launchpad.net/yade/+bug/1183402 and discussion in https://gitlab.com/yade-dev/trunk/issues/18). The option can be used multiple times, as in "yade -Llib1 -Llib2"',default=None,action='append',dest='impLibraries',type=str)
opts=par.parse_args()
args = opts.args

if opts.impLibraries:
	sys.path.append('.')
	for lib in opts.impLibraries:
		print("\n------------------------ IMPORTING ------------------------")
		# see also https://gitlab.com/yade-dev/trunk/issues/18
		print("library:", lib, "\n")
		import importlib
		importlib.import_module(str(lib))
		#__import__(lib)

if opts.version:
	print('Yade version: %s%s'%(version,debugbuild))
	sys.exit(0)

if opts.script:
	args.insert(0,opts.script) # for compatibility with userSession(), could be modified in the future

## Check if this script is exec'd or imported as part of a mpi process
# If yade is run with mpiexec prefix all processes will be in mpi_mode.
# If non-mpiexec'd yade spawns workers only workers will be in mpi_mode
opts.mpi_mode = os.getenv('OMPI_COMM_WORLD_RANK')!=None

## remove later
## python2.5 relative module imports workaround
v=sys.version_info
if v[0]==2 and v[1]<=5:
	for submodule in ('yade','gts','yade/tests'):
		sys.path.append(os.path.join(libDir,'py',submodule))

sys.path.append(os.path.join(libDir,'py'))

# run regression test suite and exit
if opts.test:
	import yade.tests
	try:
		result=yade.tests.testAll()
	except:
		# https://misc.flogisoft.com/bash/tip_colors_and_formatting
		print('\033[91m'+20*'*'+' UNEXPECTED EXCEPTION WHILE RUNNING TESTS '+20*'*'+'\033[0m')
		print(20*'*'+' '+str(sys.exc_info()[0]))
		print(20*'*'+" Please report bug at https://gitlab.com/yade-dev/trunk/issues providing the following traceback:")
		import traceback; traceback.print_exc()
		print(20*'*'+' Thank you '+20*'*')
		sys.exit(2)
	if (os.path.exists('./DirSearchYade')): os.remove('./DirSearchYade')
	if (os.path.exists('./MicroMacroAnalysis')): os.remove('./MicroMacroAnalysis')
	if result.wasSuccessful():
		print('\033[92m'+"*** ALL TESTS PASSED ***"+'\033[0m')
		sys.exit(0)
	else:
		print('\033[91m'+20*'*'+' SOME TESTS FAILED '+20*'*'+'\033[0m')
		sys.exit(1)

if not 'OPENMP' in features and (opts.cores or (opts.threads and opts.threads>1)):
	print('WARNING: compiled without OpenMP, -j/--threads/--cores have no effect.')

# OpenMP env variables must be set before loading yade libs ("import yade" below)
# changes have no effect after libgomp initializes
if opts.cores:
	if opts.threads: print('WARNING: --threads ignored, since --cores specified.')
	try:
		cores=[int(i) for i in opts.cores.split(',')]
	except ValueError:
		raise ValueError('Invalid --cores specification %s, should be a comma-separated list of non-negative integers'%opts.cores)
	opts.nthreads=len(cores)
	os.environ['GOMP_CPU_AFFINITY']=' '.join([str(c) for c in cores])
	os.environ['OMP_NUM_THREADS']=str(len(cores))
elif opts.threads: os.environ['OMP_NUM_THREADS']=str(opts.threads)
else: os.environ['OMP_NUM_THREADS']='1'

if __name__ == "__main__" and not opts.mpi_mode: # do not print this while importing yade in other python application
	sys.stdout.write('Welcome to Yade '+version+debugbuild+'\n')

# initialization and c++ plugins import
import yade
# other parts we will need soon
from yade.yexecfile import execfile
import yade.config
import yade.wrapper
import yade.log
import yade.system
import yade.runtime
# import to have quick access to function printAllVersions()
from yade.libVersions import printAllVersions

# continue option processing

if opts.updateScripts:
	yade.system.updateScripts(args)
	sys.exit(0)

if 'LOGGER' in yade.config.features:
	# by default use colors, but reading config file will override that
	yade.log.setUseColors(True)
	# set default filter level. This level is mentioned in documentation as the starting default level.
	yade.log.setLevel("Default",yade.log.WARN)
	# make sure verbosity filter is set correctly during loading the config file
	if opts.verbosity: yade.log.setLevel("Default",opts.verbosity)
	# read config file if it exists
	import os
	if(os.path.isfile( yade.log.defaultConfigFileName() )):
		yade.log.readConfigFile( yade.log.defaultConfigFileName() )
	# setting filter level from command line will override config file
	if opts.verbosity: yade.log.setLevel("Default",opts.verbosity)


# modify sys.argv in-place so that it can be handled by userSession
sys.yade_argv=sys.argv
sys.argv=yade.runtime.argv=args
yade.runtime.opts=opts

from yade import utils, pack, geom
from yade.utils import *
from yade.pack import *
from math import *

# Run the check tests listed in scripts/checks-and-tests/checks/checkList.py
if opts.check or opts.checkall:
	checksPath=libDir+'/py/yade/tests/checks'
	execfile(checksPath+'/checkList.py')

# Run performance check test
if opts.performance or opts.quickperformance or opts.stdperformance:
	checksPath=libDir+'/py/yade/tests/checks/performance'
	execfile(checksPath+'/checkPerf.py')

def checkVersions():
	# put here all yade compatibility tests with other libraries
	import yade.libVersions
	# Remember that libVersions always can return None, when library does not exist. Make sure this will not lead to a crash.
	glutVer = yade.libVersions.getVersion('freeglut')

def userSession(gui='none',qapp=None):
	# prepare nice namespace for users
	import yade.runtime
	import sys
	if __name__ != "__main__": # for importing as python module
		return

	checkVersions()
	# common ipython configuration
	# use standard terminal color codes to make F12, F11 bold, and underline the word "help", see https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences
	banner='[[ ^L clears screen, ^U kills line. '+', '.join((['\033[1mF12\033[0m controller','\033[1mF11\033[0m 3D view (press \033[4m"h"\033[0m in 3D view for \033[4mhelp\033[0m)','\033[1mF10\033[0m both','\033[1mF9\033[0m generator'] if (gui!='none') else [])+['\033[1mF8\033[0m plot'])+'. ]]'
	ipconfig=dict( # ipython options, see e.g. http://www.cv.nrao.edu/~rreid/casa/tips/ipy_user_conf.py
		prompt_in1='Yade [\#]: ',
		prompt_in2='     .\D.: ',
		prompt_out=" ->  [\#]: ",
		separate_in='',separate_out='',separate_out2='',
		#execfile=[prefix+'/lib/yade'+suffix+'/py/yade/ipython.py'],
		readline_parse_and_bind=[
			'tab: complete',
			# only with the gui; the escape codes might not work on non-linux terminals.
			]
			+(['"\e[24~": "\C-Uyade.qt.Controller();\C-M"','"\e[23~": "\C-Uyade.qt.View();\C-M"','"\e[21~": "\C-Uyade.qt.Controller(), yade.qt.View();\C-M"','"\e[20~": "\C-Uyade.qt.Generator();\C-M"'] if (gui!='none') else []) #F12,F11,F10,F9
			+['"\e[19~": "\C-Uimport yade.plot; yade.plot.plot();\C-M"', #F8
				'"\e[A": history-search-backward', '"\e[B": history-search-forward', # incremental history forward/backward
		]
	)

	# show python console
	# handle ipython version (version 10 not tested after some changes here)
	if yade.runtime.ipython_version==10:
		from IPython.Shell import IPShellEmbed
		ipshell=IPShellEmbed(banner=banner,rc_override=ipconfig)
		# save history -- a workaround for atexit handlers not being run (why?)
		# http://lists.ipython.scipy.org/pipermail/ipython-user/2008-September/005839.html
		import IPython.ipapi
		IPython.ipapi.get().IP.atexit_operations()
	elif yade.runtime.ipython_version==11:
		from IPython.frontend.terminal.embed import InteractiveShellEmbed
		# use the dict to set attributes
		for k in ipconfig: setattr(InteractiveShellEmbed,k,ipconfig[k])
		del k # leaving k as not defined (instead of ipconfig[-1]) for the user
		InteractiveShellEmbed.banner1=banner+'\n'  # called banner1 here, not banner anymore
		ipshell=InteractiveShellEmbed()
	elif yade.runtime.ipython_version>=12:
		if yade.runtime.ipython_version>=100:
			from IPython.terminal.embed import InteractiveShellEmbed
		else:
			from IPython.frontend.terminal.embed import InteractiveShellEmbed
		if yade.runtime.ipython_version>=500:
			from traitlets.config.loader import Config
			cfg = Config()
			prompt_config = cfg.TerminalInteractiveShell.prompts_class
		else:
			from IPython.config.loader import Config
			cfg = Config()
			prompt_config = cfg.PromptManager
		prompt_config.in_template = ipconfig['prompt_in1']
		prompt_config.in2_template = ipconfig['prompt_in2']
		prompt_config.out_template = ipconfig['prompt_out']
		import readline
		for k in ipconfig['readline_parse_and_bind']: readline.parse_and_bind(k)
		del k # leaving k not defined (instead of ipconfig['readline_parse_and_bind'][-1]) for the user
		InteractiveShellEmbed.config=cfg
		InteractiveShellEmbed.banner1=banner+'\n'
		while True:
			try:	#simultaneous calls to InteractiveShellEmbed() sometimes lead to race condition for /root/.ipython/* in mpi spawn (https://gitlab.com/yade-dev/trunk/issues/94)
				ipshell=InteractiveShellEmbed()
				break
			except:
				continue

		# ipython since version 5.0 stopped using libreadline, now we must make the bindings in another way
		# see: https://python-prompt-toolkit.readthedocs.io/en/1.0.14/pages/full_screen_apps.html?highlight=add_binding
		# and: https://stackoverflow.com/questions/38443907/how-does-one-set-specific-vim-bindings-in-ipython-5-0-0
		# and: https://github.com/ipython/ipython/pull/11426/commits/e08dc1b4dd30e8d2b607134693f226e8a22100ec#diff-d47c339db04f28c44a59c6662dbd53eaR228
		registry=None
		# depending on version the keybindigs registry sits in different places
		if hasattr(ipshell, 'pt_cli') and hasattr(ipshell.pt_cli, 'application') and hasattr(ipshell.pt_cli.application, 'key_bindings_registry'):
			registry = ipshell.pt_cli.application.key_bindings_registry
		if hasattr(ipshell, 'pt_app') and hasattr(ipshell.pt_app, 'key_bindings'):
			registry = ipshell.pt_app.key_bindings
		if(registry != None):
			from prompt_toolkit.keys import Keys
			@registry.add_binding(Keys.F12, eager=True)
			def invokeController_(event):
				try:
					yade.qt.Controller()
				except Exception as e:
					print(e)
			@registry.add_binding(Keys.F11, eager=True)
			def invokeView_(event):
				try:
					yade.qt.View()
				except Exception as e:
					print(e)
			@registry.add_binding(Keys.F10, eager=True)
			def invokeBoth_(event):
				try:
					yade.qt.View()
					yade.qt.Controller()
				except Exception as e:
					print(e)
			@registry.add_binding(Keys.F9, eager=True)
			def invokeGenerator_(event):
				try:
					yade.qt.Generator()
				except Exception as e:
					print(e)
			@registry.add_binding(Keys.F8, eager=True)
			def invokeGenerator_(event):
				try:
					import yade.plot
					yade.plot.plot()
				except:
					print("\nThere is no data to plot (yet), see https://yade-dem.org/doc/user.html#plotting-variables\n")

		# If IPython > 5 one need to initialize graphic gui
		if ((gui == "qt5" or gui == "qt4")and yade.runtime.ipython_version>=500):
			ipshell.enable_gui(gui)

	# start non-blocking qt4 app here; need to ask on the mailing list on how to make it functional
	## with ipython 0.11, start the even loop early (impossible with 0.10, which is thread-based)
	#if qt4 and yade.runtime.ipython_version==11:
	#	import IPython
	#	IPython.appstart_qt4(qapp)
	if len(sys.argv)>0:
		arg0=sys.argv[0]
		if (gui!='none'): yade.qt.Controller();
		if sum(bool(arg0.endswith(ext)) for ext in ('.xml','.xml.bz2','.xml.gz','.yade','.yade.gz','.yade.bz2','.bin','.bin.gz','.bin.bz2'))>0:
			if len(sys.argv)>1: raise RuntimeError('Extra arguments to saved simulation to run: '+' '.join(sys.argv[1:]))
			sys.stderr.write("Running simulation "+arg0+'\n')
		if arg0.endswith('.py'):
			def runScript(script):
				if not opts.mpi_mode: sys.stderr.write("Running script "+arg0+'\n')
				try:
					execfile(script,globals())
				except SystemExit: raise
				except: # all other exceptions
					import traceback
					traceback.print_exc()
					if yade.runtime.opts.exitAfter: sys.exit(1)
			runScript(arg0)
	## start interpreter
	if not opts.mpi_mode:
		if yade.runtime.opts.exitAfter: sys.exit(0)
		else: ipshell()
	else: # if mpi_mode start shell on master, waiting loops on workers, and exit mpi properly
		from yade import mpy as mp
		if mp.rank==0:
			if not yade.runtime.opts.exitAfter:
				mp.declareMasterInteractive()
				ipshell()
			mp.disconnect()
		else:
			# after the script execution, never exit directly. Instead, wait a command from master to do so.
			# We do this even in non-interactive mode, in that case the wait process is spawned just long enough to receive the 'disconnect' instruction from master
			mp.spawnedProcessWaitCommand()
		mp.MPI.Finalize()

## run userSession in a way corresponding to the features we use:
gui=None
yade.runtime.hasDisplay=False # this is the default initialized in the module, anyway
if 'QT5' in features: gui='qt5'
if opts.nogui or opts.mpi_mode: gui=None
if gui:
	import Xlib.display
	# PyQt4's QApplication does exit(1) if it is unable to connect to the display
	# we however want to handle this gracefully, therefore
	# we test the connection with bare xlib first, which merely raises DisplayError
	try:
		# contrary to display.Display, _BaseDisplay does not check for extensions and that avoids spurious message "Xlib.protocol.request.QueryExtension" (bug?)
		Xlib.display._BaseDisplay();
		yade.runtime.hasDisplay=True
	except:
		# usually Xlib.error.DisplayError, but there can be Xlib.error.XauthError etc as well
		# let's just pretend any exception means the display would not work
		gui=None
		print('Warning: no X rendering available (see https://bbs.archlinux.org/viewtopic.php?id=13189)')

# run remote access things, before actually starting the user session (not while imported by other python application)
if __name__ == "__main__" and not opts.mpi_mode:
	from yade import remote
	if (gui=='qt4' or gui=='qt5'):
		yade.remote.useQThread=True
		yade.remote.gui=gui
	yade.remote.runServers()


if gui==None:
	userSession()
elif gui=='qt4':
	## we already tested that DISPLAY is available and can be opened
	## otherwise Qt4 might crash at this point
	import PyQt4
	from PyQt4 import QtGui
	from PyQt4.QtCore import *
	import yade.qt # this yade.qt is different from the one that comes with qt3
	qapp=QtGui.QApplication(sys.argv)
	userSession(gui=gui,qapp=qapp)
elif gui=='qt5':
	import PyQt5
	from PyQt5 import QtGui
	from PyQt5.QtCore import *
	from PyQt5.QtWidgets import *
	import yade.qt
	qapp=QApplication(sys.argv)
	userSession(gui=gui,qapp=qapp)

if __name__ == "__main__":
	O.exitNoBacktrace()
