# script to resave all networks in the toolkit.
# Runs in SCIRun with the -S flag.
# usage: ./SCIRun_test -S resave_networks.py

import os #Import OS to build file paths later

# directory with networks
net_dir = "Networks/"

#print(os.path.dirname(os.environ['SCIRUNDATADIR']))

cwd= os.getcwd()
filepath = os.path.dirname(__file__)
#print(os.path.join(cwd,filepath))

if not os.path.samefile(os.environ['SCIRUNDATADIR'],os.path.join(cwd,filepath)):
  print("SCIRUNDATADIR does not seem to be correctly set.  Please set SCIRUNDATADIR to the FwdInvToolkit directory (location of this file) in SCIRun before running this script")
  scirun_force_quit()

files = os.walk(net_dir)
for dirpath, dirnames, filenames in os.walk(net_dir):
  for filename in [f for f in filenames if f.endswith(".srn5")]:
    print((os.path.join(dirpath, filename)))
    scirun_load_network(os.path.join(dirpath, filename))
    scirun_save_network(os.path.join(dirpath, filename))


scirun_force_quit()
