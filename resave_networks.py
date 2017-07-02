# script to resave all networks in the toolkit.
# Runs in SCIRun with the -S flag.
# usage: ./.SCIRun_test -S resave_networks.py

import os #Import OS to build file paths later
import os.path

# directory with networks
net_dir = "Networks/"

files = os.walk(net_dir)
for dirpath, dirnames, filenames in os.walk(net_dir):
  for filename in [f for f in filenames if f.endswith(".srn5")]:
    print((os.path.join(dirpath, filename)))
    scirun_load_network(os.path.join(dirpath, filename))
    scirun_save_network(os.path.join(dirpath, filename))


scirun_force_quit()
