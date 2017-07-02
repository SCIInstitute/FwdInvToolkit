import os #Import OS to build file paths later
import os.path

net_dir = "Networks/"

files = os.walk(net_dir)
for dirpath, dirnames, filenames in os.walk(net_dir):
  for filename in [f for f in filenames if f.endswith(".srn5")]:
    print((os.path.join(dirpath, filename)))
    scirun_load_network(os.path.join(dirpath, filename))
    scirun_save_network(os.path.join(dirpath, filename))


scirun_force_quit()
