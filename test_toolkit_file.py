import xml.etree.ElementTree as ET
import os

toolkit_file = 'Networks/FwdInvToolkit.toolkit'


assert os.path.isfile(toolkit_file)

tree = ET.parse(toolkit_file)
root = tree.getroot()


assert root.tag == 'boost_serialization'
assert len(root) == 1
toolkit = root[0]
assert len(toolkit) == 1
networks = toolkit[0]

for child in networks:
    print(child.tag, child.attrib)


# remove when done with real code
assert False
