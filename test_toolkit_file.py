import xml.etree.ElementTree as ET
import os

toolkit = 'Networks/FwdInvToolkit.toolkit'


assert os.path.isfile(toolkit)

tree = ET.parse(toolkit)
root = tree.getroot()


assert root.tag == 'boost_serialization'
assert len(root) == 1
for child in root[0]:
    print(child.tag, child.attrib)


# remove when done with real code
assert False
