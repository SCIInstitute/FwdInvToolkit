import xml.etree.ElementTree as ET
import os

toolkit = 'Networks/FwdInvToolkit.toolkit'


assert os.path.isfile(toolkit)

tree = ET.parse(toolkit)
root = tree.getroot()

print(root.tag)
print(root.attrib)
for child in root:
    print(child.tag, child.attrib)


# remove when done with real code
assert False
