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

# count is first:
assert networks[0].tag == 'count'
# update actual network count here:
assert networks[0].text == '11'

# item_version next:
assert networks[1].tag == 'item_version'

# get actual network items:
network_items = networks[2:]

expected_network_count = int(networks[0].text)

assert len(network_items) == expected_network_count

# remove when done with real code
#assert False
