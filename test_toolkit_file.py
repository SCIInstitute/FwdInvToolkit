import xml.etree.ElementTree as ET
import os

def verify_network_node(network_xml):
    assert len(network_xml) == 9 # update when version changes
    assert network_xml[0].tag == 'networkInfo'
    assert network_xml[1].tag == 'modulePositions'
    assert network_xml[2].tag == 'moduleNotes'
    assert network_xml[3].tag == 'connectionNotes'
    assert network_xml[4].tag == 'moduleTags'
    assert network_xml[5].tag == 'disabledModules'
    assert network_xml[6].tag == 'disabledConnections'
    assert network_xml[7].tag == 'moduleTagLabels'
    assert network_xml[8].tag == 'loadTagGroups'

def compare_toolkit_network_versions(toolkit_version, network_file):
    file = 'Networks/' + network_file
    assert(os.path.isfile(file))
    assert ET.tostring(toolkit_version) == ET.tostring(ET.parse(file).getroot()[0])

def test_toolkit_file():
    toolkit_file = 'Networks/FwdInvToolkit.toolkit'

    assert os.path.isfile(toolkit_file)

    tree = ET.parse(toolkit_file)
    root = tree.getroot()

    assert root.tag == 'boost_serialization'
    assert len(root) == 1
    toolkit = root[0]
    assert len(toolkit) == 1
    networks = toolkit[0]

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

    for child in network_items:
        assert len(child) == 2
        assert child[0].tag == 'first'
        assert child[1].tag == 'second'
        network_file = child[0].text
        network_xml = child[1]
        verify_network_node(network_xml)
        compare_toolkit_network_versions(network_xml, network_file)

    # remove when done with real code
    assert False
