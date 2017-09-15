import xml.etree.ElementTree as ET
import os
import glob

def verify_network_node(network_xml):
    #  this is testing the network version number of the networks in the toolkit.
    #  when the network format changes, the version number will change.
    #  If one network is upgrade or created with the new version, the rest must be upgraded too.
    #  load and save each network to upgrade it.
    #  Run the script resave_networks.py to upgrade all networks.
    assert len(network_xml) == 10 # update when version changes
    assert network_xml[0].tag == 'networkInfo'
    assert network_xml[1].tag == 'modulePositions'
    assert network_xml[2].tag == 'moduleNotes'
    assert network_xml[3].tag == 'connectionNotes'
    assert network_xml[4].tag == 'moduleTags'
    assert network_xml[5].tag == 'disabledModules'
    assert network_xml[6].tag == 'disabledConnections'
    assert network_xml[7].tag == 'moduleTagLabels'
    assert network_xml[8].tag == 'loadTagGroups'

def trim_all_whitespace(s):
    if s is not None:
        return "".join(s.split())
    return None

def elements_equal(e1, e2):
    if e1.tag != e2.tag:
        print("tags not equal")
        return False
    if trim_all_whitespace(e1.text) != trim_all_whitespace(e2.text):
        print("text not equal: ", e1.text, e2.text)
        return False
    if trim_all_whitespace(e1.tail) != trim_all_whitespace(e2.tail):
        print("tail not equal: ", e1.tail, e2.tail)
        return False
    #if e1.attrib != e2.attrib: return False
    if len(e1) != len(e2):
        print("child count not equal")
        return False
    return all(elements_equal(c1, c2) for c1, c2 in zip(e1, e2))

def compare_toolkit_network_versions(toolkit_version, network_file):
    file = 'Networks/' + network_file
    assert(os.path.isfile(file))
    actual_file_contents = ET.parse(file).getroot()[0]
    assert len(actual_file_contents) == len(toolkit_version)
    for a, b in zip(actual_file_contents, toolkit_version):
        assert elements_equal(a, b)

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

    network_files = glob.glob('Networks/*/*.srn5')
    # for filename in network_files:
    #     print(filename)

    assert len(network_items) == expected_network_count
    assert len(network_files) == expected_network_count

    for child in network_items:
        assert len(child) == 2
        assert child[0].tag == 'first'
        assert child[1].tag == 'second'
        network_file = child[0].text
        assert ('Networks/' + network_file) in network_files
        network_xml = child[1]
        verify_network_node(network_xml)
        compare_toolkit_network_versions(network_xml, network_file)
