import xml.etree.ElementTree as ET
import os

toolkit = 'Networks/FwdInvToolkit.toolkit'

assert os.path.isfile(toolkit)

tree = ET.parse(toolkit)
root = tree.getroot()
