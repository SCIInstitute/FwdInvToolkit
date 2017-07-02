import os #Import OS to build file paths later


#Create the lat vol and edit it's parameters
lat = scirun_add_module("CreateLatVol")
scirun_set_module_state(lat, "XSize", 20)
scirun_set_module_state(lat, "YSize", 20)
scirun_set_module_state(lat, "ZSize", 20)
scirun_set_module_state(lat, "ElementSizeNormalized", 1)


#Create edit mesh bounding box to move the lat vol
edit_box = scirun_add_module("EditMeshBoundingBox")
#scirun_set_module_state(edit_box, "UseOutputCenter", 1)
scirun_set_module_state(edit_box, "OutputCenterX", 0)
scirun_set_module_state(edit_box, "OutputCenterY", 0)
scirun_set_module_state(edit_box, "OutputCenterZ", 0)

#Connect the lat vol and the edit mesh bounding box
scirun_connect_modules(lat, 0, edit_box, 0)

#Create a calculate field data module to assign values to the lat vol
calc = scirun_add_module("CalculateFieldData")

#Assign the calc field data the sphere equation
scirun_set_module_state(calc, "FunctionString", "RESULT = sqrt((X * X) + (Y * Y) + (Z * Z))")

#connect the edit mesh bounding box and the calculate field data
scirun_connect_modules(edit_box, 0, calc, 0)


#Save a copy of your network 

#Create a module to extract the isosurface from the lat vol

iso = scirun_add_module("ExtractSimpleIsosurface")

#connect the data to the isosurface
scirun_connect_modules(calc, 0, iso, 0)

#Change the isosurface module to a list and set the values
scirun_set_module_state(iso, "IsovalueChoice", "List")
scirun_set_module_state(iso, "ListOfIsovalues", "1,5,10,15,18")

#Save the network
#scirun_save_network("/YOUR/FILENAME/PATH")

#Create a showfield for the isosurface
iso_show_field = scirun_add_module("ShowField")

#Edit the showfield parameters
scirun_set_module_state(iso_show_field, "ShowFaces", 0)

#Create a colormap and rescale color map module
color_map = scirun_add_module("CreateStandardColorMap")
rescale_color_map = scirun_add_module("RescaleColorMap")

#Connect the modules
scirun_connect_modules(iso, 0, iso_show_field, 0)
scirun_connect_modules(color_map, 0, rescale_color_map, 0)
scirun_connect_modules(calc, 0, rescale_color_map, 1)
scirun_connect_modules(rescale_color_map, 0, iso_show_field, 1)

#Create the view scene
view_scene = scirun_add_module("ViewScene")

#Connect the scene 
scirun_connect_modules(iso_show_field, 0, view_scene, 0)

#Execute the network
scirun_execute_all()

## Script

'''

iso_out = scirun_add_module("WriteField")
scirun_connect_modules(iso, 0, iso_out, 0)
scirun_dump_module_state(iso_out)

scirun_set_module_state(iso, "IsovalueChoice", "Single")


iso_vals = [1,5,10,15,18]

for val in iso_vals:
	scirun_set_module_state(iso, "SingleIsoValue", val)
	fi_out = os.path.join("/Users/Gordo/Desktop/", "iso_" + str(val) + ".fld")
	scirun_set_module_state(iso_out, "Filename", fi_out)
	scirun_execute_all()

'''


