# debrisfaSavageHutterFoam
This is an OpenFoam's solver developed by Alvaro Gonzalez Bilbao as part of his master degree Thesis. The solver can simulate the motion of debris flows over complex topographies including the erosion and deposition processes.

The solver and related libraries are located in the Applications folder.
In the Mesh generation folder there are two codes. The first one, wedge_generator.py, is used to create an .asc file with the geometric information of a wedge geometry, whereas the txt2obj.py code creates from any .asc file an .obj file.
The wedge case just mentioned can be run with the files located in the Run/wedge folder.
In Postprocessing there are three codes which were designed to analize the results of any simulation and generate secundary results, like flow discharge, total volume, etc.
