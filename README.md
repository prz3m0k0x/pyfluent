## REPOSITORY DESCRITPION ##

This repo is a container for my projects regarding automatization and optimization in ANSYS Fluent.

Repository so2-optislang contains files:
- run.py
- solution.py
- geometry_script_3d_python_v.scscrpt
- reaction-mixture.scm
- response.json
The run.py is a wrapper needed for an OptiSlang solution process to call all the other files. It passes the arguments to the main script using subprocess.
Parameters for optimization are catalytic and cooling zone lenghts and SO2 mass fraction at the inlet. They are passed automatically by Optislang after initializing them in the case file with names corresponding to the ones in run.py.
The objective funcitons are written in the response.json file, which is created by solution.py script. The responses are mass out fracitons of sulfur compounds and conversion of sulfur compounds. Geometry script is called by solution.py,
as it was much simpler in this case to just track record of geometry making in the spacecalim itself in python terminal rather than writing it using ansys-geometry. Reaction mixture is a fluent file that stores all the information about the
mixuture, species and reactions. It is loaded by solution.py.

Repository so2 is in construction.
It will contain similar files as so2-optislang, but the optimization algorithm is made by me and written in Python to autonomically optimize designs rather than using OptiSlang. 

Repository photo contains come UDFs and fluent cases, along with scirpt that sets up simulation for photocatalytic conversion over a wall surface. It was my first pyfluent project,
it contains no optimization. It however utilizes pyfluent to load, compile and build UDFs to solve reactions. 
