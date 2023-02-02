# scanpy-wdl-project
A project for wrapping some basic pre-analysis tools from Scanpy in WDL


Usage - 

Please modify the input.json file so that scanpy.channels matches the absolute path to the input files for your system.

Ensure that the Cromwell jar, scanpy.wdl, inputs.json, and default.conf files are in the working directory (alternatively, you can modify the command below with the absolute paths).

Additionally, ensure that Docker Desktop is currently running.

Afterwards,

java -Dconfig.file=default.conf -jar cromwell.jar run scanpy.wdl -i inputs.json

will run the workflow in the present directory. umap PNG files can be found in the embedding task directory. gene rank PNGs and the final h5ad files can be found in the rank_genes task directory (under cromwell-executions).

Note: the Python scripts included with the WDL are already added to the Docker image. You can access the image via atex91/scanpy-project or you can build the image from the Dockerfile provided.


The workflow itself calls these python scripts in succession, feeding the output from one task into the input from the next. Each script converts a set of scanpy operations into a simple command line tool.
