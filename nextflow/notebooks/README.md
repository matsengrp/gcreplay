To run notebooks manually in a docker container, using the image defined for a given process:

From the top-level parent directory of this repo:
```
docker pull quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks
docker run -it --rm -p 8888:8888 -v $(pwd):/workspace bb7e42655ce6
```
Then, within the container:
```
jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root --NotebookApp.token=''
```
Note that, at least for me, vscode is setup to use the ports, so it cannot be open at the same time.
