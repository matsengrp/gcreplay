# Analysis/plots of nextflow output

## Utils

### Metadata

[utils/metadata.py](utils/metadata.py) is a script to parse the metadata file and rename the mice and GCs to a more readable format.
```bash
python analysis/utils/metadata.py > analysis/output/metadata_renamed.csv
```

It can also be imported as a module.
```python
from utils.metadata import df_renamed as metadata
```

### Phenotype colormaps

[utils/phenotype_colormaps.py](utils/phenotype_colormaps.py) is a module to create colormaps for the phenotypes.
```python
import utils.phenotype_colorscales as pc

cmap = pc.affinity_dms.cmap
norm = pc.affinity_dms.norm
...
```

## Notebooks

The notebooks produce output in the `output` subdirectory.