# Description
Package grouping a set of pymol scripts together. 

## Requirements

# Command line
```
usage: pymtools.py [-h] -o OUTPUT [-c CONF_FILE] [--log] [--prefix PREFIX]
                   {align,visu,aver,cnstr,multialign} ...

PyMol Toolbox

positional arguments:
  {align,visu,aver,cnstr,multialign}

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory (default: None)
  -c CONF_FILE, --conf CONF_FILE
                        configuration file (default: None)
  --log                 Generate log files (default: False)
  --prefix PREFIX       output prefix (default: prot)
```

## Align
```
usage: pymtools.py align [-h] [--pretty] target mob

optional arguments:
  -h, --help  show this help message and exit

required arguments:
  target      target structure file [.pdb]
  mob         mobile structure file [.pdb]
  --pretty    Nice representation and colors
```

## Visu
```
usage: pymtools.py visu [-h] [--pretty] struct

optional arguments:
  -h, --help  show this help message and exit

required arguments:
  struct      structure file [.pdb]
  --pretty    Nice representation and colors

```
## Aver
```
usage: pymtools.py aver [-h] [--sel SEL] [--states STATES] pdb

optional arguments:
  -h, --help       show this help message and exit

required arguments:
  pdb              pdb structure with several states
  --sel SEL        Subset of atoms on which to operate (e.g., 'name CA and
                   resi 20-50')
  --states STATES  number of the first state to include in averaging and
                   number of the last state to include in averaging. A value
                   of '0' (zero) means use the last state in the object.

```
## Cnstr

```
usage: pymtools.py cnstr [-h] [-r REF] [-g GROUP] [-w] constrains pdb

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  constrains            NMR constrains in .tbl, .txt or in .csv format. For
                        .txt and .csv format, the file should follow the same
                        syntax as in examples/data folder.
  pdb                   PDB file
  -r REF, --ref REF     Reference pdb file
  -g GROUP, --group GROUP
                        Name of the field if csv file given
  -w, --writefiles      Write output txt files
```

## Multialign
```
usage: pymtools.py multialign [-h] [--pretty] target mob [mob ...]

optional arguments:
  -h, --help  show this help message and exit

required arguments:
  target      target structure file [.pdb]
  mob         Structure file(s) [.pdb]
  --pretty    Nice representation and colors
```
