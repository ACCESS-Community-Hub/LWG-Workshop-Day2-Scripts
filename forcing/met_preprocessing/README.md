# Met-forcings Preprocessor

## Configuration Options
### User guide

See `config.yaml` for an example configuration

1. `directories` = A list of directories from which every `.nc` file would be picked, or filename is provided, use invidual file. Recommended to use absolute paths
2. `output_file` = Output file_name for combined outputs

### Developer guide

See `param_map.yaml` for an example configuration

The top level keys are the names of all possible output preprocessor parameters. For example:

```
Qair:
  type: standard
  input_param:
    - Qair
  calc:
    - deps:
        vp,vpd,Tair
      func:
        vp_vpd_tair_sh
  unit:
    kg kg-1
```

Here, `Qair` represents specific humidity as the output.

- `type` (required) = `standard`/`optional`/`conversion` = Standard/Optional types correspond to mandatory/optional parameters used in the model. Conversion inputs are only used during calculations but not present in the final output.
- `input_param` (optional) = The input datasets would have different parameter names, they are to be renamed. Not supported for naming conflicts.
- `calc` (optional) = In case the input parameters are not provided but can be calculated with some other dependencies. Not supported for Cyclic dependencies. As of now, also not supported for ordered priority of dependencies.
- `unit` (required for optional/standard types) = Should be compatible with 

## Installation

1. On Gadi, load `analysis3` environment
2. `pip install -e .` to have editable changes for testing

## Usage

From the project root folder
`python src/met_preprocessor/met_preprocessing.py`

In case running from another directory, change the following variables to their absolute file name paths in `src/met_preprocessor/met_preprocesssing.py`

```py
CONFIG_FILE_NAME = "/path/to/config.yaml"
PARAM_MAP_FILE_NAME = "/path/to/param_map.yaml"
```

## Testing

`pytest`