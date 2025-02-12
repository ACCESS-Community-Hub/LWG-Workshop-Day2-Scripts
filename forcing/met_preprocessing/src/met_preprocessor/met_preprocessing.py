import yaml
import xarray as xr
import met_preprocessor.standard_param as standard_param
import met_preprocessor.opt_param as opt_param
from met_preprocessor.unit_conv import UnitConversion
from met_preprocessor.utils import list_nc_files

xr.set_options(keep_attrs=True)

OUTPUT_FILE_FORMAT = "NETCDF4"
CONFIG_FILE_NAME = "config.yaml"
PARAM_MAP_FILE_NAME = "param_map.yaml"


def get_rename_param_criteria(params, param_map):
    """All input_param act as keys with the original key as value."""
    param_criteria = {}
    for param, param_info in param_map.items():
        for input_value in param_info.get("input_param", []):
            if input_value in params:
                param_criteria[input_value] = param
    return param_criteria


def get_unit_conv_params(param_map):
    """Units conversions are to be done for all params having unit in mapping."""
    return [
        param
        for param, param_attrs in param_map.items()
        if param_attrs.get("unit") is not None
    ]


with open(PARAM_MAP_FILE_NAME) as file:
    param_map = yaml.safe_load(file)


def process_dependencies(param_map):
    """Lists set of possible dependencies and their callable functions
    for a param."""
    dependencies = {}
    for param, param_info in param_map.items():
        ans = []
        for pi_calc in param_info.get("calc", []):
            func = None
            if param_info["type"] == "standard":
                func = getattr(standard_param, pi_calc["func"])
            elif param_info["type"] == "optional":
                func = getattr(opt_param, pi_calc["func"])
            else:
                raise Exception("Not yet defined for just conversion params")

            ans.append((pi_calc["deps"].split(","), func))
        dependencies[param] = ans

    return dependencies


def order_load_dep(res, dependencies, input_list):
    """
    Given a DAG, convert which order to calculate values.
    TODO: Cycle check - basic idea of converting list of lists into a set and checking visited
    nodes
    TODO: Algorithm works for DAG not maintaining ordering
    Eg: input_list = {1, 3}  dependencies = {2 : [[1, 4], [1, 3]], 4 : 1}
    Here answer should be {4 : [1] , 2 : [1, 4]}, but it will give out {2 : [1, 3], 4 : 1}
    """
    for param, dep_list in dependencies.items():
        for deps, func in dep_list:
            if set(deps).issubset(set(input_list)):
                dependencies.pop(param)
                return order_load_dep(
                    res + [(param, deps, func)], dependencies, input_list + [param]
                )
    return res


def run_met():
    """Run preprocessor for meteorological forcing dataset(s)."""

    with open(CONFIG_FILE_NAME) as file:
        config = yaml.safe_load(file)

    with open(PARAM_MAP_FILE_NAME) as file:
        param_map = yaml.safe_load(file)

    ## REVIEW: Have validator like cerberus
    file_list = []
    for dir in config.get("directories"):
        file_list += list_nc_files(dir)

    ## TODO: Have to combine everything as provenance
    ## TODO: Look more into parameter options for open_mfdataset
    ## (correctness in merging data, parallel processing)
    print(f"Loading combined dataset from {file_list}")
    dataset = xr.open_mfdataset(file_list, compat="override", coords="minimal")
    print("Loaded combined dataset")

    # Rename parameters
    param_criteria = get_rename_param_criteria(list(dataset.keys()), param_map)
    dataset = dataset.rename(param_criteria)

    # Unit conversions
    ## List of all params for unit conversions
    params = get_unit_conv_params(param_map)
    param_conv = UnitConversion(params)

    for param in params:
        if dataset.get(param) is not None:
            dataset[param] = param_conv.convert_param(
                dataset[param], param_map[param]["unit"]
            )
        else:
            print(f"Standard Stage: Skipping {param}")

    # Doing all possible calculations (Params)
    ## For strict ordering, resulting graph must be DAGs
    ## Can used memoisation + greedy approach
    pd = process_dependencies(param_map)
    dep_list = order_load_dep([], pd, list(dataset.keys()))

    for param, deps, func in dep_list:
        dep_attrs = list(map(lambda x: dataset[x], deps))
        # TODO: Try just base unit conversion
        dataset[param] = func(*dep_attrs)
        dataset[param] = dataset[param].metpy.dequantify()
        # After convert to actual units needed
        dataset[param] = param_conv.convert_param(
            dataset[param], param_map[param]["unit"]
        )

    # Only keep standard/optional variables (not including index variables)
    dataset.drop_vars(
        list(
            filter(
                lambda x: param_map.get(x, {}).get("type", "")
                not in ["standard", "optional"],
                list(dataset.keys()),
            )
        )
    )

    print("Saving dataset")

    # Combine filtered params
    dataset.to_netcdf(config["output_file"], format="NETCDF4")


if __name__ == "__main__":
    run_met()
