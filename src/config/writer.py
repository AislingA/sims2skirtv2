# writer.py

# To automate the creation of SKIRT (.ski) configuration files. This 
# module takes a base XML template and performs targeted string replacements. 
# It manages simulation modes, source/medium file paths, spatial grid 
# boundaries, and instrument settings, ensuring the XML remains valid while 
# allowing for highly dynamic simulation setups.

# imports
import yaml
import os

def load_config(yaml_path):
    """
    Loads the replacements settings from YAML file
    """
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)

def get_default_replacements(src_file, gas_file, bounds):
    """
    Calculates the FOV and formats the geometric bounds 
    into strings with units
    
    Parameters
    ----------
    src_file : str
        Path to the source text file.
    gas_file : str
        Path to the gas text file.
    bounds : tuple
        (xmin, xmax, ymin, ymax, zmin, zmax)
        
    Returns
    -------
    dict
        A dictionary of runtime-specific replacements.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    return {
        "SOURCEFILE": src_file,
        "GASFILE": gas_file,
        "XMIN": f"{xmin} pc",
        "XMAX": f"{xmax} pc",
        "YMIN": f"{ymin} pc",
        "YMAX": f"{ymax} pc",
        "ZMIN": f"{zmin} pc",
        "ZMAX": f"{zmax} pc",
        "FOVX": f"{xmax - xmin} pc",
        "FOVY": f"{ymax - ymin} pc"
    }

def apply_yaml_replacements(template_path, output_path, runtime_replacements, yaml_path='replacements.yaml'):
    """
    Combines runtime data with YAML configuration to create the final .ski file.
    
    Parameters
    ----------
    template_path : str
        Path to the template.ski file.
    output_path : str
        Where to save the resulting .ski file.
    runtime_replacements : dict
        Dictionary from get_default_replacements() containing file paths and bounds.
    yaml_path : str
        Path to the configuration YAML file.
    """
    # loading in setting from yaml
    config_data = load_config(yaml_path)

    # runtime replacements
    final_replacements = runtime_replacements.copy()

    # flattening the yaml structure and merging
    if config_data:
        for section, settings in config_data.items():
            if isinstance(settings, dict):
                final_replacements.update(settings)

    # reading template and replacing
    with open(template_path, 'r') as f:
        content = f.read()

    for placeholder, value in final_replacements.items():
        content = content.replace(placeholder, str(value))

    with open(output_path, 'w') as f:
        f.write(content)

    print(f"Config written to: {output_path}")
    return output_path