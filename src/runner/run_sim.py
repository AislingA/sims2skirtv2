# run_sim.py

# This is the primary entry point for the research pipeline. It coordinates
# the workflow by connecting the I/O, transformation, formatting, and 
# configuration modules. It finally triggers the SKIRT 9 execution and 
# monitors the simulation's progress and runtime performance.

# in src as root directory: python -m runner.run_sim

# imports
import time
import os
import PTS9.simulation as sm
from file_io.loader import load_snapshot, get_header_data, get_particle_data, filter_by_id
from file_io.transform import finalize_dataset
from processing.formatter import format_source_file, format_gas_file
from config.writer import get_default_replacements, apply_yaml_replacements

def run_pipeline(snapshot_path, percentage=1, verbose=True):
    """
    Executes the full end-to-end pipeline from raw snapshot to SKIRT output.

    Parameters
    ----------
    snapshot_path : str
        Path to the raw STARFORGE .hdf5 snapshot.
    percentage : float, optional
        The fraction of the box to extract (0.0 to 1.0). Default is 1.0.
    verbose : bool, optional
        If True, prints detailed statistics and progress updates.

    Returns
    -------
    sim_result : PTS9.simulation.Result or None
        The result object from the SKIRT execution, or None if the run failed.
    """
    start_time = time.time()

    if not os.path.exists(snapshot_path):
        print(f"Error: Snapshot {snapshot_path} not found.")
        return None
    
    try:
        # 1. Interfacing and transforming
        with load_snapshot(snapshot_path) as f:
            header = get_header_data(f)
            raw_gas = get_particle_data(f, 0)
            raw_sinks = get_particle_data(f, 5)

        gas_clean = filter_by_id(raw_gas)
        pt0, pt5, center, basis = finalize_dataset(header, gas_clean, raw_sinks, percentage)

        # 2. File formatting
        src_path = snapshot_path.replace('.hdf5', '_src.txt')
        gas_path = snapshot_path.replace('.hdf5', '_gas.txt')

        format_source_file(pt5, src_path, verbose=verbose)
        gas_info = format_gas_file(pt0, gas_path, verbose=verbose)
        gas_file, *bounds = gas_info

        # 3. Configurations
        template_path = 'config/template.ski'
        yaml_path = 'config/replacements.yaml'
        ski_output = snapshot_path.replace('.hdf5', '.ski')

        replacements = get_default_replacements(src_path, gas_path, bounds)
        apply_yaml_replacements(template_path, ski_output, replacements, yaml_path)
        
        # 4. Executing
        print('--- Executing SKIRT Simulation ---')
        skirt = sm.Skirt()
        sim = skirt.execute(ski_output, console='brief')

        elapsed = time.time() - start_time
        print(f"Pipeline complete. Execution time: {elapsed:.2e} seconds")
        return sim

    except Exception as e:
        print(f"Pipeline failed with error: {e}")
        return None
    
if __name__ == "__main__":
    snapshot_file = '../data/snapshot_150.hdf5'
    run_pipeline(snapshot_file, percentage=1)