import matplotlib.pyplot as plt
import numpy as np

from wfn.propagation import PropagationConfig, assemble_polar, loop_angle, prepare_cases, run_local, to_grid
from wfn.receivers import ReceiversConfig


cfg = PropagationConfig.from_yaml("input_side.yaml")
prepare_cases(cfg)            # input.nml + 2D_WAPE in results/wape/t0/<h>/ and run.sh
run_local(cfg.work_path)      # runs every case one after the other
