from pqdm.processes import pqdm
import sys
sys.path.append('.')

def pqdm_map(f, data, n_jobs=10):
  """
  f(): function that processes the item: (idx, item) -> result
  data: list of items to process
  n_jobs: number of parallel jobs to run
  """
  return pqdm(data, f, n_jobs=n_jobs, total=len(data), smoothing=0)
