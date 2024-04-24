# Python code to recreate a simulated example from "Double Robust Variance Estimation"

Dependencies:
- NumPy: 1.25.2
- SciPy: 1.11.2
- Pandas: 1.4.1
- formulaic: 0.5.2
- Delicatessen: 2.2

## File Manifesto

`Applied-Example-Python.ipynb`
- Python notebook that walks through application of doubly robust estimators to the IPOP example (data not provided).

`efuncs.py`
- Contains the corresponding estimating equations for the three doubly robust estimators for `delicatessen`.

`helper.py`
- Helper functions for TMLE bounding step, and computing the IF variance from `delicatessen` output.

`Sim-Example-Python.ipynb`
- Python notebook that walks through application of doubly robust estimators to the simulated example (data provided).
