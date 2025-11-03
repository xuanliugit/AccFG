# Installation

AccFG can be installed from PyPI or directly from the source repository. The recommended path is via `pip`, but advanced users can bootstrap a local development copy.

## 1. Install from PyPI (recommended)

```bash
pip install accfg
```

This installs the latest published release together with its Python dependencies. AccFG targets Python 3.10; create a dedicated virtual environment if you need to isolate the package.

## 2. Install from the GitHub repository

Use this route when you need the latest commits or intend to contribute.

```bash
git clone https://github.com/xuanliugit/AccFG.git
cd AccFG

conda create --name accfg python=3.10
conda activate accfg
pip install -r requirements.txt  # or: pip install -e .
```

The editable install (`pip install -e .`) keeps the package linked to your working tree so changes are picked up immediately.

## 3. Verify the installation

Once installed, run a quick functional group extraction to confirm everything works:

```bash
python -c "from accfg import AccFG; print(AccFG().run('CCO'))"
```

If the command prints a dictionary of functional group matches without errors, the installation succeeded.
