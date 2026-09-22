# wiring-worth

A secondary analysis of the MICrONS dataset (mouse visual cortex: an EM connectome co-registered with calcium imaging), accessed through `caveclient` (datastack `minnie65_public`).

**Questions**
1. **December:** how does the like-to-like wiring rule change with reconstruction completeness (proofreading level)?
2. **Feb–Mar:** how well does a connectome-constrained model of activity do, compared against shuffled-wiring and data-driven controls?

## Reproducing the setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name microns
```

CAVE access (one time):
1. Accept the MICrONS terms of service with your Google account:
   https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept
2. Run `.venv/bin/python scripts/save_token.py` and paste the token when prompted.

All analyses use the materialization version pinned in `config.py`. Run scripts from the repo root, for example `.venv/bin/python notebooks/01_setup_and_explore.py`. Downloaded tables go to `data/`, which is not committed.
