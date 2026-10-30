"""Pinned dataset settings. Every analysis uses this version.

VERSION is the pinned release (1822) unless the environment variable MICRONS_VERSION is set: the v1926 add-on
(scripts/run_v1926.sh) runs the complete-axon pipeline on the newer release side by side. Every data file is named
by version; outputs that are not (notebook tables, caches, R output folders) take the suffix VSFX ("" for the pinned
release, "_v1926" otherwise), so the two releases never overwrite each other.
"""
import os

DATASTACK = "minnie65_public"
PINNED = 1822  # materialized 2026-06-27 14:05 UTC; newest public version when pinned (2026-09-22)
VERSION = int(os.environ.get("MICRONS_VERSION", PINNED))
VSFX = "" if VERSION == PINNED else f"_v{VERSION}"
