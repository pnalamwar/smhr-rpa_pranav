#!/usr/bin/env python
"""
inject_synth_models.py

Written: 7-5-26 (AI assisted)
Last Updated: 7-5-26

Bulk-inject SpectralSynthesisModels into an SMHR session for heavy elements.
No abundance assumptions are made — models are created clean, ready to fit.

Usage:
    python inject_synth_models.py mystar.smh heavy_elements.json [output.smh]
"""
import sys
import json
import numpy as np

# SMH APIs: session loader, model class, LineList I/O, and utility helpers.
from smh import session as smhr_session
from smh.spectral_models import SpectralSynthesisModel
from smh.linelists import LineList
from smh import utils


def get_wavelength_range(linelist_path, padding_angstrom=0.0):
    """
    Parse a MOOG-format linelist text file and return a minimal wavelength window.
    Returns (wl_min - padding, wl_max + padding).

    - Reads the first token of each non-empty/non-comment line and converts to float.
    - Raises ValueError if no wavelengths parsed.
    - Used to set a model's `wavelength_region` metadata.
    """
    wavelengths = []
    with open(linelist_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                # skip blank lines and comments
                continue
            try:
                wl = float(line.split()[0])
                wavelengths.append(wl)
            except (ValueError, IndexError):
                # If a line can't be parsed, ignore it (robust to varied linelist formats)
                continue
    if not wavelengths:
        # Fail-fast if file contained no wavelength-looking tokens
        raise ValueError(f"No wavelengths parsed from {linelist_path}")
    return (min(wavelengths) - padding_angstrom,
            max(wavelengths) + padding_angstrom)


def inject_synth_models(session_path, config_path, output_path=None):
    # Note: the script previously generated a default output filename,
    # but here `output_path` may be provided by the caller.
    # (Original code had that line commented out; keep behavior consistent.)
    if output_path is None:
        output_path = session_path.replace('.smh', '_synth.smh')

    # Load the existing SMH session from disk. This reconstructs metadata
    # and any saved spectral models (via Session.load).
    session = smhr_session.Session.load(session_path)

    # Load the JSON config mapping linelist file -> {"elements": [...]}
    with open(config_path, 'r') as f:
        config = json.load(f)

    # If the session has a `line_list` preloaded into its metadata, use it
    # as a source of transitions to match against elements.
    session_line_list = session.metadata.get("line_list", None)

    # Build a set of linelist file paths already recorded by existing models so
    # we avoid creating duplicate models for the same external file.
    existing_linelists = set()
    for m in session.spectral_models:
        # Many models carry a `metadata['linelists']` entry listing external files.
        if hasattr(m, 'metadata') and 'linelists' in m.metadata:
            for ll in m.metadata['linelists']:
                existing_linelists.add(ll)

    new_models = []  # collect newly-created SpectralSynthesisModel instances

    # Iterate config entries: keys are linelist file paths, values contain 'elements'
    for linelist_path, info in config.items():
        elements = info.get('elements', [])  # expected like ["Ba"] or ["Ba", "La"]

        # Skip if we've already injected this linelist file previously
        if linelist_path in existing_linelists:
            print(f"Skipping {linelist_path} — already loaded.")
            continue

        # Attempt to match transitions from the session's line list (if present)
        transitions = None
        if session_line_list is not None:
            # Convert element names like "Ba" to numeric `species` representations
            # tried for neutral (" I") and singly-ionized (" II")
            target_species = []
            for elem in elements:
                for ion_suffix in (" I", " II"):
                    try:
                        sp = utils.element_to_species(elem + ion_suffix)
                        target_species.append(sp)
                    except Exception:
                        # ignore conversion failures; proceed with what we have
                        pass

            if target_species:
                # Build boolean mask for rows whose 'species' column matches any target
                mask = np.zeros(len(session_line_list), dtype=bool)
                for sp in target_species:
                    # species is encoded as a float; tolerance accounts for encoding nuance
                    mask |= (np.abs(session_line_list['species'] - sp) < 0.05)
                transitions = session_line_list[mask]
                if len(transitions) == 0:
                    # If no matches found in the session line list, fall back later
                    transitions = None

        # If we couldn't extract transitions from the session line list,
        # try reading the external linelist directly. This gives us a LineList
        # table that can be assigned to the model's `transitions`.
        if transitions is None:
            try:
                transitions = LineList.read(linelist_path)
            except Exception:
                # If reading fails, fallback to an *empty* LineList slice:
                # this leaves the model with no parsed transitions, but the
                # external file path will still be available for MOOG to use
                # when synthesizing (because we add it to reconstruct_copied_paths).
                if session_line_list is not None:
                    # zero-length slice of session line list preserves table schema
                    transitions = session_line_list[0:0]
                else:
                    transitions = LineList()  # empty LineList (may lack schema)

        # Determine a wavelength fitting window from the external linelist file
        # using the raw file's wavelengths. If this fails, skip this entry.
        try:
            wl_min, wl_max = get_wavelength_range(linelist_path)
        except Exception as e:
            print(f"ERROR reading {linelist_path}: {e}")
            continue

        # Create the SpectralSynthesisModel:
        # - `session`: parent session the model attaches to
        # - `transitions`: a LineList (possibly empty) describing transitions-- Is this a good idea to have it as empty if it fails???
        # - `elements`: list of element symbols to be fit simultaneously
        model = SpectralSynthesisModel(session, transitions, elements)

        # Metadata entries used by SMH and the synth pipeline:
        # - 'linelists': explicit external file(s) MOOG should read for detailed HFS
        model.metadata['linelists'] = [linelist_path]

        # - 'elements': store the elements list (redundant with model.elements but explicit)
        model.metadata['elements'] = elements

        # - 'wavelength_region': region used to choose an observed order and mask
        model.metadata['wavelength_region'] = [wl_min, wl_max]

        # - 'rt_abundances': explicit abundances to pass to RT; leave empty to indicate
        #   "no assumption" (user will set during fits)
        model.metadata['rt_abundances'] = {}

        new_models.append(model)
        print(f"  Created: {elements} | {linelist_path} [{wl_min:.2f}–{wl_max:.2f} Å]")

        # To ensure the external linelist is included inside the saved .smh tarball,
        # append it to session.metadata["reconstruct_copied_paths"]. Session.save()
        # will tar these files into the archive so MOOG can access them later.
        session.metadata.setdefault("reconstruct_copied_paths", [])
        if linelist_path not in session.metadata["reconstruct_copied_paths"]:
            session.metadata["reconstruct_copied_paths"].append(linelist_path)

    # If nothing new was created, print and exit.
    if not new_models:
        print("No new models to add.")
        return

    # Append the created models to the session's spectral models list.
    # Note: session.spectral_models is a convenience property returning
    # session.metadata["spectral_models"].
    session.metadata.setdefault("spectral_models", [])
    session.metadata["spectral_models"].extend(new_models)

    # Save the session. `Session.save` serializes each spectral model by
    # calling its __getstate__() and includes files listed in
    # metadata["reconstruct_copied_paths"] in the .smh tarball.
    session.save(output_path, overwrite=True)
    print(f"\nAdded {len(new_models)} synth models → {output_path}")


if __name__ == '__main__':
    # CLI entrypoint: expects at least session and config JSON
    if len(sys.argv) < 3:
        print("Usage: python inject_synth_models.py <session.smh> <config.json> [output.smh]")
        sys.exit(1)
    inject_synth_models(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3] if len(sys.argv) > 3 else None
    )