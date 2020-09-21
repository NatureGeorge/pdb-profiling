import pytest
import importlib

def test_default_config():
    from pdb_profiling import default_config
    default_config()

def import_third_party_packages() -> bool:
    return importlib.import_module("unsync.unsync") is not None

def import_self_modules():
    return all((
        importlib.import_module("pdb_profiling.fetcher.webfetch") is not None,
        importlib.import_module("pdb_profiling.utils") is not None,
        importlib.import_module("pdb_profiling.processers.pdbe.api") is not None,
        importlib.import_module("pdb_profiling.processers.pdbe.record") is not None))


def test_import():
    assert import_third_party_packages()
    assert import_self_modules()

def test_retrieve():
    from pdb_profiling.processers.pdbe.record import PDB
    assert PDB('1a01').status['status_code'] == 'REL'