import pytest

def import_third_party_packages():
    from unsync import unsync
    return True

def import_self_modules():
    from pdb_profiling.fetcher.webfetch import UnsyncFetch
    from pdb_profiling.utils import DisplayPDB
    import pdb_profiling.processers.pdbe.api
    import pdb_profiling.processers.pdbe.record
    return True


def test_mytest():
    assert import_third_party_packages()
    assert import_self_modules()