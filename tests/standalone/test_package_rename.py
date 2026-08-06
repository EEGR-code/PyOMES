def test_pyomes_is_primary_import_name():
    import PyOMES
    from PyOMES.core import Simulation

    assert PyOMES.__version__ == "0.12.5"
    assert Simulation.__name__ == "Simulation"


def test_vlsim_compatibility_import_name():
    import VLsim
    from VLsim.core import Simulation

    assert VLsim.__version__ == "0.12.5"
    assert Simulation.__name__ == "Simulation"
