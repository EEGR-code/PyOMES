"""Stock ChemistryDatabase modules.

Each module exposes a top-level :class:`~PyOMES.databases.database.ChemistryDatabase`
constant.  Extend via ``.extend()`` to add model-specific species and
reactions::

    from PyOMES.databases.anaerobic_digestion import AD_BASIC
    from PyOMES.chemistry.species import Species

    MY_DB = AD_BASIC.extend(
        species={"ButyricAcid": Species(id="ButyricAcid",
                                         atoms={"C":4,"H":8,"O":2},
                                         charge=0, MW=88.106)},
    )
"""
