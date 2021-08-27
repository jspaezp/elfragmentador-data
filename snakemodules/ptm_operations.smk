REPLACEMENT_VARIABLE = "variable_mod02 = 0.0 X 0 3 -1 0 0"
REPLACEMENT_TEMPLATE = "variable_mod02 = {MASS} {AA} 0 3 -1 0 0"
REPLACEMENT_COMMAND = 'cat {{input}} | sed -e "s/{removal_variable}/{replacement_variable}/g" | tee {{output}}'

MOD_STRING_ALIASES = {
    "Kmod_Acetyl": {"AA": "K", "MASS": 42.010565},
    "Kmod_Biotinyl": {"AA": "K", "MASS": 226.077598},
    "Kmod_Butyryl": {"AA": "K", "MASS": 70.041865},
    "Kmod_Crotonyl": {"AA": "K", "MASS": 68.026215},
    "Kmod_Dimethyl": {"AA": "K", "MASS": 28.031300},
    "Kmod_Formyl": {"AA": "K", "MASS": 27.994915},
    "Kmod_Glutaryl": {"AA": "K", "MASS": 114.031694},
    "Kmod_Hydroxyisobutyryl": {"AA": "K", "MASS": 86.036779},
    "Kmod_Malonyl": {"AA": "K", "MASS": 86.000394},
    "Kmod_Methyl": {"AA": "K", "MASS": 14.015650},
    "Kmod_Propionyl": {"AA": "K", "MASS": 56.026215},
    "Kmod_Succinyl": {"AA": "K", "MASS": 100.016044},
    "Kmod_Trimethyl": {"AA": "K", "MASS": 42.046950},
    "Kmod_Ubiquitinyl": {"AA": "K", "MASS": 114.042927},
    "Pmod_Hydroxyproline": {"AA": "P", "MASS": 15.994915},
    "Rmod_Citrullin": {"AA": "R", "MASS": 0.984016},
    "Rmod_Dimethyl_asymm": {"AA": "R", "MASS": 28.031300},
    "Rmod_Dimethyl_symm": {"AA": "R", "MASS": 28.031300},
    "Rmod_Methyl": {"AA": "R", "MASS": 14.015650},
    "Ymod_Nitrotyr": {"AA": "Y", "MASS": 44.985078},
    "Ymod_Phospho": {"AA": "Y", "MASS": 79.966331},
    "Rmod_Unmod": {"AA": "R", "MASS": 0.0},
    "Kmod_Unmod": {"AA": "X", "MASS": 0.0},
    "Pmod_Unmod": {"AA": "X", "MASS": 0.0},
    "Ymod_Unmod": {"AA": "X", "MASS": 0.0},
}


rule comet_mod_parameters:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/{experiment}.params.high_high",
    run:
        ptm_params = MOD_STRING_ALIASES[wildcards.experiment]
        new_ptm = REPLACEMENT_TEMPLATE.format(
            MASS=ptm_params["MASS"], AA=ptm_params["AA"]
        )
        cmd = REPLACEMENT_COMMAND.format(
            removal_variable=REPLACEMENT_VARIABLE, replacement_variable=new_ptm
        )
        print(cmd)
        shell(cmd)
