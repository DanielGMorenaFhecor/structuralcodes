def_param_types = ["ho_la", "_mu_ndo"]
print(def_param_types)
def_param_types = [
    s.lstrip("_") for s in def_param_types
]  # remove initial "_" in parameters

print(def_param_types)
