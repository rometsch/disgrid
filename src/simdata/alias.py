def check_alias(a):
    # Allow strings or dicts containing only strings
    if isinstance(a, str):
        return True
    if isinstance(a, dict):
        if all([isinstance(s, str) for k,s in a.items()]):
            return True
    return False

class Alias:
    def __init__(self, alias_dict=None):
        if alias_dict is None:
            self.ad = dict()
        else:
            self.ad = alias_dict

    def __call__(self, name, variable_type = None):
        # return an alias if one is known
        if name in self.ad:
            # check for an alias for {name}
            alias = self.ad[name]
            if isinstance(alias, str):
                name = alias
            elif isinstance(alias, dict):
                # assume that there are separate aliases
                # for each variable type
                try:
                    name = alias[variable_type]
                except KeyError:
                    pass
        return name

    def register(self, name, alias):
        if not check_alias(alias):
            raise ValueError("Invalid alias '{}' for name '{}'".format(alias, name))
        self.ad[name] = alias

    def register_dict(self, alias_dict):
        for name, alias in alias_dict.items():
            self.register(name, alias)
