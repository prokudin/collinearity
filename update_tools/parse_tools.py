# -*- coding: utf-8 -*-
import os
import regex

# Capture groups: MOD_TERM, REL_MOD, AS_PAIR, MODULE, ALIAS, ID

subs = [
    {"<mod_term>": "(?P<MOD_TERM>(<as_pair>|<module>))"},
    {"<relative_module>":
        r"(?P<REL_MOD>(\.*<module>|\.+))",
     "<as_pair>":
         r"(?P<AS_PAIR><module> +as +<alias>)"},
    {"<module>": r"(?P<MODULE>(<identifier>\.)*<identifier>)"},
    {"<alias>": r"(?P<ALIAS>(<identifier>))"},
    {"<identifier>": r"(?P<ID><id_start><id_continue>*)"},
    {"<id_start>":
        r"[\p{Lu}\p{Ll}\p{Lt}\p{Lm}\p{Lo}\p{Nl}_\p{Other_ID_Start}]",
     "<id_continue>":
        r"[\p{Lu}\p{Ll}\p{Lt}\p{Lm}\p{Lo}\p{Nl}_\p{Other_ID_Start}"
        r"\p{Mn}\p{Mc}\p{Nd}\p{Pc}\p{Other_ID_Continue}]"}
]

import_stmt = (
    r"(import +<mod_term>( *, *<mod_term>)*)"
    r"|(from +<relative_module> +import +<mod_term>( *, *<mod_term>)*)"
    r"|(from +<relative_module> +import +\( +<mod_term>"
    r"( *, *<mod_term>)* *,? +\))"
    r"|(from +<module> +import +\*)"
)

identifier = "<identifier>"

for _replacements in subs:
    for _name, _sub in _replacements.items():
        import_stmt = import_stmt.replace(_name, _sub)
        identifier = identifier.replace(_name, _sub)

import_re = regex.compile(import_stmt)
id_re = regex.compile(identifier)


def fix_relative(code, dirpath):
    """
    Fix cases where code transformations (such as lib2to3) mess up imports.

    Specifically, if a relative import could instead be an absolute import
    when dirpath is treated as part of the module's path, then change the
    import from relative to absolute.

    Examples (where 'tools' exists in dir at dirpath):
        1) "from .tools.config import conf" ->
           "from tools.config import conf"

        2) "from . import tools.config" ->
           "import tools.config"

    Warning:
        This will remove intentional relative imports from files at dirpath.

    Returns:
        A string containing the code with fixed imports.
    """
    dirs = next(os.walk(dirpath))[1]

    lines = []

    for line in code.splitlines():
        match = import_re.match(line)

        if match:
            rel_mod = match.captures("REL_MOD")

            if rel_mod:
                rel_mod = rel_mod[0]
                s = match.starts("REL_MOD")[0]

                if rel_mod.startswith("."):
                    if rel_mod == ".":
                        # "from . import x.y" -> "import x.y"
                        if match.captures("MODULE")[0].split(".")[0] in dirs:
                            line = line[line.find("import"):]
                    else:
                        parts = rel_mod.split(".")
                        if parts[1] in dirs:
                            line = line[:s] + line[s + 1:]  # get rid of '.'

        lines.append(line)

    return "\n".join(lines)


def rename_module(code, old, new):
    """
    Find and rename a module within the given code.

    If the old module is imported directly and not given an alias, then it will
    be renamed throughout the rest of the code.

    Warning:
        If the name of the old module is to be replaced throughout the code,
        then instances of the name will be replaced even if the name is found
        within strings or has been reassigned.

    Returns:
        A string containing the updated code.
    """
    replace_throughout = False

    lines = []

    for line in code.splitlines():
        match = import_re.match(line)

        if match:
            mod_terms = match.captures("MOD_TERM")
            modules = match.captures("MODULE")

            for module in modules:
                if "." in module:
                    if old in module.split("."):
                        line = line.replace(
                            module,
                            ".".join((new if m == old else m)
                                     for m in module.split(".")),
                            1
                        )
                        break
                elif module == old:
                    # If imported directly and not aliased, replace throughout.
                    replace_throughout = (old in mod_terms)

                    for g, s in zip(match.captures("MODULE"),
                                    match.starts("MODULE")):
                        if g == old:
                            line = (line[:s] +
                                    line[s:].replace(old, new, 1))
                            break
                    break

        elif replace_throughout:
            while old in id_re.findall(line):
                for match in id_re.finditer(line):
                    s = slice(*match.span())
                    if line[s] == old:
                        line = line[:s.start] + new + line[s.stop:]
                        break

        lines.append(line)

    return "\n".join(lines)


if __name__ == "__main__":
    tests = ("import os",
             "import os, sys",
             "import os,sys",
             "import itertools as it",
             "import itertools as it, functools as ft",
             "import matplotlib.pyplot as plt",
             "from os import walk",
             "from os import walk as run",
             "from os",
             "from os import *",
             "from qcdlib.aux import AUX",
             "from aux import AUX",
             "from auxiliary import AUX",
             "from aux import AUX  # comment",
             "import aux as test",
             "from aux import *",
             "import aux as allah, auxiliary",
             "import auxiliary as allah, aux",  # XXX
             "import auxiliary as allah, aux as ahmed",
             "import aux\nprint(aux)")

    for test in tests:
        print("test:", test)
        print("match:", import_re.match(test))
        print("new:", rename_module(test, "aux", "auxiliary"))
        print("-------")
