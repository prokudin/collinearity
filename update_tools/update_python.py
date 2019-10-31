# -*- coding: utf-8 -*-
import os
from Lib.lib2to3.main import refactor
import autopep8
import nbformat
from parse_tools import fix_relative, rename_module


DEFAULT_FIXERS_2TO3 = tuple(sorted(
    set(refactor.get_fixers_from_package("Lib.lib2to3.fixes")) -
    {"Lib.lib2to3.fixes.fix_exitfunc"}
))

DEFAULT_EXPLICIT_2TO3 = ("idioms", "set_literal", "ws_comma")


# =======


def apply_nb(notebook, func, f_args=(), f_kwargs=None):
    """
    Apply func to code in notebook.

    Args:
        notebook: A path to a notebook (.ipynb) file.
        func: A function (str->str) that takes at least one argument.
        f_args, f_kwargs: Args and kwargs passed to func.

    func will be called as 'func(code, *f_args, **f_kwargs)', where 'code' is
    a string containing the code in the notebook.

    The notebook's code will be rewritten by the result of func.
    """
    if f_kwargs is None:
        f_kwargs = {}

    book = nbformat.read(notebook, 4)

    # Join code cells.
    code = "\n# XXX: SPLITHERE\n".join(
        "\n".join(line if not line.startswith("%") else "# " + line
                  for line in cell["source"].splitlines())
        for cell in book["cells"] if cell["cell_type"] == "code"
    )

    code += "\n"  # fixes problems

    # Apply func to code.
    code = func(code, *f_args, **f_kwargs)

    # Split code cells.
    cells = ["\n".join(line if not line.startswith("# %") else line[2:]
                       for line in cell.splitlines())
             for cell in code.split("\n# XXX: SPLITHERE\n")]

    # Replace old code cells with the new ones.
    for i in range(len(book["cells"])):
        if book["cells"][i]["cell_type"] == "code":
            book["cells"][i]["source"] = cells.pop(0)

    # Write changes to notebook.
    nbformat.write(book, notebook)


def apply_py(py, func, f_args=(), f_kwargs=None):
    """
    Apply func to code in python file.

    Args:
        py: A path to a python (.py) file.
        func: A function (str->str) that takes at least one argument.
        f_args, f_kwargs: Args and kwargs passed to func.

    func will be called as 'func(code, *f_args, **f_kwargs)', where 'code' is
    a string containing the code in the notebook.

    The file will be rewritten by the result of func.
    """
    if f_kwargs is None:
        f_kwargs = {}

    with open(py, "r") as f:
        code = f.read()

    code = func(code, *f_args, **f_kwargs)

    with open(py, "w") as f:
        f.write(code)


# =======


def get_files(root, ext):
    """Walks from root and yields paths of files with extension ext."""
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            if filename.endswith(ext):
                yield os.path.join(dirpath, filename)


# =======


def apply_files(files, func, f_args=(), f_kwargs=None, ignore_err=False):
    """
    Apply func to code in each file in files.

    Args:
        files: An iterable of paths to .py and/or .ipynb files.
        func: A function (str->str) that takes at least one argument.
        f_args, f_kwargs: Args and kwargs passed to func.
        ignore_err: If true, apply_files will move on to the next file if the
                    the current file results in an exception.

    The files will be rewritten with the results of func.
    """
    for file in files:
        try:
            if file.endswith(".py"):
                apply_py(file, func, f_args, f_kwargs)
            elif file.endswith(".ipynb"):
                apply_nb(file, func, f_args, f_kwargs)
            else:
                raise ValueError("file path should end in .py or .ipynb")
        except Exception as err:
            print("ERROR: func:", func.__name__)
            print("ERROR: file:", file)
            print("ERROR:", err)
            if ignore_err:
                continue
            raise err


# =======


def convert2to3(code, py2to3tool=None, fixers=None, explicit=None):
    """
    Use lib2to3 to convert code from python 2 to python 3.

    Returns:
        A string containing the code converted to python 3.
    """
    try:
        compile(code, "", "exec")  # test if already py3-friendly
        return code
    except Exception:
        pass

    if py2to3tool is None:
        if fixers is None:
            fixers = DEFAULT_FIXERS_2TO3
        if explicit is None:
            explicit = DEFAULT_EXPLICIT_2TO3

        py2to3tool = refactor.RefactoringTool(fixers, explicit=explicit)

    return str(py2to3tool.refactor_string(code, ""))


def fix_style(code, options=None):
    """
    Use autopep8 to stylize code.

    Returns:
        A string containing stylized code.
    """
    if options is None:
        options = {"select": ["E2"]}

    return str(autopep8.fix_code(code, options=options))


if __name__ == "__main__":

    # TODO: * switch griddata from plt to scipy
    #       * change line endings

    target_dir = ".."

    pys = list(get_files(target_dir, ".py"))
    nbs = list(get_files(target_dir, ".ipynb"))

    files = pys + nbs
    # files = pys
    # files = nbs

    # apply_files(files, convert2to3, ignore_err=True)
    apply_files(files, fix_relative, (target_dir,), ignore_err=True)
    apply_files(files, rename_module, ("aux", "auxiliary"), ignore_err=True)
    apply_files(files, fix_style, ignore_err=True)
