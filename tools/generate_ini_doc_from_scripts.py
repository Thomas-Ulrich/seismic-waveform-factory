import ast
import glob
import os
from collections import defaultdict

def get_type_from_func(func_name):
    return {
        "get": "string",
        "getfloat": "float",
        "getboolean": "boolean"
    }.get(func_name, "string")

def extract_config_keys_from_function(fn_node):
    accesses = []
    for node in ast.walk(fn_node):
        if isinstance(node, ast.Call):
            if hasattr(node.func, 'attr') and node.func.attr.startswith("get"):
                func_type = get_type_from_func(node.func.attr)
                args = node.args
                if len(args) >= 2:
                    section = None
                    if isinstance(args[0], ast.Constant):
                        section = args[0].value
                    elif isinstance(args[0], ast.Name):
                        section = f"${args[0].id}"
                    key = args[1].value if isinstance(args[1], ast.Constant) else None
                    required = True
                    default = None
                    for kw in node.keywords:
                        if kw.arg == "fallback":
                            required = False
                            try:
                                default = ast.literal_eval(kw.value)
                            except Exception:
                                default = None
                    accesses.append({
                        "section": section,
                        "key": key,
                        "type": func_type,
                        "required": required,
                        "default": default,
                    })
    return accesses

def extract_config_accesses(path, helper_functions={}):
    with open(path, "r") as f:
        tree = ast.parse(f.read(), filename=path)
    accesses = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            # Direct access
            if hasattr(node.func, 'attr') and node.func.attr.startswith("get"):
                func_type = get_type_from_func(node.func.attr)
                args = node.args
                if len(args) >= 2 and all(isinstance(arg, ast.Constant) for arg in args[:2]):
                    section = args[0].value
                    key = args[1].value
                    required = True
                    default = None
                    for kw in node.keywords:
                        if kw.arg == "fallback":
                            required = False
                            try:
                                default = ast.literal_eval(kw.value)
                            except Exception:
                                default = None
                    accesses.append({
                        "section": section,
                        "key": key,
                        "type": func_type,
                        "required": required,
                        "default": default,
                        "file": os.path.basename(path)
                    })
            # Helper function call
            elif isinstance(node.func, ast.Name) and node.func.id in helper_functions:
                if len(node.args) >= 2:
                    section_arg = node.args[1]
                    if isinstance(section_arg, ast.List):
                        section_names = [
                            elt.value for elt in section_arg.elts if isinstance(elt, ast.Constant)
                        ]
                        for access in helper_functions[node.func.id]:
                            if access['section'] and access['section'].startswith("$"):
                                for sname in section_names:
                                    accesses.append({
                                        "section": sname,
                                        "key": access['key'],
                                        "type": access['type'],
                                        "required": access['required'],
                                        "default": access['default'],
                                        "file": os.path.basename(path)
                                    })
    return accesses

def find_helper_functions():
    helpers = {}
    for script in glob.glob("../scripts/*.py") + glob.glob("scripts/*.py"):
        with open(script, "r") as f:
            tree = ast.parse(f.read(), filename=script)
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                if len(node.args.args) >= 2 and node.args.args[0].arg == "config":
                    accesses = extract_config_keys_from_function(node)
                    if accesses:
                        helpers[node.name] = accesses
    return helpers

def section_heading(text):
    return text + "\n" + "-" * len(text) + "\n\n"

def main():
    helper_functions = find_helper_functions()
    all_accesses = []
    for script in glob.glob("../scripts/*.py") + glob.glob("scripts/*.py"):
        all_accesses.extend(extract_config_accesses(script, helper_functions))

    grouped = defaultdict(lambda: defaultdict(list))
    for acc in all_accesses:
        grouped[acc["section"]][acc["key"]].append(acc)

    doc = "INI Parameter Usage Documentation\n"
    doc += "=================================\n\n"
    for section in sorted(grouped.keys()):
        doc += section_heading(f"[{section}]")
        for key in sorted(grouped[section].keys()):
            items = grouped[section][key]
            types = set(item['type'] for item in items)
            required = all(item['required'] for item in items)
            defaults = set(str(item['default']) for item in items if item['default'] is not None)
            files = set(item['file'] for item in items)
            doc += f"- **{key}**\n"
            doc += f"    - Type: {'/'.join(types)}\n"
            doc += f"    - Required: {'Yes' if required else 'No'}\n"
            doc += f"    - Default: {', '.join(defaults) if defaults else 'None'}\n"
            doc += f"    - Used in: {', '.join(sorted(files))}\n\n"

        # Example INI block for the section
        doc += f"Example [{section}] block:\n\n"
        doc += ".. code-block:: ini\n\n"
        doc += f"    [{section}]\n"
        for key in sorted(grouped[section].keys()):
            items = grouped[section][key]
            default = next((item['default'] for item in items if item['default'] is not None), None)
            # If default is None or string "None"
            if default is None or default == "None":
                doc += f"    {key} = <value>\n"
            else:
                doc += f"    {key} = {default}\n"
        doc += "\n"

    with open("ini_usage_documentation.rst", "w") as f:
        f.write(doc)

if __name__ == "__main__":
    main()
