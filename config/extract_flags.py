# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-09-10
# @description: Extract and document all flags used in a _utility_imports.py script
# =============================================================================

import ast
import re

class FlagVisitor(ast.NodeVisitor):
    def __init__(self, source_lines):
        self.cli_flags = []
        self.internal_flags = []
        self.yaml_flags = []
        self.dict_params = []
        self.other_vars = []
        self.source_lines = source_lines  # Store source lines for comment parsing

    # ----------------------------
    # Helpers (Py 3.8+ compatible)
    # ----------------------------
    def _const_value(self, node):
        """Return Python value if node is a literal constant, else None."""
        if isinstance(node, ast.Constant):
            return node.value
        return None

    def _infer_type_and_default(self, node):
        """
        Infer (default_str, type_str) from an AST node.
        Covers: Constant, List, Dict, Call, Name, Attribute, etc.
        """
        if isinstance(node, ast.Constant):
            v = node.value
            return str(v), type(v).__name__
        if isinstance(node, ast.List):
            return "[]", "list"
        if isinstance(node, ast.Tuple):
            return "()", "tuple"
        if isinstance(node, ast.Dict):
            return "{}", "dict"
        if isinstance(node, ast.Call):
            return "Computed", "Unknown"
        if isinstance(node, ast.Name):
            return f"Computed({node.id})", "Unknown"
        if isinstance(node, ast.Attribute):
            # e.g., module.CONST
            parts = []
            cur = node
            while isinstance(cur, ast.Attribute):
                parts.append(cur.attr)
                cur = cur.value
            if isinstance(cur, ast.Name):
                parts.append(cur.id)
            parts.reverse()
            return f"Computed({'.'.join(parts)})", "Unknown"
        return "Unknown", "Unknown"

    def _subscript_key(self, subscript_node):
        """
        Extract key from params['key'] in a version-safe way.
        Returns the literal key if available; otherwise 'Computed'.
        """
        slc = subscript_node.slice  # py3.9+ uses expression directly here
        if isinstance(slc, ast.Constant):
            return slc.value
        return "Computed"

    def _get_comments(self, lineno):
        """Extract inline or preceding comments for a given line number."""
        comments = []
        if lineno is None or lineno <= 0 or lineno > len(self.source_lines):
            return None

        # Inline comment (same line)
        line = self.source_lines[lineno - 1].strip()
        inline_match = re.search(r"#\s*(.+)$", line)
        if inline_match:
            comments.append(inline_match.group(1).strip())

        # Preceding comment block
        for i in range(lineno - 2, -1, -1):
            prev_line = self.source_lines[i].strip()
            if prev_line.startswith("#"):
                comments.append(prev_line.lstrip("#").strip())
            elif prev_line:  # Stop at first non-comment, non-empty line
                break

        return " ".join(reversed(comments)) if comments else None

    def _generate_description(self, name, source, comment=None):
        """Generate a description based on name, source, and optional comment."""
        if comment:
            return comment

        name_str = str(name)
        name_lower = name_str.lower()

        if "flag" in name_lower:
            return f"Controls {name_lower.replace('_', ' ')} behavior in script logic"
        if source == "CLI":
            return f"Command-line argument for {name_lower.replace('_', ' ')}"
        if source == "YAML":
            return f"YAML configuration parameter for {name_lower.replace('_', ' ')}"
        if source == "Dictionary":
            return f"Parameter for {name_lower.replace('_', ' ')} in dictionary"
        if source == "Variable":
            return f"Variable used for {name_lower.replace('_', ' ')} in script logic"
        if source == "Internal":
            return f"Internal parameter controlling {name_lower.replace('_', ' ')}"
        return f"Parameter {name_lower.replace('_', ' ')}"

    # ----------------------------
    # Visitors
    # ----------------------------
    def visit_Assign(self, node):
        comment = self._get_comments(getattr(node, "lineno", None))

        # Internal flags and other tracked variables
        for target in node.targets:
            # Internal flags (variables with "flag" in name)
            if isinstance(target, ast.Name) and "flag" in target.id.lower():
                default, flag_type = self._infer_type_and_default(node.value)
                self.internal_flags.append({
                    "name": target.id,
                    "description": self._generate_description(target.id, "Internal", comment),
                    "type": flag_type,
                    "default": default,
                    "required": "No",
                    "choices": "None",
                    "source": "Internal",
                    "line": getattr(node, "lineno", None),
                })

            # Other variables of interest (e.g., params_vec, observed_params)
            elif isinstance(target, ast.Name) and target.id.lower() in [
                "params_vec", "observed_params", "joint_estimated_params"
            ]:
                default, flag_type = self._infer_type_and_default(node.value)
                self.other_vars.append({
                    "name": target.id,
                    "description": self._generate_description(target.id, "Variable", comment),
                    "type": flag_type,
                    "default": default,
                    "required": "No",
                    "choices": "None",
                    "source": "Variable",
                    "line": getattr(node, "lineno", None),
                })

        # Dictionary assignments: params['key'] = value
        for tgt in node.targets:
            if isinstance(tgt, ast.Subscript) and isinstance(tgt.value, ast.Name):
                dict_name = tgt.value.id
                if dict_name in ["params", "kwargs", "execution_mode"]:
                    key = self._subscript_key(tgt)
                    default, flag_type = self._infer_type_and_default(node.value)

                    self.dict_params.append({
                        "name": key,
                        "description": self._generate_description(key, "Dictionary", comment),
                        "type": flag_type,
                        "default": default,
                        "required": "No",
                        "choices": "None",
                        "source": "Dictionary",
                        "line": getattr(node, "lineno", None),
                    })

        self.generic_visit(node)

    def visit_Call(self, node):
        comment = self._get_comments(getattr(node, "lineno", None))

        # CLI flags from parser.add_argument
        if isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument":
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "parser":
                flag_name = "Unknown"
                description = "No description provided"
                default = "None"
                arg_type = "str"
                required = "No"
                choices = "None"

                if node.args:
                    if isinstance(node.args[0], ast.Constant):
                        flag_name = node.args[0].value
                    elif isinstance(node.args[0], ast.Tuple):
                        parts = []
                        for elt in node.args[0].elts:
                            if isinstance(elt, ast.Constant):
                                parts.append(str(elt.value))
                        if parts:
                            flag_name = ", ".join(parts)

                for kw in node.keywords:
                    if kw.arg == "help":
                        v = self._const_value(kw.value)
                        if v is not None:
                            description = v
                    elif kw.arg == "default":
                        v = self._const_value(kw.value)
                        if v is not None:
                            default = str(v)
                        else:
                            default = self._infer_type_and_default(kw.value)[0]
                    elif kw.arg == "type":
                        if isinstance(kw.value, ast.Name):
                            arg_type = kw.value.id
                        else:
                            arg_type = "Unknown"
                    elif kw.arg == "required":
                        v = self._const_value(kw.value)
                        if isinstance(v, bool):
                            required = "Yes" if v else "No"
                    elif kw.arg == "choices" and isinstance(kw.value, ast.List):
                        parts = []
                        for elt in kw.value.elts:
                            if isinstance(elt, ast.Constant):
                                parts.append(str(elt.value))
                        choices = ", ".join(parts) if parts else "None"

                self.cli_flags.append({
                    "name": flag_name,
                    "description": description,
                    "type": arg_type,
                    "default": default,
                    "required": required,
                    "choices": choices,
                    "source": "CLI",
                    "line": getattr(node, "lineno", None),
                })

        # YAML parameters from .get() calls
        elif isinstance(node.func, ast.Attribute) and node.func.attr == "get":
            if isinstance(node.func.value, ast.Name) and node.func.value.id in [
                "enkf_params", "modeling_params", "physical_params", "params", "kwargs"
            ]:
                if len(node.args) >= 1 and isinstance(node.args[0], ast.Constant):
                    key = node.args[0].value
                    default = "None"
                    flag_type = "Unknown"

                    if len(node.args) > 1:
                        default, flag_type = self._infer_type_and_default(node.args[1])

                    self.yaml_flags.append({
                        "name": key,
                        "description": self._generate_description(key, "YAML", comment),
                        "type": flag_type,
                        "default": default,
                        "required": "No",
                        "choices": "None",
                        "source": "YAML",
                        "line": getattr(node, "lineno", None),
                    })

        # Dictionary updates: kwargs.update({...})
        if isinstance(node.func, ast.Attribute) and node.func.attr == "update":
            if isinstance(node.func.value, ast.Name):
                dict_name = node.func.value.id
                if dict_name in ["params", "kwargs"]:
                    if len(node.args) == 1 and isinstance(node.args[0], ast.Dict):
                        d = node.args[0]
                        for k_node, v_node in zip(d.keys, d.values):
                            key = k_node.value if isinstance(k_node, ast.Constant) else "Computed"
                            default, flag_type = self._infer_type_and_default(v_node)
                            self.dict_params.append({
                                "name": key,
                                "description": self._generate_description(key, "Dictionary", comment),
                                "type": flag_type,
                                "default": default,
                                "required": "No",
                                "choices": "None",
                                "source": "Dictionary",
                                "line": getattr(node, "lineno", None),
                            })

                    elif len(node.args) == 1 and isinstance(node.args[0], ast.Name):
                        dict_ref = node.args[0].id
                        if dict_ref in ["physical_params", "modeling_params", "enkf_params"]:
                            self.dict_params.append({
                                "name": f"{dict_ref}_keys",
                                "description": self._generate_description(f"{dict_ref}_keys", "Dictionary", comment),
                                "type": "dict",
                                "default": "Unknown",
                                "required": "No",
                                "choices": "None",
                                "source": "Dictionary",
                                "line": getattr(node, "lineno", None),
                            })

        self.generic_visit(node)


def extract_flags(script_path):
    with open(script_path, "r", encoding="utf-8") as f:
        source_lines = f.readlines()

    tree = ast.parse("".join(source_lines))
    visitor = FlagVisitor(source_lines)
    visitor.visit(tree)

    # Combine and deduplicate by name
    all_flags = (
        visitor.cli_flags
        + visitor.internal_flags
        + visitor.yaml_flags
        + visitor.dict_params
        + visitor.other_vars
    )
    unique_flags = {flag["name"]: flag for flag in all_flags}.values()
    return sorted(unique_flags, key=lambda x: str(x["name"]))


def generate_flags_markdown(flags):
    doc_lines = [
        "## All Main Flags used in ICESEE \n",
        "| Name | Description | Type | Default | Required | Choices | Source |\n",
        "|------|-------------|------|---------|----------|---------|--------|\n",
    ]
    for flag in flags:
        description = str(flag["description"]).replace("|", "&#124;").replace("\n", " ")
        doc_lines.append(
            f"| `{flag['name']}` | {description} | {flag['type']} | {flag['default']} | "
            f"{flag['required']} | {flag['choices']} | {flag['source']} |\n"
        )
    return "".join(doc_lines)


if __name__ == "__main__":
    flags = extract_flags("config/_utility_imports.py")
    markdown = generate_flags_markdown(flags)
    print(markdown)