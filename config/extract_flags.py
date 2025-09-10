# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-09-10
# @description: Class to extract and document all flags used in a _utility_imports.py scripts
# =============================================================================

import ast
import os
import re

class FlagVisitor(ast.NodeVisitor):
    def __init__(self, source_lines):
        self.cli_flags = []
        self.internal_flags = []
        self.yaml_flags = []
        self.dict_params = []
        self.other_vars = []
        self.source_lines = source_lines  # Store source lines for comment parsing

    def _get_comments(self, lineno):
        """Extract inline or preceding comments for a given line number."""
        comments = []
        # Check inline comment (same line)
        line = self.source_lines[lineno - 1].strip()
        inline_match = re.search(r'#\s*(.+)$', line)
        if inline_match:
            comments.append(inline_match.group(1).strip())
        
        # Check preceding lines for comments
        for i in range(lineno - 2, -1, -1):
            prev_line = self.source_lines[i].strip()
            if prev_line.startswith('#'):
                comments.append(prev_line.lstrip('#').strip())
            elif prev_line.strip():  # Stop at first non-comment, non-empty line
                break
        return ' '.join(reversed(comments)) if comments else None

    def _generate_description(self, name, source, comment=None):
        """Generate a description based on name, source, and optional comment."""
        if comment:
            return comment
        # Contextual descriptions based on name
        name_lower = name.lower()
        if 'flag' in name_lower:
            return f"Controls {name_lower.replace('_', ' ')} behavior in script logic"
        elif source == 'CLI':
            return f"Command-line argument for {name_lower.replace('_', ' ')}"
        elif source == 'YAML':
            return f"YAML configuration parameter for {name_lower.replace('_', ' ')}"
        elif source == 'Dictionary':
            return f"Parameter for {name_lower.replace('_', ' ')} in dictionary"
        elif source == 'Variable':
            return f"Variable used for {name_lower.replace('_', ' ')} in script logic"
        return f"Parameter {name_lower.replace('_', ' ')}"

    def visit_Assign(self, node):
        comment = self._get_comments(node.lineno)
        # Internal flags (variables with "flag" in name)
        for target in node.targets:
            if isinstance(target, ast.Name) and 'flag' in target.id.lower():
                default = 'Unknown'
                flag_type = 'Unknown'
                if isinstance(node.value, ast.Constant):
                    default = str(node.value.value)
                    flag_type = type(node.value.value).__name__
                elif isinstance(node.value, ast.NameConstant):
                    default = str(node.value.value)
                    flag_type = 'bool' if isinstance(node.value.value, bool) else 'Unknown'
                self.internal_flags.append({
                    'name': target.id,
                    'description': self._generate_description(target.id, 'Internal', comment),
                    'type': flag_type,
                    'default': default,
                    'required': 'No',
                    'choices': 'None',
                    'source': 'Internal',
                    'line': node.lineno
                })
            # Other variables (e.g., params_vec, observed_params)
            elif isinstance(target, ast.Name) and target.id.lower() in ['params_vec', 'observed_params', 'joint_estimated_params']:
                default = 'Unknown'
                flag_type = 'Unknown'
                if isinstance(node.value, ast.List):
                    default = '[]'
                    flag_type = 'list'
                self.other_vars.append({
                    'name': target.id,
                    'description': self._generate_description(target.id, 'Variable', comment),
                    'type': flag_type,
                    'default': default,
                    'required': 'No',
                    'choices': 'None',
                    'source': 'Variable',
                    'line': node.lineno
                })
        # Dictionary assignments (e.g., params['key'] = value)
        if isinstance(node.targets[0], ast.Subscript) and isinstance(node.targets[0].value, ast.Name):
            dict_name = node.targets[0].value.id
            if dict_name in ['params', 'kwargs', 'execution_mode']:
                key = node.targets[0].slice.value if isinstance(node.targets[0].slice, ast.Constant) else 'Computed'
                default = 'Unknown'
                flag_type = 'Unknown'
                if isinstance(node.value, ast.Constant):
                    default = str(node.value.value)
                    flag_type = type(node.value.value).__name__
                elif isinstance(node.value, ast.NameConstant):
                    default = str(node.value.value)
                    flag_type = 'bool' if isinstance(node.value.value, bool) else 'Unknown'
                elif isinstance(node.value, ast.List):
                    default = '[]'
                    flag_type = 'list'
                elif isinstance(node.value, ast.Call):
                    default = 'Computed'
                    flag_type = 'Unknown'
                self.dict_params.append({
                    'name': key,
                    'description': self._generate_description(key, 'Dictionary', comment),
                    'type': flag_type,
                    'default': default,
                    'required': 'No',
                    'choices': 'None',
                    'source': 'Dictionary',
                    'line': node.lineno
                })
        # Dictionary updates (e.g., kwargs.update({'key': value}))
        elif isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == 'update':
            dict_name = node.value.func.value.id
            if dict_name in ['params', 'kwargs']:
                if len(node.value.args) == 1 and isinstance(node.value.args[0], ast.Dict):
                    for key, value in zip(node.value.args[0].keys, node.value.args[0].values):
                        if isinstance(key, ast.Constant):
                            default = 'Unknown'
                            flag_type = 'Unknown'
                            if isinstance(value, ast.Constant):
                                default = str(value.value)
                                flag_type = type(value.value).__name__
                            elif isinstance(value, ast.NameConstant):
                                default = str(value.value)
                                flag_type = 'bool' if isinstance(value.value, bool) else 'Unknown'
                            elif isinstance(value, ast.List):
                                default = '[]'
                                flag_type = 'list'
                            self.dict_params.append({
                                'name': key.value,
                                'description': self._generate_description(key.value, 'Dictionary', comment),
                                'type': flag_type,
                                'default': default,
                                'required': 'No',
                                'choices': 'None',
                                'source': 'Dictionary',
                                'line': node.lineno
                            })
                elif len(node.value.args) == 1 and isinstance(node.value.args[0], ast.Name):
                    dict_ref = node.value.args[0].id
                    if dict_ref in ['physical_params', 'modeling_params', 'enkf_params']:
                        self.dict_params.append({
                            'name': f'{dict_ref}_keys',
                            'description': self._generate_description(f'{dict_ref}_keys', 'Dictionary', comment),
                            'type': 'dict',
                            'default': 'Unknown',
                            'required': 'No',
                            'choices': 'None',
                            'source': 'Dictionary',
                            'line': node.lineno
                        })
        self.generic_visit(node)

    def visit_Call(self, node):
        comment = self._get_comments(node.lineno)
        # CLI flags from add_argument
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'add_argument':
            if isinstance(node.func.value, ast.Name) and node.func.value.id == 'parser':
                flag_name = 'Unknown'
                description = 'No description provided'
                default = 'None'
                arg_type = 'str'
                required = 'No'
                choices = 'None'
                for kw in node.keywords:
                    if kw.arg == 'help' and isinstance(kw.value, ast.Constant):
                        description = kw.value.value
                    elif kw.arg == 'default' and isinstance(kw.value, (ast.Constant, ast.NameConstant)):
                        default = str(kw.value.value)
                    elif kw.arg == 'type' and isinstance(kw.value, ast.Name):
                        arg_type = kw.value.id
                    elif kw.arg == 'required' and isinstance(kw.value, ast.Constant):
                        required = 'Yes' if kw.value.value else 'No'
                    elif kw.arg == 'choices' and isinstance(kw.value, ast.List):
                        choices = ', '.join(str(elt.value) for elt in kw.value.elts if isinstance(elt, ast.Constant))
                if node.args and isinstance(node.args[0], (ast.Constant, ast.Tuple)):
                    flag_name = node.args[0].value if isinstance(node.args[0], ast.Constant) else ', '.join(str(elt.value) for elt in node.args[0].elts if isinstance(elt, ast.Constant))
                self.cli_flags.append({
                    'name': flag_name,
                    'description': description,
                    'type': arg_type,
                    'default': default,
                    'required': required,
                    'choices': choices,
                    'source': 'CLI',
                    'line': node.lineno
                })
        # YAML parameters from .get() calls
        elif isinstance(node.func, ast.Attribute) and node.func.attr == 'get':
            if isinstance(node.func.value, ast.Name) and node.func.value.id in ['enkf_params', 'modeling_params', 'physical_params', 'params', 'kwargs']:
                if len(node.args) >= 1 and isinstance(node.args[0], ast.Constant):
                    key = node.args[0].value
                    default = 'None'
                    flag_type = 'Unknown'
                    if len(node.args) > 1:
                        if isinstance(node.args[1], ast.Constant):
                            default = str(node.args[1].value)
                            flag_type = type(node.args[1].value).__name__
                        elif isinstance(node.args[1], ast.NameConstant):
                            default = str(node.args[1].value)
                            flag_type = 'bool' if isinstance(node.args[1].value, bool) else 'Unknown'
                        elif isinstance(node.args[1], ast.List):
                            default = '[]'
                            flag_type = 'list'
                    self.yaml_flags.append({
                        'name': key,
                        'description': self._generate_description(key, 'YAML', comment),
                        'type': flag_type,
                        'default': default,
                        'required': 'No',
                        'choices': 'None',
                        'source': 'YAML',
                        'line': node.lineno
                    })
        self.generic_visit(node)

def extract_flags(script_path):
    with open(script_path, 'r') as f:
        source_lines = f.readlines()
    tree = ast.parse(''.join(source_lines))
    visitor = FlagVisitor(source_lines)
    visitor.visit(tree)
    
    # Combine and deduplicate by name
    all_flags = visitor.cli_flags + visitor.internal_flags + visitor.yaml_flags + visitor.dict_params + visitor.other_vars
    unique_flags = {flag['name']: flag for flag in all_flags}.values()  # Remove duplicates
    return sorted(unique_flags, key=lambda x: x['name'])

def generate_flags_markdown(flags):
    doc_lines = [
        "## All Main Flags used in ICESEE \n",
        "| Name | Description | Type | Default | Required | Choices | Source |\n",
        "|------|-------------|------|---------|----------|---------|--------|\n"
    ]
    for flag in flags:
        # Clean description to avoid Markdown issues
        description = flag['description'].replace('|', '&#124;').replace('\n', ' ')
        doc_lines.append(
            f"| `{flag['name']}` | {description} | {flag['type']} | {flag['default']} | {flag['required']} | {flag['choices']} | {flag['source']} |\n"
        )
    return "".join(doc_lines)

# For testing
if __name__ == "__main__":
    flags = extract_flags('_utility_imports.py')  # Replace with your script's path
    markdown = generate_flags_markdown(flags)
    print(markdown)