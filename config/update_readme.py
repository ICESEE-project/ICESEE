import os
import re
from extract_flags import extract_flags, generate_flags_markdown

def get_project_root():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    while not os.path.exists(os.path.join(current_dir, 'src')) and current_dir != '/':
        current_dir = os.path.dirname(current_dir)
    return current_dir

def update_readme(script_path):
    flags = extract_flags(script_path)
    flags_doc = generate_flags_markdown(flags)
    readme_path = os.path.join(get_project_root(), 'config', 'README.md')

    try:
        with open(readme_path, 'r') as f:
            readme_content = f.read()
    except FileNotFoundError:
        readme_content = "# Project README\n\n"

    # Regex to match the entire "All Main Flags Used in ICESEE" section
    flags_section_regex = r'(## All Main Flags Used in ICESEE\n.*?)(?=\n## |\Z)'
    if re.search(flags_section_regex, readme_content, re.DOTALL):
        # Replace the existing section
        new_content = re.sub(flags_section_regex, flags_doc, readme_content, flags=re.DOTALL)
    else:
        # Append if no section exists
        new_content = readme_content.rstrip() + '\n\n' + flags_doc

    with open(readme_path, 'w') as f:
        f.write(new_content)

if __name__ == "__main__":
    script_path = os.path.join(get_project_root(), 'config', '_utility_imports.py')  # Adjust to your script's path
    print("Updating README with flags from", script_path)
    update_readme(script_path)