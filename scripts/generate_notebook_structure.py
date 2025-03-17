import os
import json
import glob

def generate_notebook_structure(docs_dir='docs'):
    """Generate a JSON structure of available notebooks for the UI."""
    notebooks = []
    
    # Get all notebook directories
    notebook_dirs = [d for d in os.listdir(docs_dir) 
                    if os.path.isdir(os.path.join(docs_dir, d))]
    
    for nb_dir in notebook_dirs:
        # Check if default variant exists, use it to get notebook names
        default_dir = os.path.join(docs_dir, nb_dir, 'default')
        if os.path.exists(default_dir):
            html_files = glob.glob(os.path.join(default_dir, '*.html'))
            
            for html_file in html_files:
                notebook_filename = os.path.basename(html_file)
                notebook_name = os.path.splitext(notebook_filename)[0].replace('-', ' ')
                
                # Check if all variants have this notebook
                all_variants_have_file = True
                for variant in ['default', 'naive_reversions_first', 'naive_reversions_no_bp']:
                    variant_path = os.path.join(docs_dir, nb_dir, variant, notebook_filename)
                    if not os.path.exists(variant_path):
                        all_variants_have_file = False
                        break
                
                if all_variants_have_file:
                    notebooks.append({
                        # 'name': f"{nb_dir.replace('-', ' ').title()}: {notebook_name}",
                        'name': nb_dir.replace('-', ' '),
                        'path': f"{nb_dir}",
                        'filename': notebook_filename
                    })
    
    # Write the structure to a JSON file
    with open(os.path.join(docs_dir, 'notebook_structure.json'), 'w') as f:
        json.dump({'notebooks': notebooks}, f, indent=2)
    
    print(f"Generated notebook structure with {len(notebooks)} notebooks")
    
if __name__ == "__main__":
    generate_notebook_structure()
