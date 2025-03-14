import nbconvert
import os
import glob

def convert_notebooks(results_dir='results/notebooks'):
    # Find all notebook directories
    notebook_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    
    for nb_dir in notebook_dirs:
        variants = ['default', 'naive_reversions_first', 'naive_reversions_no_bp']
        for variant in variants:
            variant_path = os.path.join(results_dir, nb_dir, variant)
            if os.path.exists(variant_path):
                notebooks = glob.glob(os.path.join(variant_path, '*.ipynb'))
                for notebook in notebooks:
                    # Convert to HTML
                    html_exporter = nbconvert.HTMLExporter()
                    output_html, resources = html_exporter.from_filename(notebook)
                    
                    # Create output directory
                    output_dir = os.path.join('docs', nb_dir, variant)
                    os.makedirs(output_dir, exist_ok=True)
                    
                    # Write HTML file
                    base_name = os.path.basename(notebook).replace('.ipynb', '.html')
                    with open(os.path.join(output_dir, base_name), 'w') as f:
                        f.write(output_html)
                    
                    # Copy associated resources
                    # resource_files = glob.glob(os.path.join(os.path.dirname(notebook), '*.pdf')) + \
                    #               glob.glob(os.path.join(os.path.dirname(notebook), '*.csv'))
                    # for resource in resource_files:
                    #     resource_name = os.path.basename(resource)
                    #     os.system(f"cp {resource} {output_dir}/{resource_name}")

if __name__ == "__main__":
    convert_notebooks()
