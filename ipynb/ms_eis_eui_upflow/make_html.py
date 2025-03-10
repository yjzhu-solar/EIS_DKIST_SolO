from ipy2html import ipy2html
from glob import glob
import os
import pandas as pd
import sys
from markdown_it import MarkdownIt
from pathlib import Path


files = glob('*.ipynb')
additional_files = ['../aia/aia_coalign_1024_final_rsuncorr.ipynb', '../spice_psf/ipynb/spice_1024_NeVIII_ms_version.ipynb']
files.extend(additional_files)
output_dir = '../../eis_eui_upflow_ipynb_html/'

if os.path.exists('latest_mod_date.csv'):
    latest_mod_date = pd.read_csv('latest_mod_date.csv')

    file_name_new = []
    mod_date_new = []
    file_name_updated = []
    mod_date_updated = []
    
    for file in files:
        if file not in latest_mod_date['file_name'].values:
            file_name_new.append(file)
            mod_date_new.append(os.path.getmtime(file))
            ipy2html(file, output_dir)
        else:
            index = latest_mod_date[latest_mod_date['file_name'] == file].index[0]
            if os.path.getmtime(file) - latest_mod_date.loc[index, 'mod_date'] > 0.2:
                ipy2html(file, output_dir)
                latest_mod_date.loc[index, 'mod_date'] = os.path.getmtime(file)
                file_name_updated.append(file)
                mod_date_updated.append(os.path.getmtime(file))
            else:
                pass
    
    if len(file_name_new) > 0:
        latest_mod_date_new = pd.DataFrame({'file_name': file_name_new, 'mod_date': mod_date_new})
        latest_mod_date = pd.concat([latest_mod_date, latest_mod_date_new], axis=0)
    latest_mod_date.to_csv('latest_mod_date.csv', index=False)

    print('Updated files:', [os.path.basename(file) for file in file_name_updated])
    print('New files:', [os.path.basename(file) for file in file_name_new])

    
else:
    file_name = []
    mod_date = []
    for file in files:
        file_name.append(file)
        mod_date.append(os.path.getmtime(file))
        ipy2html(file, output_dir)
    
    latest_mod_date = pd.DataFrame({'file_name': file_name, 'mod_date': mod_date})
    latest_mod_date.to_csv('latest_mod_date.csv', index=False)

    print('New files:', [os.path.basename(file) for file in file_name])

# generate index.html from README.md
# code copied from Claude

def convert_markdown_to_html(markdown_content):
    """Convert markdown content to HTML using markdown-it-py."""
    md = MarkdownIt()
    return md.render(markdown_content)

def read_markdown_file(file_path):
    """Read content from a markdown file."""
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

def write_html_file(file_path, html_content):
    """Write HTML content to a file."""
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

def make_index(input_file, output_file):
    try:
        # Read markdown content
        markdown_content = read_markdown_file(input_file)
        
        # Convert to HTML
        html_content = convert_markdown_to_html(markdown_content)
        
        # Generate a simple HTML document
        full_html = f"""<!DOCTYPE html>
                    <html>
                    <head>
                        <meta charset="utf-8">
                        <meta name="viewport" content="width=device-width, initial-scale=1">
                        <title>{Path(input_file).stem}</title>
                        <style>
                            body {{
                                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
                                line-height: 1.6;
                                max-width: 800px;
                                margin: 0 auto;
                                padding: 20px;
                            }}
                            pre {{
                                background-color: #f5f5f5;
                                padding: 12px;
                                border-radius: 4px;
                                overflow-x: auto;
                            }}
                            code {{
                                font-family: Monaco, Consolas, "Courier New", monospace;
                            }}
                        </style>
                    </head>
                    <body>
                        {html_content}
                    </body>
                    </html>"""
        
        # Write to output file
        write_html_file(output_file, full_html)
        
        print(f"Converted {input_file} to {output_file}")
    
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if os.path.exists('../../index.html'):
    # delete the existing index.html
    os.remove('../../index.html')

make_index('../../README.md', '../../index.html')


    

