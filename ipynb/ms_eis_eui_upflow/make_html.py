from ipy2html import ipy2html
from glob import glob
import os
import pandas as pd


files = glob('*.ipynb')
output_dir = '../../eis_eui_upflow_ipynb_html/'

if os.path.exists('latest_mod_date.csv'):
    latest_mod_date = pd.read_csv('latest_mod_date.csv')

    for file in files:
        file_name_new = []
        mod_date_new = []
        file_name_updated = []
        mod_date_updated = []
        if file not in latest_mod_date['file_name'].values:
            file_name_new.append(file)
            mod_date_new.append(os.path.getmtime(file))
            
            ipy2html(file, output_dir)
        else:
            index = latest_mod_date[latest_mod_date['file_name'] == file].index[0]
            if os.path.getmtime(file) > latest_mod_date.loc[index, 'mod_date']:
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
    

