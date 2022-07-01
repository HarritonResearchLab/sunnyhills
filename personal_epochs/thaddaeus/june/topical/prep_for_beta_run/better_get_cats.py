def print_them_links():
    import requests
    from bs4 import BeautifulSoup
    
    
    url = 'https://archive.stsci.edu/tess/tic_ctl.html'
    reqs = requests.get(url)
    soup = BeautifulSoup(reqs.text, 'html.parser')
    
    urls = []
    for link in soup.find_all('a'):
        print(link.get('href')) 

def gunzip(dir):
    import os
    
    os.chdir(dir)

    for zipped_file in os.listdir(dir): 
        os.system('gunzip ' + zipped_file)

gunzip('/ar1/PROJ/fjuhsd/shared/bigdata/tic_catalog_v8')