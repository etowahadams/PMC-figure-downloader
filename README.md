# PMC Figure Downloader

This is a simple script to download figures from open access papers in [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/) database. It uses the [PMC API](https://www.ncbi.nlm.nih.gov/pmc/tools/developers/) to search for articles and uses the 
[PMC Open Access Web Service API](https://www.ncbi.nlm.nih.gov/pmc/tools/oa-service/) to get the XML of each paper to determine the figure URLs.

## Installation
You can install using pip:
```bash
python -m venv env
source env/bin/activate
pip install -r requirements.txt
```


