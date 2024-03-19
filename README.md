# PMC Figure Downloader

This is a simple script to download figures from open access papers in [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/) database. It uses the [PMC API](https://www.ncbi.nlm.nih.gov/pmc/tools/developers/) to search for articles and uses the 
[PMC Open Access Web Service API](https://www.ncbi.nlm.nih.gov/pmc/tools/oa-service/) to get the XML of each paper to determine the figure information and URLs.

## Installation
You can install using pip:
```bash
python -m venv env
source env/bin/activate
pip install -r requirements.txt
```

## Usage
See `example.ipynb` for a Jupyter notebook example.

### Get paper IDs from a search query
```python
query = '"Nature Genetics"[Journal] AND "open access"[filter]'
# a list of paper IDs
result_ids = search_pmc(query, email="your-email-here@gmail.com", max_results=10)
```
### Extract figure information and URLs from a paper
```python
figure_data = extract_pmc_figures(result_ids)
figure_data.write_parquet("figure_data.parquet")
```
The `figure_data` dataframe looks like this:
```
┌──────────┬────────┬───────────┬──────────────────────┬─────────────────────┬─────────────────────┐
│ pmcid    ┆ fig_id ┆ fig_label ┆ fig_title            ┆ fig_desc            ┆ image_url           │
│ ---      ┆ ---    ┆ ---       ┆ ---                  ┆ ---                 ┆ ---                 │
│ str      ┆ str    ┆ str       ┆ str                  ┆ str                 ┆ str                 │
╞══════════╪════════╪═══════════╪══════════════════════╪═════════════════════╪═════════════════════╡
│ 10937393 ┆ Fig1   ┆ Fig. 1    ┆ FANS-based isolation ┆ a, Schematic        ┆ https://www.ncbi.nl │
│          ┆        ┆           ┆ of nuclei o…         ┆ representation of   ┆ m.nih.gov/pmc…      │
│          ┆        ┆           ┆                      ┆ t…                  ┆                     │
│ 10937393 ┆ Fig2   ┆ Fig. 2    ┆ Purity and           ┆ a, Heatmaps depict  ┆ https://www.ncbi.nl │
│          ┆        ┆           ┆ reproducibility of   ┆ log2-transfor…      ┆ m.nih.gov/pmc…      │
│          ┆        ┆           ┆ th…                  ┆                     ┆                     │
```

### Download figures to a directory
```python
output_dir = "img"
download_status = download_imgs(figure_data, output_dir)
print("Failed downloads:")
download_status.filter(pl.col("status") != 200)
```
The `download_status` dataframe will look like this:
```
┌──────────┬────────┬────────┐
│ pmcid    ┆ fig_id ┆ status │
│ ---      ┆ ---    ┆ ---    │
│ str      ┆ str    ┆ i64    │
╞══════════╪════════╪════════╡
│ 10937393 ┆ Fig1   ┆ 200    │
│ 10937393 ┆ Fig2   ┆ 200    │
```