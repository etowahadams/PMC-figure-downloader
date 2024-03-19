from Bio import Entrez
import requests
import xml.etree.ElementTree as ET
import codecs
import polars as pl
import time


def fix_cp1252_to_utf8(text):
    """
    Converts text that was originally encoded in CP-1252 but was read and sent as UTF-8,
    back to its original form. This function was written by Claude.

    Args:
        text (str): The input text that needs to be fixed.

    Returns:
        str: The fixed text, encoded in UTF-8.
    """
    try:
        # Try decoding the text as UTF-8 first
        decoded_text = text.encode("latin-1").decode("utf-8")
    except UnicodeDecodeError:
        # If decoding as UTF-8 fails, try decoding as CP-1252
        decoded_text = codecs.decode(text, "cp1252")

    # Re-encode the decoded text as UTF-8
    fixed_text = decoded_text.encode("utf-8", errors="replace")

    return fixed_text.decode("utf-8")


ns = {"": "https://jats.nlm.nih.gov/ns/archiving/1.3/"}


def get_text_or_empty(element, xpath):
    node = element.find(xpath, namespaces=ns)
    text_content = ET.tostring(node, method="text", encoding="unicode")
    return text_content if node is not None else ""

# TODO: Determine which of these are necessary
headers = {
    "sec-ch-ua": '"Chromium";v="122", "Not(A:Brand";v="24", "Google Chrome";v="122"',
    "Referer": "https://www.ncbi.nlm.nih.gov/pmc/articles/",
    "sec-ch-ua-mobile": "?0",
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36",
    "sec-ch-ua-platform": '"macOS"',
}


def extract_pmc_figures(ids: list[str]) -> pl.DataFrame:
    """
    Downloads the papers from PMC and extracts the figure data from the XML.

    Args:
        ids (list[str]): A list of PMCID's.

    Returns:
        pl.DataFrame: A DataFrame containing the figure data.
    """

    data = {
        "pmcid": [],
        "fig_id": [],
        "fig_label": [],
        "fig_title": [],
        "fig_desc": [],
        "image_url": [],
    }

    for i, pmcid in enumerate(ids):
        time.sleep(0.2) # Be nice to the server
        print(f"Getting paper {i + 1} of {len(ids)}...")

        url = f"https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:{pmcid}&metadataPrefix=pmc"
        response = requests.get(url)
        # Fix the encoding of the response text
        response_text = fix_cp1252_to_utf8(response.text)
        # Get the XML root element
        root = ET.fromstring(response_text)
        # Find all the graphic elements in the XML
        figures = root.findall(".//fig", namespaces=ns)
        for figure in figures:
            caption = get_text_or_empty(figure, ".//caption/p")
            title = get_text_or_empty(figure, ".//caption/title")
            label = get_text_or_empty(figure, ".//label")
            id = figure.attrib["id"]

            graphic = figure.find(
                ".//graphic",
                namespaces=ns,
            )
            if graphic is not None:
                xlink_href = graphic.attrib["{http://www.w3.org/1999/xlink}href"]
                image_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmcid}/bin/{xlink_href}.jpg"

                data["pmcid"].append(pmcid)
                data["fig_label"].append(label)
                data["fig_desc"].append(caption)
                data["fig_title"].append(title)
                data["image_url"].append(image_url)
                data["fig_id"].append(id)

            else:
                print("No graphic element found in the XML")

    return pl.DataFrame(data)


def download_imgs(data: pl.DataFrame, dir: str) -> pl.DataFrame:
    status = {
        "pmcid": [],
        "fig_id": [],
        "status": [],
    }
    for i, url in enumerate(data["image_url"]):
        print(f"Downloading image {i + 1} of {len(data)}...")
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            with open(f"{dir}/{data['pmcid'][i]}_{data['fig_id'][i]}.jpg", "wb") as f:
                f.write(response.content)
        else:
            print(f"Failed to download image. Status code: {response.status_code}")
        status["pmcid"].append(data["pmcid"][i])
        status["fig_id"].append(data["fig_id"][i])
        status["status"].append(response.status_code)

    return pl.DataFrame(status)


def search_pmc(query: str, email: str, max_results: int = None, retmax: int = 100):
    """
    Searches PMC for papers matching the given query.
    Returns a list of PMCID's.
    """
    if max_results and max_results < retmax:
        retmax = max_results
    
    Entrez.email = email
    handle = Entrez.esearch(db="pmc", term=query, retmax=retmax)
    record = Entrez.read(handle)
    print("Found", record["Count"], "papers. Getting their IDs...")
    ids = []
    for i in range(int(record["Count"]) // retmax + 1):
        handle = Entrez.esearch(
            db="pmc", term=query, retmax=retmax, retstart=i * retmax
        )
        record = Entrez.read(handle)
        ids.extend(record["IdList"])
        # If we've reached the maximum number of results, stop getting more IDs
        if max_results and (i + 1) * retmax >= max_results:
            break
        time.sleep(0.2)
    print("Finished getting", len(ids), "IDs.")

    return ids


if __name__ == "__main__":
    query = '"Nature Genetics"[Journal] AND "open access"[filter]'
    ids = search_pmc(query, email="etowahadams@gmail.com", max_results=10, retmax=10)
    figure_data = extract_pmc_figures(ids)
    print(figure_data)
    figure_data.write_parquet("figures.parquet")
    # data = download_pmc_papers(query, email="your-email-here@gmail.com")
    # download_imgs(data)
