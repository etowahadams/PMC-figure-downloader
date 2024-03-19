from Bio import Entrez
import requests
import xml.etree.ElementTree as ET

headers = {
    'sec-ch-ua': '"Chromium";v="122", "Not(A:Brand";v="24", "Google Chrome";v="122"',
    'Referer': 'https://www.ncbi.nlm.nih.gov/pmc/articles/',
    'sec-ch-ua-mobile': '?0',
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36',
    'sec-ch-ua-platform': '"macOS"',
}

def download_pmc_papers(query: str, retmax=10, email="your_email@example.com") -> list[list[str, bytes]]:
    """
    Downloads the XML content of PMC papers matching the given query.

    Args:
        query (str): The search query for PMC papers.
        retmax (int, optional): The maximum number of papers to retrieve (default: 10).
        email (str, optional): The email address for Entrez (required for using the API).

    Returns:
        list: A list of strings containing the XML content of the retrieved papers.
    """
    Entrez.email = email
    handle = Entrez.esearch(db="pmc", term=query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]

    image_data: list[list[str, bytes]] = []
    for pmcid in ids[:1]:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:{pmcid}&metadataPrefix=pmc"
        response = requests.get(url)
        # Get the XML root element
        root = ET.fromstring(response.text)
        # Find all the graphic elements in the XML
        graphics = root.findall(
            ".//graphic", namespaces={"": "https://jats.nlm.nih.gov/ns/archiving/1.3/"}
        )
        # Iterate over the graphic elements and download the images
        for graphic in graphics:
            xlink_href = graphic.attrib["{http://www.w3.org/1999/xlink}href"]
            image_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmcid}/bin/{xlink_href}.jpg"
            # Download the image
            response = requests.get(image_url, headers=headers)
            # Check if the request was successful (status code 200)
            if response.status_code == 200:
                image_data.append([pmcid, response.content])
                print("Image downloaded successfully!")
            else:
                print(f"Failed to download image. Status code: {response.status_code}")

    return image_data


if __name__ == "__main__":
    query = "covid-19"
    image_data = download_pmc_papers(query, email="your-email-here@gmail.com")

    for i, image in enumerate(image_data):
        with open(f"image_{i}.jpg", "wb") as f:
            f.write(image[1])
