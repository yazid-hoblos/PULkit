from bs4 import BeautifulSoup
import re

def extract_gene_data(html_content):
    """
    Extract gene information from HTML content.
    Returns a dictionary with position, CAZyme GH equivalent, link, and CGCFinder status.
    """
    soup = BeautifulSoup(html_content, 'html.parser')
    
    gene_data = {
        'position': None,
        'cazyme_gh': None,
        'link': None,
        'found_by_cgcfinder': None
    }
    
    # Extract position (highlighted number)
    # Look for the highlighted span or any number pattern
    position_element = soup.find('span', style=lambda x: x and 'background-color' in x)
    if position_element:
        gene_data['position'] = position_element.get_text().strip()
    else:
        # Alternative: look for numbers in parentheses or other patterns
        text_content = soup.get_text()
        position_match = re.search(r'(\d+)', text_content)
        if position_match:
            gene_data['position'] = position_match.group(1)
    
    # Extract CAZyme GH equivalent
    cazy_link = soup.find('a', href=lambda x: x and 'cazy.org' in x)
    if cazy_link:
        gene_data['cazyme_gh'] = cazy_link.get_text().strip()
    
    # Extract link (NCBI protein link)
    ncbi_link = soup.find('a', href=lambda x: x and 'ncbi.nlm.nih.gov' in x)
    if ncbi_link:
        gene_data['link'] = ncbi_link['href']
    
    # Extract CGCFinder status
    cgc_header = soup.find('h3', class_='ui left aligned header')
    if cgc_header:
        cgc_text = cgc_header.get_text().strip().lower()
        if 'found by cgcfinder' in cgc_text:
            gene_data['found_by_cgcfinder'] = 'Yes'
        else:
            gene_data['found_by_cgcfinder'] = 'No'
    else:
        # Alternative: look for "Yes" or "No" in table cells
        cells = soup.find_all('td')
        for cell in cells:
            cell_text = cell.get_text().strip().lower()
            if cell_text == 'yes':
                gene_data['found_by_cgcfinder'] = 'Yes'
                break
            elif cell_text == 'no':
                gene_data['found_by_cgcfinder'] = 'No'
                break
    
    return gene_data

def extract_multiple_genes(html_list):
    """
    Extract data from multiple genes.
    html_list: list of HTML content strings
    Returns: list of dictionaries with gene data
    """
    genes_data = []
    
    for i, html_content in enumerate(html_list):
        gene_info = extract_gene_data(html_content)
        gene_info['gene_index'] = i + 1  # Add gene index for reference
        genes_data.append(gene_info)
    
    return genes_data

# Alternative function for more robust extraction
def extract_gene_data_robust(html_content):
    """
    More robust extraction that handles various HTML structures.
    """
    soup = BeautifulSoup(html_content, 'html.parser')
    
    gene_data = {
        'position': None,
        'cazyme_gh': None,
        'link': None,
        'found_by_cgcfinder': None
    }
    
    # Extract all table rows for structured parsing
    rows = soup.find_all('tr')
    
    for row in rows:
        cells = row.find_all(['td', 'th'])
        if len(cells) >= 2:
            key = cells[0].get_text().strip().lower()
            value_cell = cells[1]
            
            # Check if this row contains position info
            if any(keyword in key for keyword in ['position', 'gene', 'locus']):
                # Look for highlighted numbers or just numbers
                highlight = value_cell.find('span', style=lambda x: x and 'background-color' in x)
                if highlight:
                    gene_data['position'] = highlight.get_text().strip()
                else:
                    # Extract numbers from the cell
                    numbers = re.findall(r'\d+', value_cell.get_text())
                    if numbers:
                        gene_data['position'] = numbers[0]
            
            # Check for CAZyme info
            if 'cazyme' in key:
                cazy_link = value_cell.find('a', href=lambda x: x and 'cazy.org' in x)
                if cazy_link:
                    gene_data['cazyme_gh'] = cazy_link.get_text().strip()
            
            # Check for CGCFinder status
            if 'cgc' in key or 'finder' in key:
                cgc_text = value_cell.get_text().strip().lower()
                if 'yes' in cgc_text:
                    gene_data['found_by_cgcfinder'] = 'Yes'
                elif 'no' in cgc_text:
                    gene_data['found_by_cgcfinder'] = 'No'
    
    # Extract links (prioritize NCBI)
    all_links = soup.find_all('a', href=True)
    for link in all_links:
        href = link['href']
        if 'ncbi.nlm.nih.gov' in href:
            gene_data['link'] = href
            break
    
    # If no NCBI link, take the first meaningful link
    if not gene_data['link'] and all_links:
        for link in all_links:
            href = link['href']
            if not href.startswith('#') and 'cazy.org' not in href:
                gene_data['link'] = href
                break
    
    return gene_data

def read_html_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.read()
    
file= read_html_file("y.txt")  # Replace with your actual HTML file path
result = extract_gene_data_robust(file)

for key, value in result.items():
    print(f"{key}: {value}")