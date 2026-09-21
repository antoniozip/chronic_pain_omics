"""
Check PRIDE file download mechanisms and API structure.
"""
import json
import ssl
import urllib.request

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def api_get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, context=ctx, timeout=20) as resp:
        return json.loads(resp.read().decode())

# Check specific file details for PXD013362
acc = "PXD013362"

# Check the full API response structure
data = api_get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}")
print(f"Top-level keys: {list(data.keys())}")
print(f"'_links' keys: {list(data.get('_links', {}).keys())}")
print()

# Check files with full details
files_data = api_get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files")
file_list = (files_data if isinstance(files_data, list)
             else files_data.get("_embedded", {}).get("files", []))

# Show structure of a file entry
if file_list:
    f = file_list[0]
    print(f"File entry keys: {list(f.keys())}")
    print(f"File entry sample: {json.dumps(f, indent=2)[:2000]}")

# Check if there's a download URL  
print("\n\nFile download URLs:")
for f in file_list[:5]:
    if isinstance(f, dict):
        fname = f.get("fileName", "")
        public_uri = f.get("publicFileUri", f.get("downloadLink", "N/A"))
        ftp_url = f.get("ftpUrl", "N/A")
        print(f"  {fname}:")
        print(f"    publicFileUri: {public_uri}")
        print(f"    ftpUrl: {ftp_url}")
        print(f"    _links: {f.get('_links', {})}")
        print()

# Check the PRIDE FTP
print("\n\nPRIDE FTP base:")
print("ftp://ftp.pride.ebi.ac.uk/pride/data/archive/")

# Check if there's a way to get quantified protein groups
print("\n\nChecking for quant files in PXD013362:")
for f in file_list:
    if isinstance(f, dict):
        fname = f.get("fileName", "")
        if any(x in fname for x in ["quantitation", "Quantitation", "peptide_quant", "protein"]):
            ft = f.get("fileType", {})
            if isinstance(ft, dict):
                ftype = ft.get("name", "")
            else:
                ftype = str(ft)
            print(f"  {fname}")
            print(f"    type: {ftype}")
            print(f"    size: {f.get('fileSizeBytes', 0)}")
            print(f"    download: {f.get('publicFileUri', 'N/A')}")
            print()
