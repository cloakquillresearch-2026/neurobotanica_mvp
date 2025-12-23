# Example: Run this after starting uvicorn on port 8000
import httpx

payload = {
    "monomer": {
        "name": "THC",
        "smiles": "CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1"
    },
    "linkage_type": "METHYLENE"
}

response = httpx.post("http://localhost:8000/api/dimers/predict/homodimer", json=payload)
print("Status code:", response.status_code)
print("Response:", response.json())
