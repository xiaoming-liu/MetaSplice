"""Write the licence metadata of every installed Python package (run at image build)."""
from importlib.metadata import distributions

rows = []
for d in distributions():
    m = d.metadata
    lic = (m.get("License-Expression") or m.get("License") or "").strip().splitlines()
    cls = [c.split(" :: ")[-1] for c in (m.get_all("Classifier") or []) if c.startswith("License ::")]
    rows.append(f"{m['Name']}=={d.version}\t{(lic[0] if lic else '') or '; '.join(cls) or 'see package'}")
with open("/opt/metasplice/PYTHON_PACKAGE_LICENSES.txt", "w") as f:
    f.write("# Licence metadata of every installed Python package (generated at image build)\n")
    f.write("\n".join(sorted(rows, key=str.lower)) + "\n")
