import xml.etree.ElementTree as ET
NS = "{tram:components}"
r = ET.parse("data/v2/rigidity/wt_components.xml").getroot()
 
# THE DICTIONARY: graph/components (NOT graph/states/state)
comps_dict = r.find(f"{NS}components")
print("found top-level <components> dictionary:", comps_dict is not None)
if comps_dict is not None:
    entries = comps_dict.findall(f"{NS}component")
    id_size = {}
    for c in entries:
        cid = c.get("id"); sz = c.get("size")
        if cid is not None and sz is not None:
            id_size[int(cid)] = int(sz)
    sizes = sorted(id_size.values(), reverse=True)
    total = sum(sizes)
    print(f"dictionary entries: {len(entries)}")
    print(f"entries with id+size: {len(id_size)}")
    print(f"total atoms in dictionary: {total}")
    print(f"top 8 cluster sizes: {sizes[:8]}")
    print(f"largest single cluster: {sizes[0]} ({100*sizes[0]/total:.1f}% of total)" if total else "no total")