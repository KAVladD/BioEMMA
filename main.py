
import os
import cobra
import json
import argparse  
from urllib.request import urlretrieve
import pathlib

from mapper import *
from maps import * 
from reader import *
from merger import *

OUTPUT = "output_maps/"

def url_template(pathway_id):
    return f"https://rest.kegg.jp/get/rn{pathway_id}/kgml"

if __name__ == "__main__":

    pathlib.Path(OUTPUT).mkdir(parents=True, exist_ok=True)

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--model",
        type=str
    )

    parser.add_argument(
        "--pathways",
        type=str,
        default=None
    )   
    
    args = parser.parse_args()

    merger = EscherMerger()
    merger.reset()

    maps = []
    for pathway in args.pathways:
        maps.append("/maps/rn" + pathway + ".xml")
        if not os.path.exists(maps[-1]):
            urlretrieve(url_template(pathway), maps[-1])

    for map in maps:

        m = KeggMap()
        with open(map, "r") as f:
            m.read_from_file(f)

        model = cobra.io.read_sbml_model(args.model)
        mapper = EscherMapper(model)

        escher_map = mapper.map_flux_on_kegg(m)

        merger.append(escher_map)

        map_output_path = OUTPUT + os.path.basename(map.replace("xml", "json"))

        with open(map_output_path, "w") as f:
            json.dump(escher_map, f)
            print(map_output_path, " is saved!")

    merged = merger.merge_maps()
    print(len(merged[1]["nodes"]))

    with open(OUTPUT + "merged.json", "w") as f:
        json.dump(merged, f)


    