import argparse
import json
from urllib.request import urlopen

import cobra

from bioemma.maps import KeggMap
from bioemma.mapper_base import EscherMapper


def pathway_url(pathway_id: str) -> str:
    pathway_id = pathway_id.removeprefix("rn")
    return f"https://rest.kegg.jp/get/rn{pathway_id}/kgml"


def load_kegg_map(args) -> KeggMap:
    kegg_map = KeggMap()

    if args.kgml:
        with open(args.kgml, "r", encoding="utf-8") as file:
            kegg_map.read_from_file(file)
        return kegg_map

    if args.pathway:
        with urlopen(pathway_url(args.pathway)) as response:
            kgml = response.read().decode("utf-8")
        kegg_map.read_from_file(kgml)
        return kegg_map

    raise ValueError("Either --kgml or --pathway is required.")


def build(args) -> None:
    kegg_map = load_kegg_map(args)
    model = cobra.io.read_sbml_model(args.model)

    mapper = EscherMapper(
        metabolites=kegg_map.get_metabolites(),
        reactions=kegg_map.get_reactions(),
        database=args.database,
        scaling_factor=args.scaling_factor,
        axis_epsilon=args.axis_epsilon,
    )

    escher_map = mapper.build_map(model)

    with open(args.output, "w", encoding="utf-8") as file:
        json.dump(escher_map, file)


def main() -> None:
    parser = argparse.ArgumentParser(prog="bioemma")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--model", required=True)
    source_group = build_parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--kgml")
    source_group.add_argument("--pathway")
    build_parser.add_argument("--output", required=True)
    build_parser.add_argument("--database", choices=["BIGG", "SEED", "KEGG"], default="BIGG")
    build_parser.add_argument("--scaling-factor", type=float, default=4)
    build_parser.add_argument("--axis-epsilon", type=float, default=2)
    build_parser.set_defaults(func=build)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
