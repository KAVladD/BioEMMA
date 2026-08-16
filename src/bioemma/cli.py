import argparse
from typing import Any

from bioemma.workflow import build_many_outputs, build_outputs


def _print_map_stats(summary: dict[str, Any], label: str | None = None) -> None:
    map_stats = summary.get("map_stats", {})
    stages = map_stats.get("stages", [])
    if not stages:
        return

    header = "map_stats" if label is None else f"map_stats ({label})"
    print(f"{header}:")
    for stage in stages:
        counts = stage["counts"]
        total_change = stage["change"]["total_elements"]
        nodes_change = stage["change"]["nodes"]
        reactions_change = stage["change"]["reactions"]
        segments_change = stage["change"]["segments"]
        print(
            "  "
            f"{stage['name']}: "
            f"total={counts['total_elements']} "
            f"(+{total_change['added']}/-{total_change['removed']}), "
            f"nodes={counts['nodes']} "
            f"(+{nodes_change['added']}/-{nodes_change['removed']}), "
            f"reactions={counts['reactions']} "
            f"(+{reactions_change['added']}/-{reactions_change['removed']}), "
            f"segments={counts['segments']} "
            f"(+{segments_change['added']}/-{segments_change['removed']})"
        )


def build(args) -> None:
    if not args.output and not args.output_dir:
        raise ValueError("Either --output or --output-dir is required.")

    pathways = args.pathway or []
    kgmls = args.kgml or []
    sources_count = len(pathways) + len(kgmls)

    if sources_count > 1:
        if args.output:
            raise ValueError(
                "--output can only be used with a single map. Use --output-dir for multiple maps."
            )
        result = build_many_outputs(
            model=args.model,
            pathways=pathways or None,
            kgmls=kgmls or None,
            output_dir=args.output_dir,
            merge=not args.no_merge,
            database=args.database,
            scaling_factor=args.scaling_factor,
            axis_epsilon=args.axis_epsilon,
            markers_dist=args.markers_dist,
            metabolite_label_shift=args.metabolite_label_shift,
            reaction_label_shift=args.reaction_label_shift,
            canvas_margin_x=args.canvas_margin_x,
            canvas_margin_y=args.canvas_margin_y,
            multimarker_distance_fraction=args.multimarker_distance_fraction,
            use_constant_multimarker_distance=args.use_constant_multimarker_distance,
            constant_multimarker_distance=args.constant_multimarker_distance,
            axis_offset=args.axis_offset,
            secondary_metabolite_distance=args.secondary_metabolite_distance,
            secondary_metabolite_spacing=args.secondary_metabolite_spacing,
            remove_orphan_metabolites=args.remove_orphan_metabolites,
            include_kegg_only=args.include_kegg_only,
            use_model_metabolite_ids=args.use_model_metabolite_ids,
            use_database_secondary_metabolite_ids=(
                args.use_database_secondary_metabolite_ids
            ),
            metabolite_id_compartments=args.metabolite_id_compartments,
            compartment=args.compartment,
            use_fallback_matching=args.fallback_matching,
            use_bigg_fallback_matching=args.bigg_fallback_matching,
            use_seed_fallback_matching=args.seed_fallback_matching,
            save_kegg_map=args.save_kegg_map,
            run_fba=args.run_fba,
            save_html=args.save_html,
        )
        for name, path in result.paths.items():
            print(f"{name}: {path}")
        if args.map_stats:
            for item_result in result.results:
                label = item_result.paths.get("escher_map_json")
                label = label.parent.name if label else None
                _print_map_stats(item_result.summary, label=label)
        return

    result = build_outputs(
        model=args.model,
        pathway=pathways[0] if pathways else None,
        kgml=kgmls[0] if kgmls else None,
        output_dir=args.output_dir,
        map_json_path=args.output,
        database=args.database,
        scaling_factor=args.scaling_factor,
        axis_epsilon=args.axis_epsilon,
        markers_dist=args.markers_dist,
        metabolite_label_shift=args.metabolite_label_shift,
        reaction_label_shift=args.reaction_label_shift,
        canvas_margin_x=args.canvas_margin_x,
        canvas_margin_y=args.canvas_margin_y,
        multimarker_distance_fraction=args.multimarker_distance_fraction,
        use_constant_multimarker_distance=args.use_constant_multimarker_distance,
        constant_multimarker_distance=args.constant_multimarker_distance,
        axis_offset=args.axis_offset,
        secondary_metabolite_distance=args.secondary_metabolite_distance,
        secondary_metabolite_spacing=args.secondary_metabolite_spacing,
        remove_orphan_metabolites=args.remove_orphan_metabolites,
        include_kegg_only=args.include_kegg_only,
        use_model_metabolite_ids=args.use_model_metabolite_ids,
        use_database_secondary_metabolite_ids=(
            args.use_database_secondary_metabolite_ids
        ),
        metabolite_id_compartments=args.metabolite_id_compartments,
        compartment=args.compartment,
        use_fallback_matching=args.fallback_matching,
        use_bigg_fallback_matching=args.bigg_fallback_matching,
        use_seed_fallback_matching=args.seed_fallback_matching,
        save_kegg_map=args.save_kegg_map,
        run_fba=args.run_fba,
        save_html=args.save_html,
    )

    for name, path in result.paths.items():
        print(f"{name}: {path}")
    if args.map_stats:
        _print_map_stats(result.summary)


def main() -> None:
    parser = argparse.ArgumentParser(prog="bioemma")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--model", required=True)
    source_group = build_parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--kgml", nargs="+")
    source_group.add_argument("--pathway", nargs="+")
    build_parser.add_argument("--output")
    build_parser.add_argument("--output-dir")
    build_parser.add_argument("--database", choices=["BIGG", "SEED", "KEGG"], default="BIGG")
    build_parser.add_argument("--scaling-factor", type=float)
    build_parser.add_argument("--axis-epsilon", type=float)
    build_parser.add_argument("--markers-dist", type=float)
    build_parser.add_argument("--metabolite-label-shift", nargs=2, type=float, metavar=("DX", "DY"))
    build_parser.add_argument("--reaction-label-shift", nargs=2, type=float, metavar=("DX", "DY"))
    build_parser.add_argument("--canvas-margin-x", type=float)
    build_parser.add_argument("--canvas-margin-y", type=float)
    build_parser.add_argument("--multimarker-distance-fraction", type=float)
    build_parser.add_argument(
        "--use-constant-multimarker-distance",
        action="store_true",
        default=None,
    )
    build_parser.add_argument("--constant-multimarker-distance", type=float)
    build_parser.add_argument("--axis-offset", type=float)
    build_parser.add_argument("--secondary-metabolite-distance", type=float)
    build_parser.add_argument("--secondary-metabolite-spacing", type=float)
    build_parser.add_argument(
        "--remove-orphan-metabolites",
        "--remove-free-metabolites",
        action="store_true",
        dest="remove_orphan_metabolites",
    )
    build_parser.add_argument("--include-kegg-only", action="store_true")
    build_parser.add_argument("--use-model-metabolite-ids", action="store_true")
    build_parser.add_argument(
        "--use-database-secondary-metabolite-ids",
        action="store_true",
        help=(
            "Use identifiers from the selected database for secondary metabolites "
            "instead of raw COBRA model metabolite IDs. Ignored when "
            "--use-model-metabolite-ids is set."
        ),
    )
    build_parser.add_argument(
        "--metabolite-id-compartments",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Include model compartments in metabolite output where supported. "
            "Defaults to on for BIGG model IDs and off otherwise."
        ),
    )
    build_parser.add_argument(
        "--compartment",
        default=None,
        help=(
            "Keep only model-matched reactions whose COBRA reaction uses this "
            "compartment ID, for example 'm' for mitochondrion."
        ),
    )
    build_parser.add_argument(
        "--fallback-matching",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Use fallback reaction matching through EC-derived BiGG/SEED, "
            "MetaCyc, Rhea, and EC annotations. Disabled by default to keep "
            "strict KEGG/BiGG/SEED matches only."
        ),
    )
    build_parser.add_argument(
        "--bigg-fallback-matching",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Use BiGG identifiers inferred by EC fallback. Requires "
            "--fallback-matching to be enabled."
        ),
    )
    build_parser.add_argument(
        "--seed-fallback-matching",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Use SEED identifiers inferred by EC fallback. Requires "
            "--fallback-matching to be enabled."
        ),
    )
    build_parser.add_argument("--save-kegg-map", action="store_true")
    build_parser.add_argument("--run-fba", action="store_true")
    build_parser.add_argument("--save-html", action="store_true")
    build_parser.add_argument("--no-merge", action="store_true")
    build_parser.add_argument("--map-stats", action="store_true")
    build_parser.set_defaults(func=build)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
