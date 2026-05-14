import argparse

from bioemma.workflow import build_many_outputs, build_outputs


def build(args) -> None:
    if not args.output and not args.output_dir:
        raise ValueError("Either --output or --output-dir is required.")

    pathways = args.pathway or []
    kgmls = args.kgml or []
    sources_count = len(pathways) + len(kgmls)

    if sources_count > 1:
        if args.output:
            raise ValueError("--output can only be used with a single map. Use --output-dir for multiple maps.")
        result = build_many_outputs(
            model=args.model,
            pathways=pathways or None,
            kgmls=kgmls or None,
            output_dir=args.output_dir,
            merge=not args.no_merge,
            database=args.database,
            scaling_factor=args.scaling_factor,
            axis_epsilon=args.axis_epsilon,
            remove_orphan_metabolites=args.remove_orphan_metabolites,
            run_fba=args.run_fba,
            save_html=args.save_html,
            save_png=args.save_png,
        )
        for name, path in result.paths.items():
            print(f"{name}: {path}")
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
        remove_orphan_metabolites=args.remove_orphan_metabolites,
        run_fba=args.run_fba,
        save_html=args.save_html,
        save_png=args.save_png,
    )

    for name, path in result.paths.items():
        print(f"{name}: {path}")


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
    build_parser.add_argument("--scaling-factor", type=float, default=4)
    build_parser.add_argument("--axis-epsilon", type=float, default=2)
    build_parser.add_argument("--remove-orphan-metabolites", action="store_true")
    build_parser.add_argument("--run-fba", action="store_true")
    build_parser.add_argument("--save-html", action="store_true")
    build_parser.add_argument("--save-png", action="store_true")
    build_parser.add_argument("--no-merge", action="store_true")
    build_parser.set_defaults(func=build)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
