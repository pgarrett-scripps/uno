import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser(description="Optional CNN training for Uno")
    commands = parser.add_subparsers(dest="command", required=True)
    train = commands.add_parser("train", help="adapt the bundled CNN or train from scratch")
    train.add_argument("--data", required=True)
    train.add_argument("--out", required=True)
    train.add_argument("--units", choices=["minutes", "seconds", "irt", "hi"], required=True)
    train.add_argument("--method", required=True)
    starting = train.add_mutually_exclusive_group()
    starting.add_argument("--init")
    starting.add_argument("--from-scratch", action="store_true")
    validation = train.add_mutually_exclusive_group()
    validation.add_argument("--validation")
    validation.add_argument("--holdout-run")
    train.add_argument("--seed", type=int, default=0)
    train.add_argument("--epochs", type=int, default=40)
    train.add_argument("--config")
    args = vars(parser.parse_args())
    args.pop("command")
    try:
        from .training import train as run
        print(json.dumps(run(**args), allow_nan=False, indent=2))
    except (ValueError, OSError, RuntimeError, ImportError) as error:
        parser.exit(1, f"RT training failed: {error}\n")


if __name__ == "__main__":
    main()
