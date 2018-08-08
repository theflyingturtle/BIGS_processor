#!/usr/bin/env python3

import click
import json
import logging
import pandas as pd

@click.command()
@click.argument("json_files", nargs=-1, type=click.File())
def main(json_files):
    """Mangle BIGS dumps into a more sensible csv layout."""
    for f in json_files:
        j = json.load(f)

        df = {(int(n), level): prediction for n, v in j.items() for level, prediction in enumerate(v['analysis']['taxon_prediction'])}
        df = pd.DataFrame.from_dict(df, orient='index').unstack()
        df['name'] = df.index.map({int(n): v['isolate'] for n, v in j.items()}.get)

        df.to_csv(f.name + ".csv")
        logging.info("Produced %s.csv", f.name)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
