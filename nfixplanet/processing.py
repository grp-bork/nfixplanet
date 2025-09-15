import logging
import pandas as pd

logger = logging.getLogger(__name__)


def load_hmm_output(path: str) -> pd.DataFrame:
    column_names = [
        "target_name",
        "accession",
        "query_name",
        "full_accession",
        "full_evalue",
        "full_score",
        "full_bias",
        "domain_evalue",
        "domain_score",
        "domain_bias",
        "exp",
        "reg",
        "clu",
        "ov",
        "env",
        "dom",
        "rep",
        "inc",
    ]
    columns_to_drop = [
        "accession",
        "full_accession",
    ]
    df = pd.read_csv(
        path,
        delim_whitespace=True,
        names=column_names,
        usecols=range(len(column_names)),
        comment="#",
    )
    df = df.drop(columns=columns_to_drop)
    ordered_cols = ["query_name"] + [
        col
        for col in column_names
        if col not in columns_to_drop and col != "query_name"
    ]
    df = df[ordered_cols]
    return df


def get_best_hits(path: str) -> pd.DataFrame:
    df = load_hmm_output(path)

    result = df.loc[
        df.groupby("query_name").apply(
            lambda g: g.sort_values(
                ["full_evalue", "full_score"], ascending=[True, False]
            ).index[0]
        )
    ].reset_index(drop=True)

    logger.debug(result)

    return result
