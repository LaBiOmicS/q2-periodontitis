import os
import shutil
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
import jinja2
import qiime2

ASSET_DIR = os.path.join(os.path.dirname(__file__), "assets", "web")

COLOR_MAP = {
    "Red": "#E41A1C",
    "Orange": "#FF7F00",
    "Extended": "#984EA3",
    "Yellow": "#FFFF33",
    "Blue": "#377EB8",
    "Green": "#4DAF4A",
    "Purple": "#88419D",
    "Other": "#999999"
}

FDR_ALPHA = 0.05


def _json_for_script(obj):
    """Serializes obj to JSON that is safe to embed inside a <script> tag."""
    return json.dumps(obj).replace("</", "<\\/")


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR-adjusted q-values (same order as input)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    # Enforce monotonicity from the largest p-value downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty(n)
    qvals[order] = np.minimum(ranked, 1.0)
    return qvals


def _group_comparisons(indices, metadata):
    """Kruskal-Wallis / Mann-Whitney U tests of every index against every
    categorical metadata column, with Benjamini-Hochberg correction."""
    meta_df = metadata.to_dataframe()
    categorical = [col for col in meta_df.columns
                   if metadata.columns[col].type == "categorical"]
    merged = indices.join(meta_df[categorical], how="inner")

    stats_results = []
    for col in categorical:
        valid_merged = merged.dropna(subset=[col])
        groups_dict = {k: v for k, v in valid_merged.groupby(col)
                       if len(v) >= 2}
        if len(groups_dict) < 2:
            continue
        for metric in indices.columns:
            group_data = [g[metric].dropna().values
                          for g in groups_dict.values()]
            group_data = [d for d in group_data if len(d) >= 2]
            if len(group_data) < 2:
                continue
            try:
                if len(group_data) == 2:
                    stat, pval = mannwhitneyu(group_data[0], group_data[1],
                                              alternative="two-sided")
                    method = "Mann-Whitney U"
                else:
                    stat, pval = kruskal(*group_data)
                    method = "Kruskal-Wallis"
            except ValueError as e:
                warnings.warn(
                    f"Skipping test of '{metric}' by metadata column "
                    f"'{col}': {e}", UserWarning)
                continue

            if np.isnan(pval):
                warnings.warn(
                    f"Skipping test of '{metric}' by metadata column "
                    f"'{col}': p-value is NaN (e.g. all values identical).",
                    UserWarning)
                continue

            stats_results.append({
                "metric": metric,
                "column": col,
                "method": method,
                "statistic": float(stat),
                "p_value": float(pval)
            })

    qvals = benjamini_hochberg([r["p_value"] for r in stats_results])
    for result, q in zip(stats_results, qvals):
        result["q_value"] = float(q)
        result["significant"] = bool(q < FDR_ALPHA)
    return stats_results


def summarize(
    output_dir: str,
    complexes: pd.DataFrame,
    indices: pd.DataFrame,
    mapping: pd.DataFrame = None,
    metadata: qiime2.Metadata = None
) -> None:

    # 1. Prepare data directory in output_dir
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    # Save TSV files for user download
    complexes_path = os.path.join(data_dir, "complexes.tsv")
    indices_path = os.path.join(data_dir, "indices.tsv")

    complexes.to_csv(complexes_path, sep="\t", index_label="sample-id")
    indices.to_csv(indices_path, sep="\t", index_label="sample-id")

    # 2. Copy CSS & vendored JS assets (no CDN: the .qzv must work offline)
    for sub in ("css", "js"):
        shutil.copytree(os.path.join(ASSET_DIR, sub),
                        os.path.join(output_dir, sub), dirs_exist_ok=True)

    # 3. Tax resolution audit
    total_features = len(mapping) if mapping is not None else 0
    resolution_counts = {}
    if mapping is not None and "Resolution" in mapping.columns:
        resolution_counts = mapping["Resolution"].value_counts().to_dict()

    # 4. Prepare JSON data for Vega-Lite interactive plots
    complexes_long = []
    for sample_id, row in complexes.iterrows():
        for complex_name, val in row.items():
            complexes_long.append({
                "sample_id": str(sample_id),
                "complex": str(complex_name),
                "abundance": float(val)
            })

    indices_data = []
    for sample_id, row in indices.iterrows():
        item = {"sample_id": str(sample_id)}
        for col in indices.columns:
            item[col] = float(row[col])
        indices_data.append(item)

    # 5. Statistical tests if metadata provided
    stats_results = []
    if metadata is not None:
        stats_results = _group_comparisons(indices, metadata)

    # 6. Render Jinja2 template
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(ASSET_DIR),
        autoescape=jinja2.select_autoescape(["html"])
    )
    template = env.get_template("index.html")

    indices_preview = indices.head(20).to_html(classes="table table-striped",
                                               border=0)

    html_content = template.render(
        total_samples=len(complexes),
        total_features=total_features,
        resolution_counts=resolution_counts,
        stats_results=stats_results,
        fdr_alpha=FDR_ALPHA,
        indices_table_html=indices_preview,
        complexes_json=_json_for_script(complexes_long),
        indices_json=_json_for_script(indices_data),
        color_map_json=_json_for_script(COLOR_MAP)
    )

    with open(os.path.join(output_dir, "index.html"), "w",
              encoding="utf-8") as fh:
        fh.write(html_content)
