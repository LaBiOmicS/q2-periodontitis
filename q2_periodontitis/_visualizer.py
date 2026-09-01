import os
import shutil
import json
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

    # 2. Copy CSS & web assets
    css_dir = os.path.join(output_dir, "css")
    os.makedirs(css_dir, exist_ok=True)
    shutil.copy(os.path.join(ASSET_DIR, "css", "style.css"), os.path.join(css_dir, "style.css"))

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
        meta_df = metadata.to_dataframe()
        merged = indices.join(meta_df, how="inner")

        for col in meta_df.columns:
            # Check categorical columns (between 2 and 20 unique values)
            valid_merged = merged.dropna(subset=[col])
            n_unique = valid_merged[col].nunique()
            if 2 <= n_unique <= 20:
                groups_dict = {k: v for k, v in valid_merged.groupby(col) if len(v) >= 2}
                if len(groups_dict) >= 2:
                    for metric in ["pMDI", "pMDI_clr", "pMDI_ilr", "SWL", "KRS"]:
                        group_data = [g[metric].dropna().values for g in groups_dict.values()]
                        group_data = [d for d in group_data if len(d) >= 2]
                        if len(group_data) >= 2:
                            try:
                                if len(group_data) == 2:
                                    stat, pval = mannwhitneyu(group_data[0], group_data[1])
                                    method = "Mann-Whitney U"
                                else:
                                    stat, pval = kruskal(*group_data)
                                    method = "Kruskal-Wallis"

                                if not np.isnan(pval):
                                    stats_results.append({
                                        "metric": metric,
                                        "column": col,
                                        "method": method,
                                        "statistic": float(stat),
                                        "p_value": float(pval)
                                    })
                            except Exception:
                                pass

    # 6. Render Jinja2 template
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(ASSET_DIR))
    template = env.get_template("index.html")

    indices_preview = indices.head(20).to_html(classes="table table-striped", border=0)

    html_content = template.render(
        total_samples=len(complexes),
        total_features=total_features,
        resolution_counts=resolution_counts,
        stats_results=stats_results,
        indices_table_html=indices_preview,
        complexes_json=json.dumps(complexes_long),
        indices_json=json.dumps(indices_data),
        color_map_json=json.dumps(COLOR_MAP)
    )

    with open(os.path.join(output_dir, "index.html"), "w", encoding="utf-8") as fh:
        fh.write(html_content)
