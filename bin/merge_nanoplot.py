#!/usr/bin/env python3
import os
import glob
import argparse
import re

# Ignore html if they match one of these patterns
SKIP_PATTERNS = [
    "nanoplot-report.html",
    "weighted",
    "lengthvsqualityscatterplot_kde",
]

# Improve plot titles
RENAME_PLOTS = {
    "LengthvsQualityScatterPlot_dot.html": "Length_vs_Quality_Scatterplot",
    "Yield_By_Length.html": "Yield_by_Length",
}

# CSS
CSS = """
<style>
#nanoplot-merged-report {
    font-family: Arial, sans-serif;
    padding: 20px;
    color: #222;
}

html[data-bs-theme="dark"] #nanoplot-merged-report {
    color: #ddd;
}

#nanoplot-merged-report .container {
    max-width: 1200px;
    margin: 0 auto;
}

#nanoplot-merged-report .plot-section {
    background: #ffffff;
    border-radius: 10px;
    padding: 15px;
    margin-bottom: 25px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.08);
}

html[data-bs-theme="dark"] #nanoplot-merged-report .plot-section {
    background: #1e1e1e;
}

#nanoplot-merged-report h2 {
    font-size: 18px;
    margin-bottom: 10px;
    border-bottom: 1px solid #eee;
    padding-bottom: 6px;
}

html[data-bs-theme="dark"] #nanoplot-merged-report h2 {
    border-color: #333;
}

#nanoplot-merged-report .button-row {
    margin: 10px 0;
}

#nanoplot-merged-report button {
    margin: 3px;
    padding: 6px 10px;
    font-size: 12px;
    border-radius: 6px;
    cursor: pointer;
    border: 1px solid #cfd6e2;
    background: #f8fafc;
}

#nanoplot-merged-report button.active {
    background: #1f77b4;
    color: #fff;
    border-color: #1f77b4;
}

html[data-bs-theme="dark"] #nanoplot-merged-report button {
    background: #2b2b2b;
    color: #ddd;
    border-color: #444;
}

#nanoplot-merged-report .plot-wrapper {
    display: none;
    margin-top: 15px;
    padding: 10px;
    background: #f7f9fc;
    border-radius: 8px;
}

html[data-bs-theme="dark"] #nanoplot-merged-report .plot-wrapper {
    background: #2a2a2a;
}

#nanoplot-merged-report .plot-wrapper > div {
    width: 100% !important;
    height: 720px !important;
}
</style>
"""

# HTML template 
HTML_TEMPLATE = """<!--
id: 'nanoplot_plots'
section_name: 'NanoPlot_plots'
-->
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
{css}
<script src="https://cdn.plot.ly/plotly-2.16.1.min.js"></script>
</head>
<body>
<div id="nanoplot-merged-report">
  <div class="container">
    {body}
  </div>
</div>

<script>
function showSample(plotId, sampleId) {{
    const section = document.getElementById(plotId);

    section.querySelectorAll(".plot-wrapper")
        .forEach(el => el.style.display = "none");

    section.querySelectorAll("button")
        .forEach(btn => btn.classList.remove("active"));

    const target = document.getElementById(plotId + "_" + sampleId);
    if (!target) return;

    target.style.display = "block";

    const activeBtn = section.querySelector(
        'button[data-sample="' + sampleId + '"]'
    );
    if (activeBtn) activeBtn.classList.add("active");

    const graph = target.querySelector(".plotly-graph-div");
    if (graph) setTimeout(() => Plotly.Plots.resize(graph), 50);
}}

window.addEventListener("DOMContentLoaded", () => {{
    document.querySelectorAll(".plot-section").forEach(section => {{
        const firstBtn = section.querySelector("button");
        if (firstBtn) {{
            showSample(section.id, firstBtn.dataset.sample);
        }}
    }});
}});
</script>
</body>
</html>
"""

# Remove spaces with underscore
def changes_id(text):
    return re.sub(r"[^A-Za-z0-9]", "_", text)

# Obtain sample id from directory name
def get_sample_id(dir_path):
    return re.sub(r"_nanoplot$", "", os.path.basename(dir_path))

# Remove sample id from file name
def remove_sample_prefix(filename, sample_id):
    base = os.path.basename(filename)
    if base.startswith(sample_id + "_"):
        return base[len(sample_id) + 1 :]
    return base

def main(input_dirs, output_file):
    plot_map = {}

    for nanoplot_dir in input_dirs:
        sample_id = get_sample_id(nanoplot_dir)

        html_files = [
            f for f in glob.glob(os.path.join(nanoplot_dir, "*.html"))
            if not any(pat in os.path.basename(f).lower() for pat in SKIP_PATTERNS)
        ]

        for html_file in html_files:
            plot_name = remove_sample_prefix(html_file, sample_id)
            plot_name = RENAME_PLOTS.get(plot_name, plot_name)

            with open(html_file, "r", encoding="utf-8") as fh:
                html = fh.read()

            plot_map.setdefault(plot_name, {})[sample_id] = html

    body = ""

    for plot_name, samples in plot_map.items():
        plot_id = change_id(plot_name)
        body += f"<div class='plot-section' id='{plot_id}'>\n"
        body += f"<h2>{plot_name}</h2>\n<div class='button-row'>\n"

        for sample_id in samples:
            body += (
                f"<button data-sample='{sample_id}' "
                f"onclick=\"showSample('{plot_id}', '{sample_id}')\">"
                f"{sample_id}</button>\n"
            )

        body += "</div>\n"

        for sample_id, html in samples.items():
            body += f"<div id='{plot_id}_{sample_id}' class='plot-wrapper'>\n{html}\n</div>\n"

        body += "</div>\n"

    with open(output_file, "w", encoding="utf-8") as out:
        out.write(HTML_TEMPLATE.format(css=CSS, body=body))

# Parse input arguments and call main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dirs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    main(args.input_dirs, args.output)
